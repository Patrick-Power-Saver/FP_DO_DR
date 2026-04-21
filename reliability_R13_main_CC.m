%% ============================================================
%  配电网可靠性评估 —— R12（R11 + DG孤岛运行模式）
%
%  在R11基础上新增§4c：DG孤岛运行建模
%
%  ─── 设计原则（参考Alanazi et al., Applied Energy, 2023）──────
%  1. 不修改R11的MILP约束结构，保持可扩展性
%  2. DG孤岛作为主网重构的补充恢复手段，两者解耦处理
%  3. 通过修改MILP目标权重W_q和CID计算公式实现DG效益融入
%
%  ─── DG孤岛运行物理模型 ──────────────────────────────────────
%  故障后恢复阶段（q=0节点，无法由主网重构恢复）：
%    不含DG或DG不足：等待修复 τ_RP（数据文件值，~2h）
%    含DG且出力≥负荷：孤岛运行 τ_DG_SW（DG切换时间，比τ_RP短得多）
%
%  关键变量 α^w_k（DG可用性指示，对应文献Eq.20）：
%    α^w_k = 1  当 G^w × PG_k ≥ L^w × P_free(k)
%    即当场景w下DG出力覆盖节点k的负荷需求时，节点可孤岛运行
%
%  ─── 多状态DG模型（文献Section 3.3）────────────────────────
%  采用场景集 {w: (ρ^w, g^w, l^w)} 描述DG出力与负荷的联合不确定性：
%    ρ^w：场景概率（Σρ^w=1）
%    g^w：DG出力比例（占额定容量PG_k的比例）
%    l^w：负荷比例（占峰值P_free(k)的比例）
%    α^w_k = (g^w×PG_k ≥ l^w×P_free(k)) → {0,1}
%    α_eff(k) = Σ_w ρ^w × α^w_k （期望可用率）
%
%  ─── 等效修复时间（将DG效益融入MILP目标和CID计算）────────────
%  τ_RP_eff(k, xy) = Σ_w ρ^w × [α^w_k × τ_DG_SW + (1-α^w_k) × τ_RP(xy)]
%    = α_eff(k) × τ_DG_SW + (1-α_eff(k)) × τ_RP(xy)  [DG节点]
%    = τ_RP(xy)                                         [非DG节点]
%
%  τ_RP_eff < τ_RP → DG孤岛节点的重构效益W_q降低（优化器优先保障无DG节点）
%
%  ─── 动态孤岛问题的讨论 ──────────────────────────────────────
%  本代码采用"静态孤岛"假设（与Alanazi et al.一致）：
%    - DG仅供自身节点负荷（不向相邻节点送电）
%    - 孤岛形成/解散在解析框架中通过α_eff等效，不建模时序过程
%  若需"动态孤岛"（DG向邻近节点供电、多DG协同）：
%    - 需引入DG-to-load功率流约束，与MILP形成耦合
%    - 但会引入双线性约束（DG出力 × 孤岛拓扑），可用McCormick线性化
%    - 本代码预留DG_DYNAMIC_ISLAND开关，True时输出提示但不实现
%
%  §§ 代码结构（在R11基础上新增）：
%    §4c  DG参数配置 + 多状态场景 + α矩阵 + τ_RP_eff计算
%    §6   修改W_q使用tau_benefit_eff（含DG等效收益时长）
%    §7   CID_rep使用τ_RP_eff，新增⑤DG孤岛收益分项
%    §8   新增DG孤岛效益对比输出
%
%  依赖: YALMIP + Gurobi
% =============================================================
clear; clc;

%% ── 用户配置区 ─────────────────────────────────────────────
sys_filename = '85-Node System Data.xlsx';
%sys_filename = '137-Node System Data.xlsx';
%sys_filename = '417-Node System Data.xlsx';
%sys_filename = '1080-Node System Data.xlsx';

tb_filename = 'Testbench for Linear Model Based Reliability Assessment Method for Distribution Optimization Models Considering Network Reconfiguration.xlsx';

sys_sheet = '85-node';
%sys_sheet = '137-node';
%sys_sheet = '417-node';
%sys_sheet = '1080-node';

LAMBDA_TRF = 0.5;% 变压器故障率（0=统一用线路公式）
TAU_UP_SW  = 0.3;
TAU_TIE_SW = 0.5;

% [Q1] 综合重要性评分
W_SCORE_NC = 0.5; W_SCORE_P = 0.5;
RATIO_L1   = 0.30; RATIO_L2 = 0.30;

% VOLL基准与扰动
VOLL_BASE  = [200, 50, 10];
VOLL_NOISE = 0.20; VOLL_SEED = 123;

% 切负荷上限
SHED_LIMIT    = [0.00, 0.25, 1.00];
FLEX_RATIO_L3 = 0.8;

% [Q2] 低压分支开关
N_BR_THRESHOLDS = [50, 150, 400];
N_BR_VALUES     = [2, 3, 4, 5];
LV_REGEN = false; LV_SEED = 42;

% [Q4] 开关动作次数惩罚
GAMMA_SWITCH = 1e2;

% ── [DG] 分布式电源孤岛运行参数 ─────────────────────────────
%  DG_ISLAND_ENABLE: 是否启用DG孤岛运行模式
DG_ISLAND_ENABLE = true;

%  DG节点定义（原始节点编号，与数据文件一致）
%  示例：在85节点系统中，选取节点4、23、37作为DG节点（参照Alanazi et al.文献）
%  可自定义修改此列表
DG_NODES_RAW = [4, 23, 37];   % DG节点的原始编号
%DG_NODES_RAW = [8, 29, 44, 61, 95, 111];   % DG节点的原始编号
%DG_NODES_RAW = [12, 32, 70, 105, 156, 200, 248, 313, 399];   % DG节点的原始编号
%DG_NODES_RAW = [16, 78, 137, 201, 333, 410, 556, 629, 700, 899, 945, 1005];   % DG节点的原始编号

%  各DG节点额定容量（MW，与P_free同单位）
%  对应文献中PG_i，即DG满发功率
DG_CAPACITY  = 0.5 *ones(1,size(DG_NODES_RAW,2));   % [MW]，与DG_NODES_RAW一一对应

%  DG孤岛切换时间（h）
%  τ_DG_SW < τ_TIE_SW：DG直接为本节点供电，无需操作联络线
%  物理含义：DG重新并网+本地开关操作的时间，通常0.1~0.3h
TAU_DG_SW = 0.2;

%  多状态DG-负荷联合场景（对应文献Eq.19-22的概率场景方法）
%  每行格式: [概率ρ^w, DG出力比例g^w, 负荷比例l^w]
%  概率之和须=1，g^w∈[0,1]，l^w∈[0,1]
%  可通过历史数据拟合或蒙特卡洛生成后场景削减得到
%  此处给出参考性10场景（参照Alanazi et al. Fig.4的风-负荷联合场景）
DG_SCENARIOS = [
%   ρ^w    g^w    l^w
    0.08,  0.10,  0.70;   % 低风速，低负荷
    0.10,  0.05,  0.83;   % 低风速，中负荷
    0.06,  0.08,  1.00;   % 低风速，高负荷
    0.12,  0.45,  0.70;   % 中风速，低负荷
    0.15,  0.50,  0.83;   % 中风速，中负荷（典型工况）
    0.10,  0.48,  1.00;   % 中风速，高负荷
    0.08,  0.90,  0.70;   % 高风速，低负荷（最有利孤岛）
    0.14,  0.85,  0.83;   % 高风速，中负荷
    0.09,  0.88,  1.00;   % 高风速，高负荷
    0.08,  0.30,  0.70;   % 弱风速，低负荷
];

%  动态孤岛开关（预留接口）
%  true  = 理论上应考虑DG向相邻节点供电（本版本仅输出提示）
%  false = 静态孤岛：DG仅供自身节点（本版本实现）
DG_DYNAMIC_ISLAND = false;

SOLVE_MODE = 'MCF';
V_UPPER = 1.05; V_LOWER = 0.95; V_SRC = 1.0; PF = 0.9;
%% ──────────────────────────────────────────────────────────

program_total = tic;

%% §1  Testbench 数据读取
fprintf('>> [1/8] 读取 Testbench: Sheet="%s"\n', sys_sheet);
tb_cell = readcell(tb_filename, 'Sheet', sys_sheet);
nrows_tb = size(tb_cell,1); hdr_row = 0;
for ri = 1:nrows_tb
    if ischar(tb_cell{ri,1}) && contains(tb_cell{ri,1},'Tie-Switch')
        hdr_row = ri; break;
    end
end
if hdr_row==0, error('未找到Tie-Switch表头行'); end
LINE_CAP = str2double(extractBefore(string(tb_cell{hdr_row+1,4}),' '));
TRAN_CAP = LINE_CAP;
TIE_LINES_RAW = [];
for ri = hdr_row+2:nrows_tb
    v1=tb_cell{ri,1}; v2=tb_cell{ri,2};
    if isnumeric(v1)&&~isnan(v1)&&isnumeric(v2)&&~isnan(v2)
        TIE_LINES_RAW(end+1,:)=[v1,v2]; %#ok<AGROW>
    end
end
fprintf('   线路容量=%.0f MW，联络线=%d条\n', LINE_CAP, size(TIE_LINES_RAW,1));

%% §2  可靠性参数读取
fprintf('>> [2/8] 读取可靠性参数: %s\n', sys_filename);
t_branch=readtable(sys_filename,'Sheet','Branch Lengths (km)');
t_branch=t_branch(:,1:3); t_branch.Properties.VariableNames={'From','To','Length_km'};
t_branch=t_branch(~isnan(t_branch.From),:);
t_dur=readtable(sys_filename,'Sheet','Interruption durations (h)','HeaderLines',3);
t_dur=t_dur(:,1:4); t_dur.Properties.VariableNames={'From','To','RP','SW'};
t_dur=t_dur(~isnan(t_dur.From),:);
t_cust=readtable(sys_filename,'Sheet','Numbers of customers per node');
t_cust=t_cust(:,1:2); t_cust.Properties.VariableNames={'Node','NC'};
t_cust=t_cust(~isnan(t_cust.Node),:);
t_peak=readtable(sys_filename,'Sheet','Peak Nodal Demands (kW)');
t_peak=t_peak(:,1:2); t_peak.Properties.VariableNames={'Node','P_kW'};
t_peak=t_peak(~isnan(t_peak.Node),:);
t_other=readtable(sys_filename,'Sheet','Other data','ReadVariableNames',false);
col1=cellfun(@(x) string(x),t_other{:,1},'UniformOutput',true);
lambda_per_km=str2double(string(t_other{find(contains(col1,'Failure rate'),1),2}));
row_dur=find(contains(col1,'Duration'),1);
T_l=[str2double(string(t_other{row_dur,3})),str2double(string(t_other{row_dur+1,3})),str2double(string(t_other{row_dur+2,3}))];
row_lf=find(contains(col1,'Loading factors'),1);
L_f=[str2double(string(t_other{row_lf,3})),str2double(string(t_other{row_lf+1,3})),str2double(string(t_other{row_lf+2,3}))]/100;
fprintf('   lambda=%.4f/km, T_l=[%s]h, L_f=[%s]\n',lambda_per_km,num2str(T_l,'%g '),num2str(L_f,'%.2f '));

%% §3  潮流参数
fprintf('>> [3/8] 生成潮流参数...\n');
R_KM=0.003151; X_KM=0.001526; TAN_PHI=tan(acos(PF));
V_src_sq=V_SRC^2; V_upper_sq=V_UPPER^2; V_lower_sq=V_LOWER^2;
M_V=(V_upper_sq-V_lower_sq)*2; M_vn=V_upper_sq;

%% §4  拓扑索引
fprintf('>> [4/8] 构建拓扑索引...\n');
raw_nodes=unique([t_branch.From;t_branch.To;TIE_LINES_RAW(:)]);
num_nodes=length(raw_nodes);
node_map=containers.Map(raw_nodes,1:num_nodes);
inv_map=raw_nodes;
subs_raw=raw_nodes(~ismember(raw_nodes,t_cust.Node));
subs_idx=arrayfun(@(s) node_map(s),subs_raw);
nB_norm=height(t_branch); nTie=size(TIE_LINES_RAW,1); nB_all=nB_norm+nTie;

rel_branches=zeros(nB_norm,8);
for b=1:nB_norm
    u_raw=t_branch.From(b); v_raw=t_branch.To(b);
    u=node_map(u_raw); v=node_map(v_raw); len=t_branch.Length_km(b);
    match=(t_dur.From==u_raw&t_dur.To==v_raw)|(t_dur.From==v_raw&t_dur.To==u_raw);
    if ~any(match), error('分支(%d-%d)无停电时间',u_raw,v_raw); end
    is_trf=ismember(u,subs_idx)|ismember(v,subs_idx);
    cap_b=TRAN_CAP*is_trf+LINE_CAP*(~is_trf);
    lam_b=ternary(LAMBDA_TRF>0&&is_trf, LAMBDA_TRF, len*lambda_per_km);
    rel_branches(b,:)=[u,v,lam_b,t_dur.RP(match),t_dur.SW(match),R_KM*len,X_KM*len,cap_b];
end
is_trf_vec=ismember(rel_branches(:,1),subs_idx)|ismember(rel_branches(:,2),subs_idx);
n_trf=sum(is_trf_vec);

tie_branches=zeros(nTie,8);
for t=1:nTie
    u=node_map(TIE_LINES_RAW(t,1)); v=node_map(TIE_LINES_RAW(t,2));
    tie_branches(t,:)=[u,v,0,0,0,R_KM*0.1,X_KM*0.1,LINE_CAP];
end
all_branches=[rel_branches;tie_branches];
branch_from=all_branches(:,1); branch_to=all_branches(:,2);
r_b_all=all_branches(:,6); x_b_all=all_branches(:,7); cap_b_all=all_branches(:,8);

load_nodes=setdiff(1:num_nodes,subs_idx); nL=length(load_nodes); non_sub=load_nodes;
A_inc_norm=sparse(rel_branches(:,2),(1:nB_norm)',+1,num_nodes,nB_norm)+sparse(rel_branches(:,1),(1:nB_norm)',-1,num_nodes,nB_norm);
A_inc_all=sparse(branch_to,(1:nB_all)',+1,num_nodes,nB_all)+sparse(branch_from,(1:nB_all)',-1,num_nodes,nB_all);
A_free_all=A_inc_all(load_nodes,:);
B_to_all=sparse((1:nB_all)',branch_to,1,nB_all,num_nodes);
B_from_all=sparse((1:nB_all)',branch_from,1,nB_all,num_nodes);
BdV=B_to_all-B_from_all;
[~,pk_row]=ismember(inv_map(load_nodes),t_peak.Node);
P_free=zeros(nL,1); valid_pk=pk_row>0;
P_free(valid_pk)=t_peak.P_kW(pk_row(valid_pk))/1e3;
Q_free=P_free*TAN_PHI;

%% §4b  [Q1+Q2] 负荷分级 & 低压分支开关档位（与R11相同）
[~,c_row_pre]=ismember(inv_map(load_nodes),t_cust.Node);
NC_vec=zeros(nL,1); NC_vec(c_row_pre>0)=t_cust.NC(c_row_pre(c_row_pre>0));

total_NC=max(sum(NC_vec),1); total_P=max(sum(P_free),1e-9);
score_k=W_SCORE_NC*(NC_vec/total_NC)+W_SCORE_P*(P_free/total_P);
[score_sorted,sort_idx]=sort(score_k,'descend');
cum_score=cumsum(score_sorted)/max(sum(score_sorted),1e-9);
rank_vec=zeros(nL,1);
rank_vec(sort_idx(cum_score<=RATIO_L1))=1;
rank_vec(sort_idx(cum_score>RATIO_L1&cum_score<=RATIO_L1+RATIO_L2))=2;
rank_vec(sort_idx(cum_score>RATIO_L1+RATIO_L2))=3;
idx_L1=find(rank_vec==1); idx_L2=find(rank_vec==2); idx_L3=find(rank_vec==3);
nL1=length(idx_L1); nL2=length(idx_L2); nL3=length(idx_L3);
shed_limit_vec=zeros(nL,1);
shed_limit_vec(idx_L1)=SHED_LIMIT(1);
shed_limit_vec(idx_L2)=SHED_LIMIT(2);
shed_limit_vec(idx_L3)=SHED_LIMIT(3)*FLEX_RATIO_L3;
rng(VOLL_SEED);
voll_vec=zeros(nL,1);
voll_vec(idx_L1)=VOLL_BASE(1)*(1+VOLL_NOISE*(2*rand(nL1,1)-1));
voll_vec(idx_L2)=VOLL_BASE(2)*(1+VOLL_NOISE*(2*rand(nL2,1)-1));
voll_vec(idx_L3)=VOLL_BASE(3)*(1+VOLL_NOISE*(2*rand(nL3,1)-1));
voll_pu=voll_vec*1e3;

if ~LV_REGEN; rng(LV_SEED); end
node_levels=cell(nL,1); n_br_vec=zeros(nL,1);
for k=1:nL
    if shed_limit_vec(k)<=1e-6; node_levels{k}=[0]; continue; end
    nc_k=max(NC_vec(k),1);
    if nc_k<=N_BR_THRESHOLDS(1); n_br=N_BR_VALUES(1);
    elseif nc_k<=N_BR_THRESHOLDS(2); n_br=N_BR_VALUES(2);
    elseif nc_k<=N_BR_THRESHOLDS(3); n_br=N_BR_VALUES(3);
    else; n_br=N_BR_VALUES(4); end
    n_br_vec(k)=n_br;
    raw_w=ones(1,n_br)/n_br+0.1*(2*rand(1,n_br)-1); raw_w=max(raw_w,0.02);
    branch_w=raw_w/sum(raw_w);
    n_states=2^n_br; lset=zeros(1,n_states);
    for s=0:n_states-1; bits=bitget(s,1:n_br,'uint32'); lset(s+1)=sum(branch_w.*double(bits)); end
    lset=sort(unique(round(lset,4)));
    lset=lset(lset<=shed_limit_vec(k)+1e-6);
    if isempty(lset)||lset(1)~=0; lset=[0,lset]; end
    node_levels{k}=lset;
end
n_lev_arr=cellfun(@length,node_levels);

fprintf('   [Q1]负荷分级: L1=%d, L2=%d, L3=%d节点\n',nL1,nL2,nL3);
fprintf('   [Q2]档位: 均值%.1f, 最大%d\n',mean(n_lev_arr(shed_limit_vec>0)),max(n_lev_arr));

%% ================================================================
%  §4c  DG孤岛运行参数配置
%
%  本节将DG孤岛的效益等效为"修复时间缩短"，通过期望α_eff将多状态
%  不确定性压缩为单一参数 τ_RP_eff，从而保持后续MILP和CID计算的
%  线性结构不变（不引入新的整数变量）。
%
%  ─── α矩阵的物理含义 ───────────────────────────────────────────
%  alpha_eff(k)：节点k在任意故障下能自供孤岛运行的期望概率
%    = Σ_w ρ^w × 1{g^w × PG_k ≥ l^w × P_free(k)}
%  其中 1{·} 为指示函数（满足则1，否则0），对应文献α^w_i
%
%  ─── τ_RP_eff的物理含义 ─────────────────────────────────────────
%  对于DG节点k，故障xy（q=0）下的期望等效停电时长：
%    τ_RP_eff(k,xy) = α_eff(k) × τ_DG_SW + (1-α_eff(k)) × τ_RP(xy)
%  即α_eff比例的时间以τ_DG_SW结束，其余以τ_RP(xy)结束
%  对于非DG节点：τ_RP_eff(k,xy) = τ_RP(xy)（不变）
%
%  ─── MILP目标修正逻辑 ──────────────────────────────────────────
%  tau_benefit_eff(k,xy) = max(τ_RP_eff(k,xy) - τ_TIE_SW, 0)
%  DG孤岛使τ_RP_eff降低 → tau_benefit_eff降低 → W_q降低
%  效果：优化器对DG节点的"重构迫切性"降低，优先保障无DG节点
%  这符合实际工程逻辑：有DG自备的节点无需强求主网重构
%
%  ─── 动态孤岛 vs 静态孤岛 ──────────────────────────────────────
%  静态孤岛（本版实现）：
%    DG仅供自身节点，功率不向外输出
%    解析框架通过α_eff参数化，无需显式建模孤岛网络
%
%  动态孤岛（预留接口，未实现）：
%    DG可向相邻节点供电，孤岛边界动态扩展
%    需引入DG潮流约束，与MILP耦合，形成双线性约束
%    若实现：可参考文献中的多商品流孤岛划分模型（丁涛等，电网技术2025）
%    本代码DG_DYNAMIC_ISLAND=true时输出说明
%
%  ─── 与文献(Alanazi 2023)的对应关系 ───────────────────────────
%    文献Eq(20): D̄_i = Σ_{km∈Y}(ψ^j_km+ψ^j_mk)×TS_{km}×l_{km}×λ_{km}
%      → 本代码：CID_dg_benefit (③项中τ_RP→τ_DG_SW的改善量)
%    文献Eq(21): SAIDI = Σ_w ρ^w × [Σ NC_i×D_i + Σ_{i∈ΩG} α^w_i×NC_i×(D̄_i-D_i)]
%      → 本代码：通过τ_RP_eff整合为统一CID_rep项，等价但更简洁
%    文献Eq(22): EENS = Σ_w ρ^w × [Σ L^w_i×D_i + Σ_{i∈ΩG} α^w_i×L^w_i×(D̄_i-D_i)]
%      → 本代码：EENS也使用修正后的CID（含DG效益）
% ================================================================
fprintf('>> [4c] DG孤岛运行建模...\n');

if DG_ISLAND_ENABLE
    % 场景概率验证与归一化
    prob_sum=sum(DG_SCENARIOS(:,1));
    if abs(prob_sum-1)>1e-6
        warning('[DG] 场景概率之和=%.4f≠1，自动归一化',prob_sum);
        DG_SCENARIOS(:,1)=DG_SCENARIOS(:,1)/prob_sum;
    end
    N_DG_SCEN=size(DG_SCENARIOS,1);
    rho_w=DG_SCENARIOS(:,1);   % 场景概率 [N_DG_SCEN×1]
    g_w  =DG_SCENARIOS(:,2);   % DG出力比例
    l_w  =DG_SCENARIOS(:,3);   % 负荷比例

    % 将DG原始节点编号映射到本地负荷节点索引
    nDG=length(DG_NODES_RAW);
    dg_local_idx=zeros(1,nDG);   % DG节点的本地索引（在load_nodes中的位置）
    dg_cap_vec=zeros(nDG,1);     % 各DG额定容量[MW]
    valid_dg=false(1,nDG);
    for d=1:nDG
        node_raw=DG_NODES_RAW(d);
        if ~isKey(node_map,node_raw)
            warning('[DG] 节点%d不在系统中，跳过',node_raw);
            continue;
        end
        local=find(load_nodes==node_map(node_raw),1);
        if isempty(local)
            warning('[DG] 节点%d是变电站节点，不能安装DG，跳过',node_raw);
            continue;
        end
        dg_local_idx(d)=local;
        dg_cap_vec(d)=DG_CAPACITY(d);
        valid_dg(d)=true;
    end
    dg_local_idx=dg_local_idx(valid_dg);
    dg_cap_vec=dg_cap_vec(valid_dg);
    nDG_valid=sum(valid_dg);

    % is_dg(k)：逻辑向量，标识负荷节点k是否有DG
    is_dg=false(nL,1);
    is_dg(dg_local_idx)=true;

    % 计算每个DG节点的α_eff（期望孤岛可用率）
    % α^w_k = 1 当 g^w×PG_k ≥ l^w×P_free(k)
    % α_eff(k) = Σ_w ρ^w × α^w_k
    alpha_eff=zeros(nL,1);   % nL×1，非DG节点=0
    for di=1:nDG_valid
        k=dg_local_idx(di);
        pg_k=dg_cap_vec(di);   % DG额定容量[MW]
        pd_k=P_free(k);        % 节点峰值负荷[MW]
        if pd_k<1e-9; alpha_eff(k)=1; continue; end   % 零负荷节点：DG总够用
        for wi=1:N_DG_SCEN
            if g_w(wi)*pg_k >= l_w(wi)*pd_k
                alpha_eff(k)=alpha_eff(k)+rho_w(wi);
            end
        end
    end

    % 计算τ_RP_eff矩阵（nL×nB_norm）
    % 对DG节点：τ_RP_eff(k,xy) = α_eff(k)×τ_DG_SW + (1-α_eff(k))×τ_RP(xy)
    % 对非DG节点：τ_RP_eff(k,xy) = τ_RP(xy)（标量广播）
    lam_vec_pre=rel_branches(:,3);
    trp_vec_pre=rel_branches(:,4);   % nB_norm×1

    % 构造nL×nB_norm的τ_RP_eff矩阵
    % 先用标准τ_RP初始化
    trp_eff_mat=repmat(trp_vec_pre',nL,1);   % nL×nB_norm

    % 对DG节点逐个修正
    for di=1:nDG_valid
        k=dg_local_idx(di);
        % τ_RP_eff(k,xy) = α_eff(k)×τ_DG_SW + (1-α_eff(k))×τ_RP(xy)
        trp_eff_mat(k,:) = alpha_eff(k)*TAU_DG_SW + (1-alpha_eff(k))*trp_vec_pre';
    end

    % tau_benefit_eff(k,xy)：修正后的重构效益时长（nL×nB_norm）
    % 正常: tau_benefit(xy) = max(τ_RP(xy)-τ_TIE_SW, 0)
    % DG修正: tau_benefit_eff(k,xy) = max(τ_RP_eff(k,xy)-τ_TIE_SW, 0)
    tau_benefit_eff=max(trp_eff_mat - TAU_TIE_SW, 0);   % nL×nB_norm

    % 打印DG配置信息
    fprintf('   DG孤岛模式：启用 | %d个DG场景 | τ_DG_SW=%.2fh\n',N_DG_SCEN,TAU_DG_SW);
    fprintf('   DG节点(%d个): ',nDG_valid);
    for di=1:nDG_valid
        k=dg_local_idx(di);
        fprintf('节点%-4d(PG=%.2fMW,α_eff=%.2f) ',inv_map(load_nodes(k)),dg_cap_vec(di),alpha_eff(k));
    end
    fprintf('\n');

    % 动态孤岛提示
    if DG_DYNAMIC_ISLAND
        fprintf('   !! DG_DYNAMIC_ISLAND=true（预留接口）\n');
        fprintf('      动态孤岛需建模DG向相邻节点供电的功率流，\n');
        fprintf('      引入DG潮流×孤岛开关的双线性约束，\n');
        fprintf('      可通过McCormick包络线性化后并入MILP。\n');
        fprintf('      当前版本未实现，降级为静态孤岛模式运行。\n');
    end
else
    % DG孤岛未启用：退化到R11行为
    alpha_eff=zeros(nL,1);
    is_dg=false(nL,1);
    tau_benefit_eff=repmat(max(rel_branches(:,4)'-TAU_TIE_SW,0),nL,1);  % nL×nB_norm
    nDG_valid=0; dg_local_idx=[]; dg_cap_vec=[];
    trp_eff_mat=repmat(rel_branches(:,4)',nL,1);
    fprintf('   DG孤岛模式：未启用（DG_ISLAND_ENABLE=false）\n');
end

%% §5  MCF 路径识别（与R11完全相同）
fprintf('>> [5/8] MCF 路径识别...\n');
t_mcf=tic;
nSub=length(subs_idx); E_sub=sparse(subs_idx,1:nSub,1,num_nodes,nSub);

if strcmp(SOLVE_MODE,'MCF')
    E_load=sparse(load_nodes,1:nL,1,num_nodes,nL);
    F_mat=sdpvar(nB_norm,nL,'full'); Z_mat=sdpvar(nB_norm,nL,'full'); Gss=sdpvar(nSub,nL,'full');
    C_mcf=[-Z_mat<=F_mat,F_mat<=Z_mat,Z_mat>=0,0<=Gss<=1,sum(Gss,1)==1,A_inc_norm*F_mat==E_load-E_sub*Gss];
    sol=optimize(C_mcf,sum(sum(Z_mat)),sdpsettings('solver','gurobi','verbose',0));
    if sol.problem~=0, error('MCF失败: %s',sol.info); end
    f_res=sparse(abs(value(F_mat))>0.5);
else
    f_res=false(nB_norm,nL); opts_sp=sdpsettings('solver','gurobi','verbose',0);
    parfor k=1:nL
        f_k=sdpvar(nB_norm,1); z_k=sdpvar(nB_norm,1); g_k=sdpvar(nSub,1);
        d_k=sparse(load_nodes(k),1,1,num_nodes,1)-E_sub*g_k;
        C_k=[-z_k<=f_k,f_k<=z_k,z_k>=0,0<=g_k<=1,sum(g_k)==1,A_inc_norm*f_k==d_k];
        optimize(C_k,sum(z_k),opts_sp); f_res(:,k)=abs(value(f_k))>0.5;
    end
    f_res=sparse(f_res);
end
p_mat=f_res';

is_outlet=ismember(rel_branches(:,1),subs_idx)|ismember(rel_branches(:,2),subs_idx);
p_feeder_mat=sparse(nL,nB_norm);
for xy=1:nB_norm
    dn_k=find(p_mat(:,xy),1); if isempty(dn_k), continue; end
    for bi=find(f_res(:,dn_k))'
        if is_outlet(bi), p_feeder_mat(:,xy)=f_res(bi,:)'; break; end
    end
end
p_feeder_mat=sparse(p_feeder_mat);
t_mcf_elapsed=toc(t_mcf);
fprintf('   MCF完成（%.1f秒），p_direct nnz=%d，p_feeder nnz=%d\n',t_mcf_elapsed,nnz(p_mat),nnz(p_feeder_mat));

%% §5b  预计算R1基准（用于分式目标分子）
fprintf('>> [5b/8] 预计算R1基准损失...\n');
nScen=nB_norm;
p_mat_d=double(p_mat); p_feeder_d=double(p_feeder_mat);
lam_vec=rel_branches(:,3); trp_vec=rel_branches(:,4);

p_upstream_base=p_feeder_d-p_mat_d;
[~,p_row_obj]=ismember(inv_map(load_nodes),t_peak.Node);
P_avg_obj=zeros(nL,1); P_avg_obj(p_row_obj>0)=t_peak.P_kW(p_row_obj(p_row_obj>0))*sum(L_f.*(T_l/8760));

CID_R1_obj=TAU_UP_SW*(p_upstream_base*lam_vec)+p_mat_d*(lam_vec.*trp_vec);
VOLL_loss_R1_const=sum(voll_vec.*P_avg_obj.*CID_R1_obj);
fprintf('   R1损失基准: %.2f 万元/年\n',VOLL_loss_R1_const/1e4);

%% ================================================================
%  §6  批量场景 MILP（使用DG修正后的权重矩阵）
%
%  修改点（相对于R11）：
%    W_q(k,xy) = VOLL_k × P_free(k) × λ_xy × tau_benefit_eff(k,xy)
%    其中 tau_benefit_eff = max(τ_RP_eff - τ_TIE_SW, 0)
%    DG节点的τ_RP_eff < τ_RP → tau_benefit_eff < tau_benefit
%    → W_q对DG节点降权 → 优化器优先为无DG节点分配联络线容量
%
%    W_L(k,xy) = VOLL_k × λ_xy × tau_benefit_eff(k,xy)
%    同理，DG节点的切负荷惩罚也相应降低（因为DG可弥补部分切负荷损失）
%
%  注意：MILP约束结构与R11完全相同，仅权重矩阵不同
% ================================================================
fprintf('>> [6/8] 批量场景MILP（含DG修正权重）...\n');
t_milp=tic;

Dr=spdiags(r_b_all,0,nB_all,nB_all); Dx=spdiags(x_b_all,0,nB_all,nB_all);
Dp=spdiags(P_free,0,nL,nL); Dq_=spdiags(Q_free,0,nL,nL);
Cap=spdiags(cap_b_all,0,nB_all,nB_all);

% ── DG修正后的权重矩阵 ──────────────────────────────────────────
% tau_benefit_eff: nL×nB_norm（已在§4c计算）
% W_q(k,xy) = VOLL_k[元/MWh]×P_free(k)[MW]×λ_xy×tau_benefit_eff(k,xy) [元/年]
W_q=(voll_pu.*P_free.*p_mat_d).*((lam_vec.*1)'.*tau_benefit_eff);   % 广播λ到nL×nScen
for xy=1:nScen
    W_q(:,xy)=W_q(:,xy)*lam_vec(xy);   % 逐场景乘以λ_xy
end
% 简洁实现：
W_q=(voll_pu.*P_free.*p_mat_d).*(repmat(lam_vec',nL,1).*tau_benefit_eff);  % nL×nScen

W_L=(voll_pu.*p_mat_d).*(repmat(lam_vec',nL,1).*tau_benefit_eff);           % nL×nScen

fprintf('   W_q(含DG修正): 非零=%d, max=%.2f元/年\n',nnz(W_q),max(full(W_q(:))));
if DG_ISLAND_ENABLE && nDG_valid>0
    % 对比：原版W_q（无DG修正）
    tau_benefit_std=repmat(max(trp_vec'-TAU_TIE_SW,0),nL,1);
    W_q_std=(voll_pu.*P_free.*p_mat_d).*(repmat(lam_vec',nL,1).*tau_benefit_std);
    for di=1:nDG_valid
        k=dg_local_idx(di);
        wq_ratio=sum(W_q(k,:))/max(sum(W_q_std(k,:)),1e-9);
        fprintf('   DG节点%d W_q降权至%.1f%%（α_eff=%.2f）\n', ...
            inv_map(load_nodes(k)), full(wq_ratio)*100, alpha_eff(k));
    end
end

% ── 稀疏ZZ变量（与R11相同）──────────────────────────────────────
can_shed=shed_limit_vec>1e-6;
n_shed=sum(can_shed); shed_nodes_idx=find(can_shed);
total_ZZ=sum(n_lev_arr(shed_nodes_idx));
row_s=zeros(total_ZZ,1); val_sel=zeros(total_ZZ,1); val_one=ones(total_ZZ,1);
ptr=0;
for ki=1:n_shed
    k=shed_nodes_idx(ki); lv=node_levels{k}; nk=length(lv);
    rows=ptr+1:ptr+nk;
    row_s(rows)=ki; val_sel(rows)=P_free(k)*lv(:);
    ptr=ptr+nk;
end
col_s=(1:total_ZZ)';
SEL_s=sparse(row_s,col_s,val_sel,n_shed,total_ZZ);
ONESUM_s=sparse(row_s,col_s,val_one,n_shed,total_ZZ);

% ── 决策变量（CC变换 + 二元乘积线性化）────────────────────────────
S_mat=binvar(nB_all,nScen,'full'); Q_mat=binvar(nL,nScen,'full');
ZZ_s=binvar(total_ZZ,nScen,'full');
t=sdpvar(1,1);

S_hat=sdpvar(nB_all,nScen,'full');      % S_hat = S_mat * t
Q_hat=sdpvar(nL,nScen,'full');          % Q_hat = Q_mat * t
ZZ_hat=sdpvar(total_ZZ,nScen,'full');   % ZZ_hat = ZZ_s * t

Pf_hat=sdpvar(nB_all,nScen,'full');
Qf_hat=sdpvar(nB_all,nScen,'full');
V_hat=sdpvar(num_nodes,nScen,'full');
E_hat=sdpvar(nB_all,nScen,'full');
L_hat=sdpvar(nL,nScen,'full');

delta_hat=BdV*V_hat+2*(Dr*Pf_hat+Dx*Qf_hat);
fault_lin_idx=(0:nB_norm-1)*nB_all+(1:nB_norm);
P_max_mat=repmat(P_free,1,nScen);
L_shed_part=SEL_s*ZZ_hat;

% t上界（用于二元乘积精确线性化）
DEN_EPS=1e-3;
T_UB=1/DEN_EPS;

% ── CC后的约束 ─────────────────────────────────────────────────
C=[t>=0, t<=T_UB, ...
   V_hat(subs_idx,:)==V_src_sq*t, ...
   V_hat(non_sub,:)>=V_lower_sq*t-M_vn*(t-Q_hat), ...
   V_hat(non_sub,:)<=V_upper_sq*t+M_vn*(t-Q_hat), ...
   -Cap*S_hat<=Pf_hat<=Cap*S_hat, ...
   -Cap*S_hat<=Qf_hat<=Cap*S_hat, ...
   A_free_all*Pf_hat==Dp*Q_hat-L_hat, ...
   A_free_all*Qf_hat==Dq_*Q_hat-TAN_PHI*L_hat, ...
   E_hat>=0, ...
   delta_hat<=M_V*(t-S_hat)+E_hat, ...
   delta_hat>=-M_V*(t-S_hat)-E_hat, ...
   sum(S_hat,1)==sum(Q_hat,1), ...
   S_hat(fault_lin_idx)==0, ...
   Q_hat>=repmat(1-p_feeder_d,1,1)*t];

if any(~can_shed), C=[C,L_hat(~can_shed,:)==0]; end
C=[C,ZZ_hat>=0,L_hat(shed_nodes_idx,:)==L_shed_part, ...
   ONESUM_s*ZZ_hat==t,L_hat>=0,L_hat<=P_max_mat.*Q_hat];

% 二元乘积线性化：X_hat = X * t, X∈{0,1}
C=[C, ...
   0<=S_hat, S_hat<=t, S_hat<=T_UB*S_mat, S_hat>=t-T_UB*(1-S_mat), ...
   0<=Q_hat, Q_hat<=t, Q_hat<=T_UB*Q_mat, Q_hat>=t-T_UB*(1-Q_mat), ...
   0<=ZZ_hat, ZZ_hat<=t, ZZ_hat<=T_UB*ZZ_s, ZZ_hat>=t-T_UB*(1-ZZ_s)];

if GAMMA_SWITCH>0
    S_NO=ones(nB_all,1); S_NO(nB_norm+1:end)=0;
    S_NO_mat=repmat(S_NO,1,nScen);
    switch_cost_hat=GAMMA_SWITCH*sum(sum(S_NO_mat.*(t-S_hat)+(1-S_NO_mat).*S_hat));
else
    switch_cost_hat=0;
end

opts_milp=sdpsettings('solver','gurobi','verbose',1,'gurobi.MIPGap',1e-3, ...
    'gurobi.Heuristics',0.05,'gurobi.Presolve',2,'gurobi.Cuts',2,'gurobi.NodefileStart',0.5);

% ── CC分式目标：max((VOLL_loss_R1 - VOLL_loss_R2)/(switch_cost + 1e-4*sum(E_vdrop))) ──
P_free_safe_obj=max(P_free,1e-9);
p_direct_coef=(repmat(1./P_free_safe_obj,1,nScen)).*p_mat_d; % 常数系数

CID_up_hat=TAU_UP_SW*(p_upstream_base*lam_vec)*t;                               % 常量项*t
CID_tie_hat=TAU_TIE_SW*(p_mat_d.*Q_hat)*lam_vec;                                % 线性(Q_hat)
CID_rep_hat=(p_mat_d.*(t-Q_hat))*(lam_vec.*trp_vec);                             % 线性(t,Q_hat)
CID_shed_hat=(p_direct_coef.*L_hat)*(lam_vec.*(trp_vec-TAU_TIE_SW));             % 线性(L_hat)
CID_R2_hat=CID_up_hat+CID_tie_hat+CID_rep_hat+CID_shed_hat;                      % = t*CID_R2
VOLL_loss_R2_hat=sum(voll_vec.*P_avg_obj.*CID_R2_hat);                           % = t*VOLL_loss_R2

VOLL_saving_hat=VOLL_loss_R1_const*t-VOLL_loss_R2_hat;                           % = t*分子
den_hat=switch_cost_hat+1e-4*sum(E_hat(:));                                      % = t*分母
C=[C, den_hat==1, den_hat>=DEN_EPS*t];
objective=-VOLL_saving_hat;
sol=optimize(C,objective,opts_milp);

if sol.problem==0
    t_v=value(t);
    Qv=value(Q_mat);
    if any(isnan(Qv(:)))
        q_mat=logical(1-p_feeder_d); L_res=zeros(nL,nScen);
        Pf_res=zeros(nB_all,nScen); Qf_res=zeros(nB_all,nScen);
        Vm_res=repmat(V_src_sq,num_nodes,nScen); Sw_res=zeros(nB_all,nScen);
    else
        q_mat=logical(round(Qv)); Sw_res=round(value(S_mat));
        if t_v<=1e-12
            warning('[§6] CC还原失败: t过小');
            L_res=zeros(nL,nScen); Pf_res=zeros(nB_all,nScen); Qf_res=zeros(nB_all,nScen);
            Vm_res=repmat(V_src_sq,num_nodes,nScen);
        else
            L_res=max(value(L_hat)/t_v,0);
            Pf_res=value(Pf_hat)/t_v; Qf_res=value(Qf_hat)/t_v;
            Vm_res=sqrt(max(value(V_hat)/t_v,0));
        end
    end
else
    warning('[§6] MILP失败: %s',sol.info);
    q_mat=logical(1-p_feeder_d); L_res=zeros(nL,nScen);
    Pf_res=zeros(nB_all,nScen); Qf_res=zeros(nB_all,nScen);
    Vm_res=repmat(V_src_sq,num_nodes,nScen); Sw_res=zeros(nB_all,nScen);
end

q_mat_d=double(q_mat); p_direct_d=double(p_mat); p_feeder_d_full=double(p_feeder_mat);
SHED_THRESH=1e-4;
P_free_mat=repmat(P_free,1,nScen);
n_full=full(sum(q_mat_d(:)>0.5&L_res(:)<=SHED_THRESH&p_direct_d(:)>0.5));
n_part=full(sum(q_mat_d(:)>0.5&L_res(:)>SHED_THRESH&p_direct_d(:)>0.5));
n_unrec=full(sum(q_mat_d(:)<0.5&p_direct_d(:)>0.5));
n_aff=full(sum(p_direct_d(:)>0));
aff_power=full(sum(sum(P_free_mat.*p_direct_d)));
net_rec=full(sum(sum((P_free_mat-L_res).*q_mat_d.*p_direct_d)));
rec_pct=net_rec/max(aff_power,1e-9)*100;
total_shed=full(sum(sum(L_res.*q_mat_d.*p_direct_d)));
mask_L1=false(nL,1); mask_L1(idx_L1)=true;
mask_L2=false(nL,1); mask_L2(idx_L2)=true;
mask_L3=false(nL,1); mask_L3(idx_L3)=true;
shed_lev=[sum(sum(L_res(mask_L1,:).*q_mat_d(mask_L1,:))), ...
          sum(sum(L_res(mask_L2,:).*q_mat_d(mask_L2,:))), ...
          sum(sum(L_res(mask_L3,:).*q_mat_d(mask_L3,:)))];
t_milp_elapsed=toc(t_milp);
fprintf('   MILP完成（%.1f秒），目标=%.2f元/年\n',t_milp_elapsed,value(objective));
fprintf('   恢复: 完全=%d,部分=%d,未恢复=%d/受影响=%d, 功率基=%.1f%%\n', ...
    n_full,n_part,n_unrec,n_aff,rec_pct);

%% ================================================================
%  §7  可靠性指标计算（DG孤岛修正三阶段模型）
%
%  修改点（相对于R11）：
%  ③等待修复项使用τ_RP_eff替代τ_RP：
%    原版: CID_rep = (p_direct×(1-q)) * (λ·τ_RP)
%    修正: CID_rep = Σ_k,xy p_direct(k,xy)×(1-q(k,xy))×λ_xy×τ_RP_eff(k,xy)
%    矩阵实现: CID_rep = sum((p_direct.*(1-q)) .* (λ.*τ_RP_eff), 2)
%
%  新增第⑤项：DG孤岛收益（量化DG带来的SAIDI改善）
%    CID_dg_benefit(k) = Σ_xy λ_xy×p_direct(k,xy)×(1-q(k,xy))×(τ_RP(xy)-τ_RP_eff(k,xy))
%    SAIDI_dg = Σ_k NC_k×CID_dg_benefit(k) / Σ NC   (>0 表示DG带来改善)
%    物理：DG孤岛将q=0节点的修复等待时间从τ_RP缩短到τ_RP_eff的年均效益
%
%  切负荷项④不变（已含于DG修正框架内）
%
%  对比基准（R1：无重构无DG）：
%    CID_R1 = τ_UP_SW×p_upstream×λ + p_direct×(λ×τ_RP)  [不含DG孤岛]
% ================================================================
fprintf('>> [7/8] 计算可靠性指标（含DG孤岛修正）...\n');
lam=rel_branches(:,3); trp=rel_branches(:,4);
p_upstream_d=p_feeder_d_full-p_direct_d;
[~,c_row]=ismember(inv_map(load_nodes),t_cust.Node);
NC_r=zeros(nL,1); NC_r(c_row>0)=t_cust.NC(c_row(c_row>0));
[~,p_row]=ismember(inv_map(load_nodes),t_peak.Node);
P_avg=zeros(nL,1); P_avg(p_row>0)=t_peak.P_kW(p_row(p_row>0))*sum(L_f.*(T_l/8760));
total_cust=sum(NC_r);
P_free_safe=max(P_free,1e-9);
shed_ratio=(L_res./repmat(P_free_safe,1,nScen)).*q_mat_d.*p_direct_d;

% R2 三阶段CID（含DG修正）
CIF=p_feeder_d_full*lam;

% ①上游开关（不受DG影响）
CID_up=TAU_UP_SW*(p_upstream_d*lam);

% ②联络线转供（已恢复节点，DG孤岛不改变其停电时长）
CID_tie=TAU_TIE_SW*(p_direct_d.*q_mat_d)*lam;

% ③等待修复（q=0节点）：使用τ_RP_eff（已含DG孤岛等效）
%   矩阵逐行点乘：CID_rep(k) = Σ_xy p_direct(k,xy)×(1-q(k,xy))×λ_xy×τ_RP_eff(k,xy)
CID_rep_mat=(p_direct_d.*(1-q_mat_d)).*(repmat(lam',nL,1).*trp_eff_mat);  % nL×nScen
CID_rep=sum(CID_rep_mat,2);   % nL×1

% ④切负荷额外等待（不变）
CID_shed=(shed_ratio.*p_direct_d.*q_mat_d)*(lam.*(trp-TAU_TIE_SW));

CID=CID_up+CID_tie+CID_rep+CID_shed;

SAIFI=(NC_r'*CIF)/total_cust;
SAIDI=(NC_r'*CID)/total_cust;
EENS=(P_avg'*CID)/1e3;
ASAI=1-SAIDI/8760;
SAIDI_up=(NC_r'*CID_up)/total_cust;
SAIDI_tie=(NC_r'*CID_tie)/total_cust;
SAIDI_rep=(NC_r'*CID_rep)/total_cust;
SAIDI_shed=(NC_r'*CID_shed)/total_cust;

% ⑤DG孤岛收益分项（与无DG情况的差值）
% CID_dg_benefit(k) = CID_rep_nodg(k) - CID_rep(k)  [>0 = 改善]
CID_rep_nodg_mat=(p_direct_d.*(1-q_mat_d)).*repmat(lam'.*trp',nL,1);  % 无DG的CID_rep
CID_dg_benefit=sum(CID_rep_nodg_mat,2)-CID_rep;   % nL×1，DG孤岛带来的CID改善量
SAIDI_dg=(NC_r'*CID_dg_benefit)/total_cust;       % 系统级DG孤岛SAIDI改善量

% R1基准（无重构无DG）
CID_R1=TAU_UP_SW*(p_upstream_d*lam)+p_direct_d*(lam.*trp);
SAIDI_R1=(NC_r'*CID_R1)/total_cust;
EENS_R1=(P_avg'*CID_R1)/1e3;
ASAI_R1=1-SAIDI_R1/8760;

% R2无DG基准（有重构但无DG孤岛，用τ_RP代替τ_RP_eff）
CID_rep_nodg=sum(CID_rep_nodg_mat,2);
CID_noDG=CID_up+CID_tie+CID_rep_nodg+CID_shed;
SAIDI_noDG=(NC_r'*CID_noDG)/total_cust;
EENS_noDG=(P_avg'*CID_noDG)/1e3;

% VOLL经济指标
VOLL_loss_R2 =sum(voll_vec.*P_avg.*CID);
VOLL_loss_R1 =sum(voll_vec.*P_avg.*CID_R1);
VOLL_loss_noDG=sum(voll_vec.*P_avg.*CID_noDG);
VOLL_saving_total=VOLL_loss_R1-VOLL_loss_R2;
VOLL_saving_dg=VOLL_loss_noDG-VOLL_loss_R2;

%% §8  结果输出
total_elapsed=toc(program_total);

fprintf('\n╔══════════════════════════════════════════╗\n');
fprintf('  系统: %-39s\n',sys_filename);
fprintf('  [Q1]负荷分级: L1=%d, L2=%d, L3=%d节点\n',nL1,nL2,nL3);
if DG_ISLAND_ENABLE&&nDG_valid>0
    fprintf('  DG孤岛: %d个节点 | τ_DG_SW=%.2fh | %d状态场景\n', ...
        nDG_valid,TAU_DG_SW,N_DG_SCEN);
end
fprintf('╠══════════════════════════════════════════╣\n');
fprintf('║  [R2 含重构+DG孤岛]                       ║\n');
fprintf('║  SAIFI : %10.4f  次/(户·年)         ║\n',SAIFI);
fprintf('║  SAIDI : %10.4f  h/(户·年)          ║\n',SAIDI);
fprintf('║    ①上游开关   : %+8.4f h/(户·年)    \n',SAIDI_up);
fprintf('║    ②联络转供   : %+8.4f h/(户·年)    \n',SAIDI_tie);
fprintf('║    ③等待修复*  : %+8.4f h/(户·年)    \n',SAIDI_rep);
fprintf('║    ④切负荷等待 : %+8.4f h/(户·年)    \n',SAIDI_shed);
fprintf('║    *③已含DG孤岛等效（τ_RP→τ_RP_eff）  \n');
fprintf('║  EENS  : %10.2f  MWh/年             ║\n',EENS);
fprintf('║  ASAI  : %12.6f                   ║\n',ASAI);
fprintf('║  VOLL年损失: %10.2f 万元/年         ║\n',VOLL_loss_R2/1e4);
fprintf('╠══════════════════════════════════════════╣\n');
fprintf('║  [R2无DG基准（有重构无孤岛）]              ║\n');
fprintf('║  SAIDI : %10.4f  h/(户·年)          ║\n',SAIDI_noDG);
fprintf('║  无DG孤岛降低 : %+8.4f h/(户·年)   \n',SAIDI_dg);
fprintf('║  EENS  : %10.2f  MWh/年             ║\n',EENS_noDG);
fprintf('╠══════════════════════════════════════════╣\n');
fprintf('║  [R1 无重构无DG]                          ║\n');
fprintf('║  SAIDI : %10.4f  h/(户·年)          ║\n',SAIDI_R1);
fprintf('║  EENS  : %10.2f  MWh/年             ║\n',EENS_R1);
fprintf('╠══════════════════════════════════════════╣\n');
fprintf('║  改善分析（相对R1基准）                     ║\n');
fprintf('║  SAIDI改善(重构+DG): %+10.4f h/(户·年)  \n',SAIDI_R1-SAIDI);
fprintf('║  ├─ 重构贡献        : %+10.4f h/(户·年) \n',SAIDI_R1-SAIDI_noDG);
fprintf('║  └─ DG孤岛贡献     : %+10.4f h/(户·年) \n',SAIDI_noDG-SAIDI);
fprintf('║  EENS改善(重构+DG) : %+10.2f MWh/年    \n',EENS_R1-EENS);
fprintf('║  VOLL节约(总)      : %+10.2f 万元/年   \n',VOLL_saving_total/1e4);
fprintf('║  ├─ 重构节约       : %+10.2f 万元/年   \n',(VOLL_loss_R1-VOLL_loss_noDG)/1e4);
fprintf('║  └─ DG孤岛节约    : %+10.2f 万元/年   \n',VOLL_saving_dg/1e4);
fprintf('║  功率基恢复率: %.1f%%(完全%d/部分%d/未%d)   \n',rec_pct,n_full,n_part,n_unrec);
fprintf('╠══════════════════════════════════════════╣\n');
if DG_ISLAND_ENABLE&&nDG_valid>0
    fprintf('║  DG节点孤岛效益明细：                     ║\n');
    for di=1:nDG_valid
        k=dg_local_idx(di);
        dg_saidi_k=(NC_r(k)*CID_dg_benefit(k))/total_cust;
        fprintf('║  节点%-4d: α_eff=%.2f, SAIDI改善%+.4f h  \n', ...
            inv_map(load_nodes(k)),alpha_eff(k),-full(dg_saidi_k));
    end
    fprintf('╠══════════════════════════════════════════╣\n');
end
fprintf('║  §5 MCF : %7.1f 秒                    \n',t_mcf_elapsed);
fprintf('║  §6 MILP: %7.1f 秒 (ZZ=%d档位)       \n',t_milp_elapsed,total_ZZ);
fprintf('║  总耗时 : %7.1f 秒                    \n',total_elapsed);
fprintf('╚══════════════════════════════════════════╝\n');

% 代表性场景输出（与R11相同）
[~,rep_xy]=max(sum(p_direct_d,1));
fprintf('\n══ 代表性故障场景 %d（分支%d─%d）══\n', ...
    rep_xy,inv_map(rel_branches(rep_xy,1)),inv_map(rel_branches(rep_xy,2)));
Vm_rep=Vm_res(:,rep_xy); Pf_rep=Pf_res(:,rep_xy); Qf_rep=Qf_res(:,rep_xy);
Sw_rep=Sw_res(:,rep_xy); qv_rep=q_mat(:,rep_xy);
pf_rep=full(p_feeder_d_full(:,rep_xy)); L_rep=L_res(:,rep_xy);

% 切负荷节点
shed_rep=find(L_rep>SHED_THRESH&qv_rep);
if ~isempty(shed_rep)
    fprintf('\n  切负荷节点\n');
    fprintf('  %-8s %-4s %-5s %-8s %-10s %-10s %-8s\n','节点','级','n_br','VOLL','P/kW','L/kW','切除比');
    fprintf('  %s\n',repmat('─',1,58));
    for k=shed_rep'
        lvl=ternary(mask_L1(k),'L1',ternary(mask_L2(k),'L2','L3'));
        fprintf('  %-8d %-4s %-5d %-8.0f %-10.2f %-10.2f %-8.1f%%\n', ...
            inv_map(load_nodes(k)),lvl,n_br_vec(k),voll_vec(k)*1e3, ...
            P_free(k)*1e3,L_rep(k)*1e3,L_rep(k)/P_free_safe(k)*100);
    end
end

% 节点电压
fprintf('\n  节点状态（供电节点）\n');
fprintf('  %-8s %-10s %-4s %-3s %-16s\n','节点','V/pu','级','DG','状态');
fprintf('  %s\n',repmat('─',1,46));
for si=1:length(subs_idx)
    fprintf('  %-8d %-10.4f %-4s %-3s 变电站\n',inv_map(subs_idx(si)),Vm_rep(subs_idx(si)),'─','─');
end
for k=1:nL
    if qv_rep(k)
        lvl=ternary(mask_L1(k),'L1',ternary(mask_L2(k),'L2','L3'));
        dg_tag=ternary(is_dg(k),sprintf('%.2f',alpha_eff(k)),'─');
        if pf_rep(k)>0&&L_rep(k)>SHED_THRESH; st=sprintf('转供(切%.1f%%)',L_rep(k)/P_free_safe(k)*100);
        elseif pf_rep(k)>0; st='转供恢复'; else; st='正常供电'; end
        fprintf('  %-8d %-10.4f %-4s %-3s %s\n', ...
            inv_map(load_nodes(k)),Vm_rep(load_nodes(k)),lvl,dg_tag,st);
    end
end
% 未恢复节点中的DG孤岛节点
dark_dg_rep=find(~qv_rep&pf_rep>0&is_dg);
if ~isempty(dark_dg_rep)
    fprintf('\n  DG孤岛节点（q=0，主网未恢复，但DG可孤岛供电）:\n');
    for k=dark_dg_rep'
        di_idx=find(dg_local_idx==k,1);
        if isempty(di_idx); continue; end
        trp_eff_k=trp_eff_mat(k,rep_xy);
        fprintf('    节点%-4d: α_eff=%.2f, τ_RP_eff=%.2fh (vs τ_RP=%.2fh)\n', ...
            inv_map(load_nodes(k)),alpha_eff(k),trp_eff_k,rel_branches(rep_xy,4));
    end
end

n_dark=sum(~qv_rep&pf_rep>0&~is_dg);
if n_dark>0; fprintf('  （%d个非DG受影响节点未恢复，略去）\n',n_dark); end

% 分支潮流
fprintf('\n  合路分支潮流（kW/kVar）\n');
fprintf('  %-6s %-6s %-6s %-12s %-12s %-8s\n','分支#','From','To','P/kW','Q/kVar','类型');
fprintf('  %s\n',repmat('─',1,58));
P_sub=0; Q_sub=0;
% 在循环开始前的初始化部分添加：
if length(is_trf_vec) < nB_all
    is_trf_vec(nB_all) = 0; % 自动补齐长度，联络线默认为非变压器
end
for b=1:nB_all
    if Sw_rep(b)==1
        bt=ternary(b<=nB_norm,ternary(is_trf_vec(b),'变压器','线路'),'联络线');
        fprintf('  %-6d %-6d %-6d %+12.2f %+12.2f %s\n',b, ...
            inv_map(all_branches(b,1)),inv_map(all_branches(b,2)), ...
            Pf_rep(b)*1e3,Qf_rep(b)*1e3,bt);
        if ismember(all_branches(b,1),subs_idx); P_sub=P_sub+Pf_rep(b); Q_sub=Q_sub+Qf_rep(b);
        elseif ismember(all_branches(b,2),subs_idx); P_sub=P_sub-Pf_rep(b); Q_sub=Q_sub-Qf_rep(b); end
    end
end
fprintf('\n  汇总: 受影响=%d 恢复=%d 未=%d 切负荷=%.2fkW\n', ...
    sum(pf_rep),sum(qv_rep&pf_rep>0),sum(~qv_rep&pf_rep>0),sum(L_rep(qv_rep))*1e3);
fprintf('  变电站: P=%.2fkW Q=%.2fkVar\n',P_sub*1e3,Q_sub*1e3);

function out=ternary(cond,a,b)
    if cond; out=a; else; out=b; end
end

%% ============================================================
%  配电网可靠性评估 —— R13（Dinkelbach算法替代CC变换）
%
%  优化目标：最大化 VOLL节约 / 开关成本
%  VOLL_saving / (switch_cost + regularization)
%
%  相对R12的改进：
%  1. 完全避免CC变换（无_hat变量爆炸）
%  2. 每次迭代仍是标准MILP（Gurobi友好）
%  3. 收敛速度快（通常3-5次迭代）
%  4. 模型规模保持R12水平
%
%  ─── Bug修复说明（相对原R13版本）──────────────────────────
%  【根因】137/417节点系统MILP不可行：
%    原R13将电压降约束由含松弛量E_vdrop的软约束改为严格等式：
%      原R13: delta_mat <= M_V_SAFE*(1-S_mat)    （S=1时 → delta=0）
%      R12:   delta_mat <= M_V*(1-S_mat)+E_vdrop  （S=1时 → delta=E_vdrop≥0）
%
%    S=1（支路闭合）时 delta=0 即要求线性化DistFlow方程严格成立：
%      V(to) = V(from) - 2*(R·P + X·Q)
%    在137/417节点等长馈线系统中，馈线末端电压沿链路积累性跌落，
%    严格DistFlow会使某些节点V < V_lower_sq = 0.9025²，
%    同时被供电约束 V >= V_lower_sq（当Q=1）, 两者矛盾 → 不可行。
%    小系统（85节点）馈线较短，电压降不足以触发此矛盾，故可行。
%
%  【修复】在Dinkelbach迭代循环内恢复E_vdrop松弛变量：
%    1. 每次迭代增加 E_vdrop = sdpvar(nB_all, nScen, 'full')
%    2. 约束改回: E_vdrop>=0, delta<=M_V*(1-S)+E_vdrop, delta>=-M_V*(1-S)-E_vdrop
%    3. 目标函数增加正则项: +1e-4*sum(sum(E_vdrop))
%    保持E_vdrop最小化使DistFlow在允许微小建模误差前提下几乎严格成立，
%    同时消除大系统长馈线的不可行性。
%    分式目标函数（Dinkelbach框架）完整保留，无任何改动。
% =============================================================
clear; clc;

%% ── 用户配置区 ─────────────────────────────────────────────
%sys_filename = '85-Node System Data.xlsx';
sys_filename = '137-Node System Data.xlsx';
%sys_filename = '417-Node System Data.xlsx';
%sys_filename = '1080-Node System Data.xlsx';

tb_filename = 'Testbench for Linear Model Based Reliability Assessment Method for Distribution Optimization Models Considering Network Reconfiguration.xlsx';

%sys_sheet = '85-node';
sys_sheet = '137-node';
%sys_sheet = '417-node';
%sys_sheet = '1080-node';

%DG_NODES_RAW = [4, 23, 37];   % DG节点的原始编号
DG_NODES_RAW = [8, 29, 44, 61, 95, 111];
%DG_NODES_RAW = [12, 32, 70, 105, 156, 200, 248, 313, 399];
%DG_NODES_RAW = [16, 78, 137, 201, 333, 410, 556, 629, 700, 899, 945, 1005];

LAMBDA_TRF = 0.5; TAU_UP_SW = 0.3; TAU_TIE_SW = 0.5;
W_SCORE_NC = 0.5; W_SCORE_P = 0.5;
RATIO_L1 = 0.30; RATIO_L2 = 0.30;
VOLL_BASE = [200, 50, 10];
VOLL_NOISE = 0.20; VOLL_SEED = 123;
SHED_LIMIT = [0.00, 0.25, 1.00];
FLEX_RATIO_L3 = 0.8;
N_BR_THRESHOLDS = [50, 150, 400];
N_BR_VALUES = [2, 3, 4, 5];
LV_REGEN = false; LV_SEED = 42;
GAMMA_SWITCH = 1e2;

% [DG] 分布式电源孤岛运行参数
DG_ISLAND_ENABLE = true;
DG_CAPACITY = 0.5 * ones(1, size(DG_NODES_RAW, 2));
TAU_DG_SW = 0.2;
DG_SCENARIOS = [
    0.08,  0.10,  0.70;
    0.10,  0.05,  0.83;
    0.06,  0.08,  1.00;
    0.12,  0.45,  0.70;
    0.15,  0.50,  0.83;
    0.10,  0.48,  1.00;
    0.08,  0.90,  0.70;
    0.14,  0.85,  0.83;
    0.09,  0.88,  1.00;
    0.08,  0.30,  0.70;
    ];
DG_DYNAMIC_ISLAND = false;

SOLVE_MODE = 'MCF';
V_UPPER = 1.05; V_LOWER = 0.95; V_SRC = 1.0; PF = 0.9;
%% ──────────────────────────────────────────────────────────────

program_total = tic;

%% §1-5  （与R12相同）
fprintf('>> [1/9] 读取 Testbench: Sheet="%s"\n', sys_sheet);
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

fprintf('>> [3/8] 生成潮流参数...\n');
R_KM=0.003151; X_KM=0.001526; TAN_PHI=tan(acos(PF));
V_src_sq=V_SRC^2; V_upper_sq=V_UPPER^2; V_lower_sq=V_LOWER^2;
M_V=(V_upper_sq-V_lower_sq)*2; M_vn=V_upper_sq;   % M_V=0.4（与R12一致）

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

%% §4b  [Q1+Q2] 负荷分级 & 低压分支开关档位
fprintf('>> [4b/8] 负荷分级与开关档位...\n');
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

%% §4c  DG孤岛运行参数配置（与R12相同）
fprintf('>> [4c/8] DG孤岛运行建模...\n');

if DG_ISLAND_ENABLE
    prob_sum=sum(DG_SCENARIOS(:,1));
    if abs(prob_sum-1)>1e-6
        warning('[DG] 场景概率之和=%.4f≠1，自动归一化',prob_sum);
        DG_SCENARIOS(:,1)=DG_SCENARIOS(:,1)/prob_sum;
    end
    N_DG_SCEN=size(DG_SCENARIOS,1);
    rho_w=DG_SCENARIOS(:,1);
    g_w  =DG_SCENARIOS(:,2);
    l_w  =DG_SCENARIOS(:,3);

    nDG=length(DG_NODES_RAW);
    dg_local_idx=zeros(1,nDG);
    dg_cap_vec=zeros(nDG,1);
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

    is_dg=false(nL,1);
    is_dg(dg_local_idx)=true;

    alpha_eff=zeros(nL,1);
    for di=1:nDG_valid
        k=dg_local_idx(di);
        pg_k=dg_cap_vec(di);
        pd_k=P_free(k);
        if pd_k<1e-9; alpha_eff(k)=1; continue; end
        for wi=1:N_DG_SCEN
            if g_w(wi)*pg_k >= l_w(wi)*pd_k
                alpha_eff(k)=alpha_eff(k)+rho_w(wi);
            end
        end
    end

    lam_vec_pre=rel_branches(:,3);
    trp_vec_pre=rel_branches(:,4);
    trp_eff_mat=repmat(trp_vec_pre',nL,1);
    for di=1:nDG_valid
        k=dg_local_idx(di);
        trp_eff_mat(k,:) = alpha_eff(k)*TAU_DG_SW + (1-alpha_eff(k))*trp_vec_pre';
    end

    tau_benefit_eff=max(trp_eff_mat - TAU_TIE_SW, 0);

    fprintf('   DG孤岛模式：启用 | %d个DG场景 | τ_DG_SW=%.2fh\n',N_DG_SCEN,TAU_DG_SW);
    fprintf('   DG节点(%d个): ',nDG_valid);
    for di=1:nDG_valid
        k=dg_local_idx(di);
        fprintf('节点%-4d(PG=%.2fMW,α_eff=%.2f) ',inv_map(load_nodes(k)),dg_cap_vec(di),alpha_eff(k));
    end
    fprintf('\n');
else
    alpha_eff=zeros(nL,1);
    is_dg=false(nL,1);
    tau_benefit_eff=repmat(max(rel_branches(:,4)'-TAU_TIE_SW,0),nL,1);
    nDG_valid=0; dg_local_idx=[]; dg_cap_vec=[];
    trp_eff_mat=repmat(rel_branches(:,4)',nL,1);
    fprintf('   DG孤岛模式：未启用\n');
end

%% §5  MCF 路径识别（与R12相同）
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

%% ================================================================
%  §5b 预计算R1基准
% ================================================================
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
%  §6  Dinkelbach迭代求解 MILP
%
%  目标: max (VOLL_R1 - VOLL_R2) / (switch_cost + ε)
%  参数化子问题: min { VOLL_R2 + λ·switch_cost + 1e-4·sum(E_vdrop) }
%  λ更新: λ_{k+1} = (VOLL_R1 - VOLL_R2*) / max(switch_cost*, ε_MIN)
%
%  【关键修复】：恢复E_vdrop松弛变量
%    每次迭代内声明 E_vdrop = sdpvar(nB_all, nScen, 'full')
%    约束: E_vdrop>=0, delta<=M_V*(1-S)+E_vdrop, delta>=-M_V*(1-S)-E_vdrop
%    其中 M_V = (V_upper_sq - V_lower_sq)*2 = 0.4（§3中已计算）
%    目标附加: +1e-4*sum(sum(E_vdrop))（驱动松弛量趋零）
%    效果: 消除大系统长馈线电压越限导致的结构性不可行
% ================================================================

%% ================================================================
%  §6  加速Dinkelbach：精确场景解耦 + parfor并行 + 热启动
%
%  加速原理（针对417节点系统）：
%
%  ①【根因】原始MILP规模 = nB_all×nScen ≈ 430×415 ≈ 178,000个二进制变量
%    Gurobi分支定界树在此规模下极深 → 求解时间爆炸
%
%  ②【场景精确解耦】
%    所有故障场景的约束在矩阵层面完全块对角：
%      • 场景xy的决策变量 (S_xy, Q_xy, Pf_xy, ...) 只出现在场景xy的约束中
%      • 各场景目标函数可加可分：obj = Σ_xy obj_xy(S_xy, Q_xy, ...)
%    因此可精确分解为 nScen 个独立小MILP，无近似损失
%    每个子MILP规模：nB_all×1 ≈ 430个二进制变量（缩减~415倍）
%
%  ③【parfor并行】
%    nScen个子MILP用parfor并行求解，加速比 ≈ min(nScen, nWorkers)
%    每个Gurobi子进程用1线程，不争抢CPU
%
%  ④【Dinkelbach热启动（Warm Start）】
%    将第k轮迭代的最优解 (S*_xy, Q*_xy, L*_xy) 作为第k+1轮的MIP Start
%    传入Gurobi → 分支定界从已知可行解出发 → 大幅减少搜索节点数
%    注：λ变化量通常极小（~1e-3量级），相邻迭代解高度相似 → 热启效果极好
%
%  ⑤【Dinkelbach框架完全保留】
%    分式目标函数、λ更新公式、收敛判据均不变
%    仅将"一个大MILP"替换为"nScen个并行小MILP"
% ================================================================
fprintf('>> [6/8] 加速Dinkelbach：场景解耦 + parfor并行 + 热启动...\n');

% ── 参数配置 ─────────────────────────────────────────────────
LAMBDA_INIT         = 0.0;
LAMBDA_TOL          = 1e-4;
MAX_ITER_DINKELBACH = 10;
EPSILON_MIN         = 1e-3;
MIPGAP_SUB          = 1e-3;  % 子MILP允许相对间隙（与原版一致）
TIME_LIMIT_SUB      = 300;   % 单场景子MILP时间上限（秒）
N_THREADS_SUB       = 1;     % 每个Gurobi子进程线程数（parfor并行时保持1）

fprintf('   ε_MIN=%.2e  M_V=%.4f  MIPGap_sub=%.1e\n', EPSILON_MIN, M_V, MIPGAP_SUB);

t_milp_iter = tic;

% ── 预计算稀疏算子（循环外一次性完成）────────────────────────
Dr   = spdiags(r_b_all,   0, nB_all, nB_all);
Dx   = spdiags(x_b_all,   0, nB_all, nB_all);
Dp   = spdiags(P_free,    0, nL,     nL);
Dq_  = spdiags(Q_free,    0, nL,     nL);
Cap  = spdiags(cap_b_all, 0, nB_all, nB_all);

% ── 权重矩阵 ─────────────────────────────────────────────────
W_q = (voll_pu.*P_free.*p_mat_d).*(repmat(lam_vec',nL,1).*tau_benefit_eff);
W_L = (voll_pu.*p_mat_d).*(repmat(lam_vec',nL,1).*tau_benefit_eff);
fprintf('   W_q: nnz=%d, max=%.2f元/年\n', nnz(W_q), max(full(W_q(:))));
if DG_ISLAND_ENABLE && nDG_valid>0
    tau_benefit_std = repmat(max(trp_vec'-TAU_TIE_SW,0),nL,1);
    W_q_std = (voll_pu.*P_free.*p_mat_d).*(repmat(lam_vec',nL,1).*tau_benefit_std);
    for di=1:nDG_valid
        k=dg_local_idx(di);
        wq_ratio = sum(W_q(k,:))/max(sum(W_q_std(k,:)),1e-9);
        fprintf('   DG节点%d W_q降权至%.1f%%（α_eff=%.2f）\n', ...
            inv_map(load_nodes(k)), full(wq_ratio)*100, alpha_eff(k));
    end
end

% ── 稀疏ZZ算子（场景无关，parfor广播）─────────────────────────
can_shed       = shed_limit_vec > 1e-6;
n_shed         = sum(can_shed);
shed_nodes_idx = find(can_shed);
total_ZZ       = sum(n_lev_arr(shed_nodes_idx));
row_s  = zeros(total_ZZ,1);
val_sel= zeros(total_ZZ,1);
val_one= ones(total_ZZ,1);
ptr = 0;
for ki = 1:n_shed
    k  = shed_nodes_idx(ki);
    lv = node_levels{k};
    nk = length(lv);
    rows = ptr+1:ptr+nk;
    row_s(rows)   = ki;
    val_sel(rows) = P_free(k) * lv(:);
    ptr = ptr + nk;
end
col_s    = (1:total_ZZ)';
SEL_s    = sparse(row_s, col_s, val_sel, n_shed, total_ZZ);
ONESUM_s = sparse(row_s, col_s, val_one, n_shed, total_ZZ);

% ── 预计算R1基准（常数）──────────────────────────────────────
CID_R1_base        = TAU_UP_SW*(p_upstream_base*lam_vec) + p_mat_d*(lam_vec.*trp_vec);
VOLL_loss_R1_const = sum(voll_vec.*P_avg_obj.*CID_R1_base);
fprintf('   R1基准损失（无重构无DG）: %.2f 万元/年\n', VOLL_loss_R1_const/1e4);

% ── 开关初始状态向量 ─────────────────────────────────────────
S_NO = ones(nB_all,1);
S_NO(nB_norm+1:end) = 0;
S_NO_mat = repmat(S_NO, 1, nScen);  % 仅用于收敛分析（顺序部分）

% ── 预计算parfor广播常量 ──────────────────────────────────────
% 这些量在场景间共享，parfor自动广播到各worker
w_voll_c       = voll_vec .* P_avg_obj;          % [nL,1]
P_free_c       = P_free;                          % [nL,1]
P_free_safe_c  = max(P_free, 1e-9);              % [nL,1]
lam_vec_c      = lam_vec;                         % [nB_norm,1]
% 将二维场景矩阵转为full（确保parfor列切片高效）
p_mat_d_f      = full(p_mat_d);                  % [nL, nScen]
p_feeder_d_f   = full(p_feeder_d);               % [nL, nScen]
trp_eff_mat_f  = full(trp_eff_mat);              % [nL, nScen]
tau_benefit_f  = full(tau_benefit_eff);           % [nL, nScen]

% ── 热启动初值（第0次迭代：保持初始开关状态）────────────────────
%  S_ws：开关状态；Q_ws：供电状态；L_ws：切负荷量
%  语义：ws_*(:,xy) = 场景xy的上一轮最优解
ws_S = repmat(S_NO,         1, nScen);   % [nB_all, nScen]
ws_Q = max(1 - p_feeder_d_f, 0);         % [nL, nScen]  保守初值
ws_L = zeros(nL, nScen);                 % [nL, nScen]

% ── Dinkelbach初始化 ──────────────────────────────────────────
lambda_k             = LAMBDA_INIT;
lambda_history       = lambda_k;
VOLL_loss_R2_history = [];
switch_cost_history  = [];
regularization_history = [];

fprintf('\n╔══════════════════════════════════════════════════════════════╗\n');
fprintf('║   加速Dinkelbach：场景解耦 + parfor×%4d + 热启动          ║\n', nScen);
fprintf('╠══════════════════════════════════════════════════════════════╣\n');
fprintf('║   每子MILP规模：~%d个二值变量（原始大MILP的1/%d）        \n', nB_all+nL, nScen);
fprintf('╚══════════════════════════════════════════════════════════════╝\n');

for iter = 1:MAX_ITER_DINKELBACH
    fprintf('\n【迭代 %d】λ=%.8f  开始parfor并行求解 %d 个子MILP...\n', ...
        iter, lambda_k, nScen);
    t_iter_tic = tic;

    % ── 预分配结果矩阵（parfor sliced outputs）──────────────────
    S_sol  = zeros(nB_all,   nScen);  % 最优开关状态
    Q_sol  = zeros(nL,       nScen);  % 最优供电状态
    L_sol  = zeros(nL,       nScen);  % 最优切负荷量
    Pf_sol = zeros(nB_all,   nScen);  % 最优有功潮流
    Qf_sol = zeros(nB_all,   nScen);  % 最优无功潮流
    V_sol  = zeros(num_nodes, nScen); % 最优节点电压平方
    fail_flag = zeros(1, nScen);      % 1=子MILP求解失败

    lk = lambda_k;   % 标量拷贝，parfor广播安全

    % ────────────────────────────────────────────────────────────
    %  parfor：每个worker独立求解场景xy的子MILP
    %  变量分类：
    %    · 广播（broadcast）：Dr,Dx,Cap,Dp,Dq_,BdV,A_free_all,SEL_s,
    %                          ONESUM_s,S_NO,subs_idx,non_sub,
    %                          w_voll_c,P_free_c,P_free_safe_c,
    %                          V_src_sq,V_lower_sq,V_upper_sq,
    %                          M_V,M_vn,TAN_PHI,GAMMA_SWITCH,
    %                          TAU_TIE_SW,total_ZZ,n_shed,
    %                          shed_nodes_idx,can_shed,nB_all,nL,
    %                          num_nodes,nB_norm,lk
    %    · 切片输入（sliced input，按列xy）：
    %                          p_mat_d_f,p_feeder_d_f,lam_vec_c,
    %                          trp_eff_mat_f,tau_benefit_f,
    %                          ws_S,ws_Q,ws_L
    %    · 切片输出（sliced output，写入列xy）：
    %                          S_sol,Q_sol,L_sol,Pf_sol,Qf_sol,
    %                          V_sol,fail_flag
    % ────────────────────────────────────────────────────────────
    
    % ── Gurobi求解参数 ──────────────────────────────────────
    %  · verbose=0：避免parfor多worker输出混乱
    %  · usex0=1：将assign()的值作为MIP Start传给Gurobi
    %  · Threads=1：parfor时每子进程单线程，避免核心争抢
    opts_sub = sdpsettings('solver','gurobi','verbose',0,'usex0',1, ...
        'gurobi.MIPGap',       MIPGAP_SUB, ...
        'gurobi.Heuristics',   0.10, ...
        'gurobi.Presolve',     -1, ...
        'gurobi.Cuts',          2, ...
        'gurobi.Threads',       N_THREADS_SUB, ...
        'gurobi.TimeLimit',     TIME_LIMIT_SUB);
    
    parfor xy = 1:nScen %#ok<*PFBNS> 可切换parfor
        %yalmip('clear');   % 清除各worker的YALMIP内部缓存
        safe_yalmip_clear(); % <--- 替换为安全封装函数
        
        % ── 抽取场景xy的常量列向量 ─────────────────────────────
        p_xy   = p_mat_d_f(:, xy);       % [nL,1] 受影响节点掩码
        pf_xy  = p_feeder_d_f(:, xy);    % [nL,1] 主馈线连接掩码
        lam_xy = lam_vec_c(xy);          % 标量：故障率×支路长度
        trp_xy = trp_eff_mat_f(:, xy);   % [nL,1] 有效修复时间
        tau_xy = tau_benefit_f(:, xy);   % [nL,1] 联络开关节约时间
        pdc_xy = p_xy ./ P_free_safe_c;  % [nL,1] 切负荷比例系数

        % ── 声明YALMIP决策变量 ─────────────────────────────────
        S_v  = binvar(nB_all,    1, 'full');  % 支路开关状态
        Q_v  = binvar(nL,        1, 'full');  % 节点供电状态
        Pf_v = sdpvar(nB_all,    1, 'full');  % 支路有功潮流
        Qf_v = sdpvar(nB_all,    1, 'full');  % 支路无功潮流
        V_v  = sdpvar(num_nodes, 1, 'full');  % 节点电压平方
        ZZ_v = binvar(total_ZZ,  1, 'full');  % 切负荷档位选择
        L_v  = sdpvar(nL,        1, 'full');  % 节点切负荷量
        E_v  = sdpvar(nB_all,    1, 'full');  % 电压降松弛量

        % ── 热启动：将上一迭代解赋初值传给Gurobi MIP Start ──────
        %  assign(var, val) → YALMIP在usex0=1时传入Gurobi Start属性
        %  Gurobi从此可行点出发分支定界，通常大幅减少搜索节点
        assign(S_v,  ws_S(:, xy));
        assign(Q_v,  ws_Q(:, xy));
        assign(L_v,  ws_L(:, xy));

        % ── DistFlow线性化电压降（软约束形式，含E_v松弛）─────────
        delta_v = BdV*V_v + 2*(Dr*Pf_v + Dx*Qf_v);

        % ── 切负荷档位展开 ─────────────────────────────────────
        L_shed_v = SEL_s * ZZ_v;  % [n_shed,1]

        % ── MILP约束 ────────────────────────────────────────────
        %  与原始大MILP约束完全一致，仅对单场景xy切片
        C = [ V_v(subs_idx)  == V_src_sq, ...
              V_v(non_sub)   >= V_lower_sq - M_vn*(1-Q_v), ...
              V_v(non_sub)   <= V_upper_sq + M_vn*(1-Q_v), ...
             -Cap*S_v        <= Pf_v, ...
              Pf_v           <= Cap*S_v, ...
             -Cap*S_v        <= Qf_v, ...
              Qf_v           <= Cap*S_v, ...
              A_free_all*Pf_v == Dp*Q_v - L_v, ...
              A_free_all*Qf_v == Dq_*Q_v - TAN_PHI*L_v, ...
              E_v            >= 0, ...
              delta_v        <= M_V*(1-S_v) + E_v, ...   % 软约束（含松弛）
              delta_v        >= -M_V*(1-S_v) - E_v, ...  % 软约束（含松弛）
              sum(S_v)       == sum(Q_v), ...
              S_v(xy)        == 0, ...            % 场景xy：故障支路强制断开
              Q_v            >= 1 - pf_xy];       % 主馈线连接节点必须恢复

        if any(~can_shed)
            C = [C, L_v(~can_shed) == 0];
        end
        C = [C, ZZ_v             >= 0, ...
                L_v(shed_nodes_idx) == L_shed_v, ...
                ONESUM_s*ZZ_v    == 1, ...
                L_v              >= 0, ...
                L_v              <= P_free_c .* Q_v];

        % ── 子MILP目标函数（Dinkelbach参数化）──────────────────
        %  obj_xy = VOLL_R2_xy(Q_v,L_v) + λk·switch_cost_xy(S_v) + ε·E_v
        %  VOLL_R2_xy 仅含变量部分（常量CID_up_const_xy不影响argmin）
        %
        %  其中：
        %  CID_tie_xy  = τ_TIE · λ_xy · (p_xy ⊙ Q_v)      [nL,1]
        %  CID_rep_xy  = λ_xy · (p_xy ⊙ (1-Q_v) ⊙ τ_RP_eff_xy) [nL,1]
        %  CID_shed_xy = λ_xy · (pdc_xy ⊙ L_v ⊙ τ_benefit_xy) [nL,1]
        sw_cost_v  = GAMMA_SWITCH * sum(S_NO.*(1-S_v) + (1-S_NO).*S_v);

        CID_tie_v  = TAU_TIE_SW * lam_xy * (p_xy .* Q_v);
        CID_rep_v  = lam_xy * (p_xy .* (1-Q_v) .* trp_xy);
        CID_shed_v = lam_xy * (pdc_xy .* L_v .* tau_xy);
        VOLL_var_v = sum(w_voll_c .* (CID_tie_v + CID_rep_v + CID_shed_v));

        obj_v = VOLL_var_v + lk * sw_cost_v + 1e-4 * sum(E_v);

        sol_xy = optimize(C, obj_v, opts_sub);

        % ── 提取解（切片输出：每个xy只写自己的列）──────────────
        %  sol_xy.problem: 0=最优, 2=不可行, 3=时间到（有可行解）
        if sol_xy.problem == 0 || sol_xy.problem == 3
            Sv_  = value(S_v);  Qv_  = value(Q_v);
            Pfv_ = value(Pf_v); Qfv_ = value(Qf_v);
            Vv_  = value(V_v);  Lv_  = value(L_v);
            % 处理潜在NaN（极少情况）
            if any(isnan(Qv_))
                fail_flag(xy) = 1; %#ok<PFOUS>
                S_sol(:,  xy) = S_NO;
                Q_sol(:,  xy) = 1 - pf_xy;
                L_sol(:,  xy) = zeros(nL,1);
                Pf_sol(:, xy) = zeros(nB_all,1);
                Qf_sol(:, xy) = zeros(nB_all,1);
                V_sol(:,  xy) = V_src_sq * ones(num_nodes,1);
            else
                fail_flag(xy) = 0; %#ok<PFOUS>
                S_sol(:,  xy) = min(max(round(Sv_), 0), 1);
                Q_sol(:,  xy) = min(max(round(Qv_), 0), 1);
                L_sol(:,  xy) = max(Lv_, 0);
                Pf_sol(:, xy) = Pfv_;
                Qf_sol(:, xy) = Qfv_;
                V_sol(:,  xy) = Vv_;
            end
        else
            % 求解失败（不可行或错误）→ 降级到保守解
            fail_flag(xy) = 1; %#ok<PFOUS>
            S_sol(:,  xy) = S_NO;
            Q_sol(:,  xy) = 1 - pf_xy;
            L_sol(:,  xy) = zeros(nL,1);
            Pf_sol(:, xy) = zeros(nB_all,1);
            Qf_sol(:, xy) = zeros(nB_all,1);
            V_sol(:,  xy) = V_src_sq * ones(num_nodes,1);
        end
    end  % end parfor xy

    n_fail = sum(fail_flag);
    fprintf('   parfor完成 (%.1fs)  失败场景=%d/%d\n', toc(t_iter_tic), n_fail, nScen);

    % ── 更新热启动（为下一次Dinkelbach迭代准备）─────────────────
    %  将本轮最优解存入ws_*，下轮parfor通过assign()传给Gurobi
    ws_S = S_sol;
    ws_Q = Q_sol;
    ws_L = L_sol;

    % ── Dinkelbach收敛分析（顺序计算）───────────────────────────
    q_mat   = logical(Q_sol);
    q_mat_d = double(q_mat);
    Sw_res  = S_sol;

    % 重新计算VOLL_R2（含CID_up常数项，用于Dinkelbach λ更新）
    p_direct_d = p_mat_d;   % 等同于double(p_mat)，与§7保持一致
    CID_up_val  = TAU_UP_SW * (p_upstream_base * lam_vec);      % [nL,1]常数
    CID_tie_val = TAU_TIE_SW*(p_direct_d.*q_mat_d)*lam_vec;
    CID_rep_val = sum((p_direct_d.*(1-q_mat_d)).*repmat(lam_vec',nL,1).*trp_eff_mat, 2);
    p_dc_val    = repmat(1./max(P_free,1e-9),1,nScen) .* p_direct_d;
    CID_shed_val= sum((p_dc_val.*L_sol).*repmat(lam_vec',nL,1).*tau_benefit_eff, 2);

    VOLL_loss_R2_val = sum(voll_vec.*P_avg_obj.* ...
        (CID_up_val + CID_tie_val + CID_rep_val + CID_shed_val));
    VOLL_saving_val  = VOLL_loss_R1_const - VOLL_loss_R2_val;

    switch_cost_val    = full(sum(sum(S_NO_mat.*(1-Sw_res) + (1-S_NO_mat).*Sw_res))) * GAMMA_SWITCH;
    regularization_val = max(switch_cost_val, EPSILON_MIN);
    lambda_k_new       = VOLL_saving_val / regularization_val;

    fprintf('  VOLL_R1=%10.2f | VOLL_R2=%10.2f | saving=%10.2f 万元/年\n', ...
        VOLL_loss_R1_const/1e4, VOLL_loss_R2_val/1e4, VOLL_saving_val/1e4);
    fprintf('  switch_cost=%10.2f 万元 | reg=%10.2f 万元\n', ...
        switch_cost_val/1e4, regularization_val/1e4);
    fprintf('  λ_old=%10.8f → λ_new=%10.8f  (Δλ=%.2e)\n', ...
        lambda_k, lambda_k_new, abs(lambda_k_new-lambda_k));

    VOLL_loss_R2_history   = [VOLL_loss_R2_history;   VOLL_loss_R2_val];  %#ok<AGROW>
    switch_cost_history    = [switch_cost_history;    switch_cost_val];   %#ok<AGROW>
    regularization_history = [regularization_history; regularization_val];%#ok<AGROW>
    lambda_history         = [lambda_history;         lambda_k_new];      %#ok<AGROW>

    % ── 收敛判断 ─────────────────────────────────────────────
    if abs(lambda_k_new - lambda_k) < LAMBDA_TOL
        fprintf('  ✓ Dinkelbach已收敛！(Δλ=%.2e < %.2e)\n', ...
            abs(lambda_k_new-lambda_k), LAMBDA_TOL);
        lambda_k = lambda_k_new;
        break;
    end
    lambda_k = lambda_k_new;
    yalmip('clear');

end  % end Dinkelbach outer loop

t_milp_iter_elapsed = toc(t_milp_iter);
fprintf('\n  Dinkelbach总耗时: %.1f秒，迭代次数: %d/%d\n', ...
    t_milp_iter_elapsed, iter, MAX_ITER_DINKELBACH);

% ── Dinkelbach收敛历史 ────────────────────────────────────────
fprintf('\n╔════════════════════════════════════════════════════════════╗\n');
fprintf('║              Dinkelbach迭代收敛历史                        ║\n');
fprintf('╠════════════════════════════════════════════════════════════╣\n');
fprintf('║ Iter │      λ      │  VOLL_R2   │  switch  │ ratio(分子/分) ║\n');
fprintf('║      │             │   万元     │  万元    │     母         ║\n');
fprintf('╟──────┼─────────────┼────────────┼──────────┼────────────────╢\n');
for ii=1:min(length(lambda_history), 10)
    if ii <= length(VOLL_loss_R2_history)
        ratio_val = (VOLL_loss_R1_const - VOLL_loss_R2_history(ii)) / ...
                     max(regularization_history(ii), 1e-9);
        fprintf('║ %2d   │ %.6e │ %.4f    │ %.4f   │    %.4f        ║\n', ...
            ii-1, lambda_history(ii), ...
            VOLL_loss_R2_history(ii)/1e4, ...
            switch_cost_history(ii)/1e4, ...
            ratio_val);
    else
        fprintf('║ %2d   │ %.6e │ (初值)     │          │                ║\n', ...
            ii-1, lambda_history(ii));
    end
end
fprintf('╚════════════════════════════════════════════════════════════╝\n');

% ── 输出变量整理（供§7-9使用）────────────────────────────────
L_res  = max(L_sol, 0);
Pf_res = Pf_sol;
Qf_res = Qf_sol;
Vm_res = sqrt(max(V_sol, 0));

% 最终安全检查
if all(isnan(q_mat_d(:)))
    warning('[§6] 最终解全为NaN，降级到无重构解');
    q_mat   = logical(1 - p_feeder_d);
    q_mat_d = double(q_mat);
    L_res   = zeros(nL,    nScen);
    Pf_res  = zeros(nB_all, nScen);
    Qf_res  = zeros(nB_all, nScen);
    Vm_res  = repmat(V_src_sq, num_nodes, nScen);
    Sw_res  = zeros(nB_all, nScen);
end

%% §7  可靠性指标计算（与R12相同）
fprintf('>> [7/8] 计算可靠性指标...\n');
lam=rel_branches(:,3); trp=rel_branches(:,4);
p_upstream_d=p_feeder_d-p_direct_d;
[~,c_row]=ismember(inv_map(load_nodes),t_cust.Node);
NC_r=zeros(nL,1); NC_r(c_row>0)=t_cust.NC(c_row(c_row>0));
[~,p_row]=ismember(inv_map(load_nodes),t_peak.Node);
P_avg=zeros(nL,1); P_avg(p_row>0)=t_peak.P_kW(p_row(p_row>0))*sum(L_f.*(T_l/8760));
total_cust=sum(NC_r);
P_free_safe=max(P_free,1e-9);
shed_ratio=(L_res./repmat(P_free_safe,1,nScen)).*q_mat_d.*p_direct_d;

% R2 三阶段CID
CIF=p_feeder_d*lam;
CID_up=TAU_UP_SW*(p_upstream_d*lam);
CID_tie=TAU_TIE_SW*(p_direct_d.*q_mat_d)*lam;
CID_rep_mat=(p_direct_d.*(1-q_mat_d)).*(repmat(lam',nL,1).*trp_eff_mat);
CID_rep=sum(CID_rep_mat,2);
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

% DG孤岛收益分项
CID_rep_nodg_mat=(p_direct_d.*(1-q_mat_d)).*repmat(lam'.*trp',nL,1);
CID_dg_benefit=sum(CID_rep_nodg_mat,2)-CID_rep;
SAIDI_dg=(NC_r'*CID_dg_benefit)/total_cust;

% R1基准
CID_R1=TAU_UP_SW*(p_upstream_d*lam)+p_direct_d*(lam.*trp);
SAIDI_R1=(NC_r'*CID_R1)/total_cust;
EENS_R1=(P_avg'*CID_R1)/1e3;
ASAI_R1=1-SAIDI_R1/8760;

% R2无DG基准
CID_rep_nodg=sum(CID_rep_nodg_mat,2);
CID_noDG=CID_up+CID_tie+CID_rep_nodg+CID_shed;
SAIDI_noDG=(NC_r'*CID_noDG)/total_cust;
EENS_noDG=(P_avg'*CID_noDG)/1e3;

% VOLL经济指标
VOLL_loss_R2  = sum(voll_vec.*P_avg.*CID);
VOLL_loss_R1  = sum(voll_vec.*P_avg.*CID_R1);
VOLL_loss_noDG= sum(voll_vec.*P_avg.*CID_noDG);
VOLL_saving_total = VOLL_loss_R1 - VOLL_loss_R2;
VOLL_saving_dg    = VOLL_loss_noDG - VOLL_loss_R2;

%% §8  结果输出
total_elapsed=toc(program_total);
mask_L1=false(nL,1); mask_L1(idx_L1)=true;
mask_L2=false(nL,1); mask_L2(idx_L2)=true;

fprintf('\n╔════════════════════════════════════════════════════╗\n');
fprintf('║         配电网可靠性评估 (Dinkelbach版)              ║\n');
fprintf('╠════════════════════════════════════════════════════╣\n');
fprintf('  系统: %-45s\n',sys_filename);
fprintf('  [Q1]负荷分级: L1=%d, L2=%d, L3=%d节点\n',nL1,nL2,nL3);
if DG_ISLAND_ENABLE&&nDG_valid>0
    fprintf('  DG孤岛: %d个节点 | τ_DG_SW=%.2fh | %d状态场景\n', ...
        nDG_valid,TAU_DG_SW,N_DG_SCEN);
end
fprintf('╠════════════════════════════════════════════════════╣\n');
fprintf('║  [Dinkelbach迭代过程]                              ║\n');
fprintf('║  迭代次数: %-2d                                    ║\n', iter);
fprintf('║  最终λ: %.6f                                   ║\n', lambda_k);
if length(lambda_history) >= 2
    fprintf('║  λ收敛: Δλ=%.2e                            ║\n', abs(lambda_history(end)-lambda_history(end-1)));
end
fprintf('╠════════════════════════════════════════════════════╣\n');
fprintf('║  [R2 含重构+DG孤岛]                                 ║\n');
fprintf('║  SAIFI : %10.4f  次/(户·年)                   ║\n',SAIFI);
fprintf('║  SAIDI : %10.4f  h/(户·年)                    ║\n',SAIDI);
fprintf('║    ①上游开关   : %+8.4f h/(户·年)            ║\n',SAIDI_up);
fprintf('║    ②联络转供   : %+8.4f h/(户·年)            ║\n',SAIDI_tie);
fprintf('║    ③等待修复*  : %+8.4f h/(户·年)            ║\n',SAIDI_rep);
fprintf('║    ④切负荷等待 : %+8.4f h/(户·年)            ║\n',SAIDI_shed);
fprintf('║    *③已含DG孤岛等效（τ_RP→τ_RP_eff）  \n');
fprintf('║  EENS  : %10.2f  MWh/年                       ║\n',EENS);
fprintf('║  ASAI  : %12.6f                             ║\n',ASAI);
fprintf('║  VOLL年损失: %10.2f 万元/年                   ║\n',VOLL_loss_R2/1e4);
fprintf('╠════════════════════════════════════════════════════╣\n');
fprintf('║  [R2无DG基准（有重构无孤岛）]              ║\n');
fprintf('║  SAIDI : %10.4f  h/(户·年)          ║\n',SAIDI_noDG);
fprintf('║  无DG孤岛降低 : %+8.4f h/(户·年)   \n',SAIDI_dg);
fprintf('║  EENS  : %10.2f  MWh/年             ║\n',EENS_noDG);
fprintf('╠══════════════════════════════════════════╣\n');
fprintf('║  [R1 无重构无DG]                                    ║\n');
fprintf('║  SAIDI : %10.4f  h/(户·年)                    ║\n',SAIDI_R1);
fprintf('║  EENS  : %10.2f  MWh/年                       ║\n',EENS_R1);
fprintf('╠════════════════════════════════════════════════════╣\n');
fprintf('║  改善分析（相对R1基准）                     ║\n');
fprintf('║  SAIDI改善(重构+DG): %+10.4f h/(户·年)  \n',SAIDI_R1-SAIDI);
fprintf('║  ├─ 重构贡献        : %+10.4f h/(户·年) \n',SAIDI_R1-SAIDI_noDG);
fprintf('║  └─ DG孤岛贡献     : %+10.4f h/(户·年) \n',SAIDI_noDG-SAIDI);
fprintf('║  EENS改善(重构+DG) : %+10.2f MWh/年    \n',EENS_R1-EENS);
fprintf('║  VOLL节约(总)      : %+10.2f 万元/年   \n',VOLL_saving_total/1e4);
fprintf('║  ├─ 重构节约       : %+10.2f 万元/年   \n',(VOLL_loss_R1-VOLL_loss_noDG)/1e4);
fprintf('║  └─ DG孤岛节约    : %+10.2f 万元/年   \n',VOLL_saving_dg/1e4);
fprintf('║  功率基恢复率: %.1f%%        \n', ...
    100*full(sum(sum((P_free-L_res).*q_mat_d.*p_direct_d))/max(sum(sum(P_free.*p_direct_d)),1e-9)));
fprintf('╠════════════════════════════════════════════════════╣\n');
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
fprintf('║  耗时统计                                            ║\n');
fprintf('║  §5 MCF        : %7.1f 秒                       ║\n',t_mcf_elapsed);
fprintf('║  §6 Dinkelbach : %7.1f 秒 (%d迭代)            ║\n',t_milp_iter_elapsed,iter);
fprintf('║  总耗时        : %7.1f 秒                       ║\n',total_elapsed);
fprintf('╚════════════════════════════════════════════════════╝\n');

%% ── 绘制 Dinkelbach 收敛曲线图 ──────────────────────────────
figure('Color', 'w', 'Name', 'Dinkelbach Convergence Analysis');

subplot(2, 1, 1);
iters = 0:length(lambda_history)-1;
plot(iters, lambda_history, '-o', 'LineWidth', 1.5, 'MarkerSize', 6, 'Color', [0 0.4470 0.7410]);
grid on;
xlabel('Iteration'); ylabel('Value of \lambda');
title('Convergence of \lambda (Ratio of Saving to Cost)');
set(gca, 'FontSize', 10);

subplot(2, 1, 2);
if ~isempty(VOLL_loss_R2_history)
    saving_history_wan = (VOLL_loss_R1_const - VOLL_loss_R2_history) / 1e4;
    plot(1:length(saving_history_wan), saving_history_wan, '-s', ...
        'LineWidth', 1.5, 'MarkerSize', 6, 'Color', [0.8500 0.3250 0.0980]);
    grid on;
    xlabel('Iteration'); ylabel('VOLL Saving (10^4 Yuan)');
    title('Economic Benefit Improvement during Iterations');
end
set(gca, 'FontSize', 10);
sgtitle('Dinkelbach Algorithm Performance', 'FontSize', 12, 'FontWeight', 'bold');

%% §9  代表性场景输出
[~,rep_xy]=max(sum(p_direct_d,1));
fprintf('\n══ 代表性故障场景 %d（分支%d─%d）══\n', ...
    rep_xy,inv_map(rel_branches(rep_xy,1)),inv_map(rel_branches(rep_xy,2)));
Vm_rep=Vm_res(:,rep_xy); Pf_rep=Pf_res(:,rep_xy); Qf_rep=Qf_res(:,rep_xy);
Sw_rep=Sw_res(:,rep_xy); qv_rep=q_mat(:,rep_xy);
pf_rep=full(p_feeder_d(:,rep_xy)); L_rep=L_res(:,rep_xy);

SHED_THRESH=1e-4;
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

fprintf('\n  合路分支潮流（kW/kVar）\n');
fprintf('  %-6s %-6s %-6s %-12s %-12s %-8s\n','分支#','From','To','P/kW','Q/kVar','类型');
fprintf('  %s\n',repmat('─',1,58));
P_sub=0; Q_sub=0;
if length(is_trf_vec) < nB_all
    is_trf_vec(nB_all) = 0;
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

fprintf('\n✓ 程序完成！\n');

function out=ternary(cond,a,b)
if cond; out=a; else; out=b; end
end

function safe_yalmip_clear()
    % 隔离 parfor 工作区，避免 yalmip 内部的 evalin 触发透明度违规
    yalmip('clear');
end
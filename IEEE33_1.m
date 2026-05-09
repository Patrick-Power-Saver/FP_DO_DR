%多时段socp-opf,Sb=100MVA,Vb=12.66KV
clear 
clc
tic
warning off
%% 1.设参
mpc = IEEE33BW;
pload = mpc.Pload;%节点有功负荷
qload = mpc.Qload;%节点无功负荷
branch = mpc.branch;
branch(:,3) = branch(:,3)*100/(12.66^2);%求阻抗标幺值
r = real(branch(:,3));
x = imag(branch(:,3));
T = 24;%时段数为24小时
nb = 33;%节点数
nl = 32;%支路数
nwt = 3;%3个风机
npv = 4;%4个光伏
ness = 2;%ESS数
upstream = zeros(nb,nl);
dnstream = zeros(nb,nl);
for i = 1:nl
    upstream(i,i) = 1;
end
for i=[1:16,18:20,22:23,25:31]
    dnstream(i,i+1) = 1;
end
dnstream(1,18) = 1;
dnstream(2,22) = 1;
dnstream(5,25) = 1;
dnstream(33,1) = 1;

Vmax = [1.05*1.05*ones(nb-1,T)
        1.05*1.05*ones(1,T)];
Vmin = [0.95*0.95*ones(nb-1,T)
        1.05*1.05*ones(1,T)];
Pgmax = [zeros(nb-1,T)
         ones(1,T)];
Qgmax = [zeros(nb-1,T)
         ones(1,T)];
%% 2.设变量
V = sdpvar(nb,T);%电压的平方
I = sdpvar(nl,T);%电流的平方
P = sdpvar(nl,T);%线路有功
Q = sdpvar(nl,T);%线路无功
Pg = sdpvar(nb,T);%发电机有功
Qg = sdpvar(nb,T);%发电机无功
p_wt = sdpvar(nwt,T);%风机有功
p_pv = sdpvar(npv,T);%光伏有功
p_dch = sdpvar(ness,T);%ESS放电功率
p_ch = sdpvar(ness,T);%ESS充电功率
u_dch = binvar(ness,T);%ESS放电状态
u_ch = binvar(ness,T);%ESS充电状态
E_ess = sdpvar(ness,T);%ESS的电量

%% 3.设约束
C = [];
%% 储能装置（ESS）约束
%充放电状态约束
C = [C, u_dch + u_ch <= 1];%表示充电，放电，不充不放三种状态
%功率约束
C = [C, 0 <= p_dch(1,:) <= u_dch(1,:)*0.3];
C = [C, 0 <= p_dch(2,:) <= u_dch(2,:)*0.2];
C = [C, 0 <= p_ch(1,:) <= u_ch(1,:)*0.3];
C = [C, 0 <= p_ch(2,:) <= u_ch(2,:)*0.2];
%容量约束
for t = 1:23
        C = [C, E_ess(:,t+1) == E_ess(:,t) + 0.9*p_ch(:,t) - 1.11*p_dch(:,t)]; 
end
C = [C, 0.18 <= E_ess(1,:) <= 1.8];
C = [C, 0.10 <= E_ess(2,:) <= 1.0];
%投入节点选择
P_dch = [zeros(14,T);p_dch(1,:);zeros(16,T);p_dch(2,:);zeros(1,T)];
P_ch = [zeros(14,T);p_ch(1,:);zeros(16,T);p_ch(2,:);zeros(1,T)];
%% 风机+光伏约束
C = [C, 0 <= p_wt,p_wt <= 1];
P_wt = [zeros(2,T);p_wt(1,:);zeros(13,T);p_wt(2,:);zeros(3,T);p_wt(3,:);zeros(12,T)];
C = [C, 0 <= p_pv,p_pv <= 0.5];
P_pv = [zeros(15,T);p_pv(1,:);zeros(2,T);p_pv(2,:);zeros(10,T);p_pv(3,:);p_pv(4,:);zeros(2,T)];

%% 潮流约束
%节点功率约束
Pin = -upstream*P + upstream*(I.*(r*ones(1,T))) + dnstream*P;%节点注入有功
Qin = -upstream*Q + upstream*(I.*(x*ones(1,T))) + dnstream*Q;%节点注入无功
C = [C, Pin + pload - Pg - P_wt - P_pv - P_dch + P_ch == 0];
C = [C, Qin + qload - Qg == 0];
%欧姆定律约束
C = [C, V(branch(:,2),:) == V(branch(:,1),:) - 2*(r*ones(1,T)).*P - 2*(x*ones(1,T)).*Q + ((r.^2+x.^2)*ones(1,T)).*I];
%二阶锥约束
C = [C, V(branch(:,1),:).*I >= P.^2 + Q.^2];
%% 通用约束
%节点电压约束
C = [C, Vmin <= V,V <= Vmax];
%发电机功率约束
C = [C, -Pgmax <= Pg,Pg <= Pgmax,-Qgmax <= Qg,Qg <= Qgmax];
%支路电流约束
C = [C, 0 <= I,I <= 0.11];
%支路功率约束
C = [C, -0.11<=P,P <= 0.11,-0.11 <= Q,Q <= 0.11];
%% 4.设目标函数
objective = 500*sum(sum(I.*(r*ones(1,T))))...%网损费用
          + 500*sum(sum(Pg))...%主网购电成本
          + 500*sum(sum(p_wt));%弃风成本
          + 500*sum(sum(p_pv));%弃光成本
toc%建模时间
%% 5.设求解器
ops=sdpsettings('verbose', 1, 'solver', 'gurobi');
sol=optimize(C,objective,ops);
toc%求解时间
objective=100*value(objective)

clear branch Constraints dnstream upstream  i mpc nb nl ops Pgmax Qgmax Pin Qin r x Vmax Vmin T 
%% 6.输出AMPL模型
% saveampl(Constraints,objective,'mymodel');

%% 7.分析错误标志
if sol.problem == 0
    disp('succcessful solved');
else
    disp('error');
    yalmiperror(sol.problem)
end
%% 8.保存运行结果到XLSX
% Doc_name='OPF_data';
% V = value(V);I = value(I);P = value(P);Q = value(Q);
% Pg = value(Pg);Qg = value(Qg);
% xlswrite(Doc_name,V,'V');
% xlswrite(Doc_name,I,'I');
% xlswrite(Doc_name,P,'P');
% xlswrite(Doc_name,Q,'Q');
% xlswrite(Doc_name,Pg,'Pg');
% xlswrite(Doc_name,Qg,'Qg');



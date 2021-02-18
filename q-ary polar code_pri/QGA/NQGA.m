function [BestP,best_iter] = NQGA(N,K,q,ebno1,ebno2,M,GA_itermax,max_runs,alpha)
% addpath('../Decoding_Index/')
% addpath('../GA/')
% addpath('../chartsave/')
% addpath('../')
load(['q',num2str(q),'.mat'],'addchart','mulchart','divchart');

% donot delete this code!
temp = addchart;
addchart = temp;
temp = mulchart;
mulchart = temp;
temp = divchart;
divchart = temp;

max_iter = 0;
max_err = 0;
PopNum=M;   %种群规模
ChomLength=N;  %染色体位长度

Iteration=GA_itermax;   %迭代次数
t=0;    %当前代数
BestFit_Num=zeros(Iteration,1);   %迭代每步的最优解

Q=InitPop(PopNum,ChomLength); %初始化染色体
P=Watch(Q,K);
Fitness = zeros(M,1);
for i = 1:M
[~, ber1] = Simulation_given_construction(max_iter, max_err, max_runs, find(P(i,:)==1), ebno1, N, K,q,addchart,mulchart,divchart,alpha);
[~, ber2] = Simulation_given_construction(max_iter, max_err, max_runs, find(P(i,:)==1), ebno2, N, K,q,addchart,mulchart,divchart,alpha);
Fitness(i,1) = 1-ber1*ber2;
end
%Fitness = FitMeasure(P);  %测量其适应度
[BestFit,BestIndex]=max(Fitness);   %记录最佳个体状态极其适应度值
BestP=P(BestIndex,:);       %最优个体的二进制编码（状态）
best_iter = 1;
for t=1:1:Iteration
    disp(t);
    P=Watch(Q,K);
    Fitness = zeros(M,1);
    for i = 1:M
    [~, ber1] = Simulation_given_construction(max_iter, max_err, max_runs, find(P(i,:)==1), ebno1, N, K,q,addchart,mulchart,divchart,alpha);
    [~, ber2] = Simulation_given_construction(max_iter, max_err, max_runs, find(P(i,:)==1), ebno2, N, K,q,addchart,mulchart,divchart,alpha);
    Fitness(i,1) = 1-ber1*ber2;
    end
    %Fitness = FitMeasure(P);  %测量其适应度
    Q=QUpdate(Fitness, BestFit,P, BestP, Q);    %利用旋转门更新Q
    P=Watch(Q,K);
    %%%%%%更新最优B
    Fitness = zeros(M,1);
    for i = 1:M
    [~, ber1] = Simulation_given_construction(max_iter, max_err, max_runs, find(P(i,:)==1), ebno1, N, K,q,addchart,mulchart,divchart,alpha);
    [~, ber2] = Simulation_given_construction(max_iter, max_err, max_runs, find(P(i,:)==1), ebno2, N, K,q,addchart,mulchart,divchart,alpha);
    Fitness(i,1) = 1-ber1*ber2;
    end
    %Fitness = FitMeasure(P);  %测量其适应度
    [BestFit_temp,BestIndex]=max(Fitness);   %记录最佳个体状态极其适应度值    
    BestP_temp=P(BestIndex,:);       %最优个体的二进制编码（状态）   
    if BestFit_temp> BestFit
        BestFit=BestFit_temp;
        BestP=BestP_temp;
        best_iter = t;
    end        
    BestFit_Num(t)=BestFit;
end


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
PopNum=M;   %��Ⱥ��ģ
ChomLength=N;  %Ⱦɫ��λ����

Iteration=GA_itermax;   %��������
t=0;    %��ǰ����
BestFit_Num=zeros(Iteration,1);   %����ÿ�������Ž�

Q=InitPop(PopNum,ChomLength); %��ʼ��Ⱦɫ��
P=Watch(Q,K);
Fitness = zeros(M,1);
for i = 1:M
[~, ber1] = Simulation_given_construction(max_iter, max_err, max_runs, find(P(i,:)==1), ebno1, N, K,q,addchart,mulchart,divchart,alpha);
[~, ber2] = Simulation_given_construction(max_iter, max_err, max_runs, find(P(i,:)==1), ebno2, N, K,q,addchart,mulchart,divchart,alpha);
Fitness(i,1) = 1-ber1*ber2;
end
%Fitness = FitMeasure(P);  %��������Ӧ��
[BestFit,BestIndex]=max(Fitness);   %��¼��Ѹ���״̬������Ӧ��ֵ
BestP=P(BestIndex,:);       %���Ÿ���Ķ����Ʊ��루״̬��
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
    %Fitness = FitMeasure(P);  %��������Ӧ��
    Q=QUpdate(Fitness, BestFit,P, BestP, Q);    %������ת�Ÿ���Q
    P=Watch(Q,K);
    %%%%%%��������B
    Fitness = zeros(M,1);
    for i = 1:M
    [~, ber1] = Simulation_given_construction(max_iter, max_err, max_runs, find(P(i,:)==1), ebno1, N, K,q,addchart,mulchart,divchart,alpha);
    [~, ber2] = Simulation_given_construction(max_iter, max_err, max_runs, find(P(i,:)==1), ebno2, N, K,q,addchart,mulchart,divchart,alpha);
    Fitness(i,1) = 1-ber1*ber2;
    end
    %Fitness = FitMeasure(P);  %��������Ӧ��
    [BestFit_temp,BestIndex]=max(Fitness);   %��¼��Ѹ���״̬������Ӧ��ֵ    
    BestP_temp=P(BestIndex,:);       %���Ÿ���Ķ����Ʊ��루״̬��   
    if BestFit_temp> BestFit
        BestFit=BestFit_temp;
        BestP=BestP_temp;
        best_iter = t;
    end        
    BestFit_Num(t)=BestFit;
end


clear;
addpath('Decoding_Index/')
addpath('GA/')
addpath('chartsave/')
n = 2;
q = 2;
N = 2^n;
K = 2^(n - 1);
max_iter = 50;
max_err = 800;
max_runs = 1e5;
resolution = 1e2;
ebno_vec = 1 : 0.5 : 3.5;

% sampleM = GeneticA(N, K,q, 1,3,1000,2000,0.03,0.01,2);
% 
% [bler1, ber1] = Simulation_given_construction_for_main(max_iter, max_err, max_runs, resolution, find(sampleM(1,:)==1), ebno_vec, N, K,q,2);
% 
% [bler, ber] = Simulation(max_iter, max_err, max_runs, resolution, ebno_vec, N, K,q,2);

[bler, ber] = Simulation_for2to4(max_iter, max_err, max_runs, resolution, ebno_vec, N, K,q,1);

% [bler1, ber1] = Simulation(max_iter, max_err, max_runs, resolution, ebno_vec, N, K,2,1);

%info_esti = info_chosen_by_GenieSC(max_iter, max_err, 10,10000, resolution, [1.5,3.5], N, K,q);

%[bler2, ber2] = Simulation_given_construction_for_main(max_iter, max_err, max_runs, resolution, info_esti, ebno_vec, N, K,q);

% info = combntns(1:N,K);
% bler_t = zeros(length(bler),nchoosek(N,K));
% ber_t = zeros(length(bler),nchoosek(N,K));
% for i = 1:nchoosek(N,K)
%     [bler_tmp, ber_tmp] = Simulation_given_construction(max_iter, max_err, max_runs, info(i,:), ebno_vec, N, K,q);
%     bler_t(:,i) = bler_tmp;
%     ber_t(:,i) = ber_tmp;
% end





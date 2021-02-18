clear;
addpath('Decoding_Index/')
addpath('GA/')
addpath('chartsave/')
addpath('QGA/')
%addpath('D:/iTsinghua/Major/github_syncs/Encoding/PolarCpp/PolarCpp/x64/Release');

n = 6;
q = 4;
N = 2^n;
K = 2^(n - 1);
max_iter = 50;
max_err = 800;
max_runs = 1e5;
resolution = 1e2;
ebno_vec = 1 : 0.5 : 3.5;

% info5 = [12    14    15    16    20    22    23    24    25    26    27    28    29    30    31    32];
% prior5 = zeros(1,32);
% prior5(info5) = 1;

info6 = [16    24    27    28    29    30    31    32    39    40    42    43    44    45    46    47    48    50    51    52    53    54  55    56    57    58    59    60    61    62    63    64];
prior6 = zeros(1,64);
prior6(info6) = 1;

%[info,best_iter]=NQGA(32,16,2,1,3,100,100,1e3,1);
% sampleM = GeneticA(N, K,q, 1,3,1000,19000,0.03,0.01,2,0,0);
% % 
% [bler1, ber1] = Simulation_given_construction_for_main(max_iter, max_err, max_runs, resolution, find(sampleM(1,:)==1), ebno_vec, N, K,q,2);
% 
% sampleM = GeneticA(N, K,q, 1,3,1000,2000,0.03,0.01,2,1,prior6);
% 
% [bler2, ber2] = Simulation_given_construction_for_main(max_iter, max_err, max_runs, resolution, find(sampleM(1,:)==1), ebno_vec, N, K,q,2);
% 
%[bler, ber] = Simulation(max_iter, max_err, max_runs, resolution, ebno_vec, N, K,q,2);

%[bler1, ber1] = Simulation_GenieSC(max_iter, max_err, max_runs, resolution, ebno_vec, N, K,q,2);

%[bler, ber] = Simulation(max_iter, max_err, max_runs, resolution, ebno_vec, N, K,2,1);

info_esti = info_chosen_by_GenieSC(max_iter, max_err, 1,1e5, 8, [1,3], N, K,q,2);

[bler3, ber3] = Simulation_given_construction_for_main(max_iter, max_err, max_runs, resolution, info_esti, ebno_vec, N, K,q,2);

info_esti1 = info_chosen_by_GenieSC(max_iter, max_err, 1,1e5, 8, [1,3], N, N,q,2);
[bler4, ber4] = Simulation_given_construction_for_main(max_iter, max_err, max_runs, resolution, info_esti, ebno_vec, N, K,q,2);

% info = combntns(1:N,K);
% bler_t = zeros(length(bler),nchoosek(N,K));
% ber_t = zeros(length(bler),nchoosek(N,K));
% for i = 1:nchoosek(N,K)
%     [bler_tmp, ber_tmp] = Simulation_given_construction(max_iter, max_err, max_runs, info(i,:), ebno_vec, N, K,q);
%     bler_t(:,i) = bler_tmp;
%     ber_t(:,i) = ber_tmp;
% end





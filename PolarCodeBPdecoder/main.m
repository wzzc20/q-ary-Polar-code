clear;
addpath('Decoding_Index/')
addpath('GA/')
n = 6;
N = 2^n;
K = 2^(n - 1);
max_iter = 50;
max_err = 100;
max_runs = 1e5;
resolution = 1e5;
ebno_vec = 1 : 0.5 : 3.5;
[bler, ber] = Simulation(max_iter, max_err, max_runs, resolution, ebno_vec, N, K);

sampleM = GeneticA_pool_expand(N, K, 1.5,3,1000,2000,0.03,0.01);

[bler1, ber1] = Simulation_given_construction(max_iter, max_err, max_runs, find(sampleM(1,:)==1), ebno_vec, N, K);

% sampleM = GeneticA_pool_expand(N, K, 1.5,3,1000,3000,0.03,0.01);
% 
% [bler2, ber2] = Simulation_given_construction(max_iter, max_err, max_runs, find(sampleM(1,:)==1), ebno_vec, N, K);
% 
% sampleM = GeneticA_pool_expand(N, K, 1.5,3,2000,20000,0.03,0.01);
% 
% [bler3, ber3] = Simulation_given_construction(max_iter, max_err, max_runs, find(sampleM(1,:)==1), ebno_vec, N, K);

% sampleM = GeneticA_pool_expand(N, K, 1.5,3,1000,6000,0.03,0.01);
% 
% [bler2, ber2] = Simulation_given_construction(max_iter, max_err, max_runs, find(sampleM(1,:)==1), ebno_vec, N, K);
% 
% sampleM = GeneticA_pool_expand(N, K, 1.5,3,1000,6000,0.05,0.01);
% 
% [bler3, ber3] = Simulation_given_construction(max_iter, max_err, max_runs, find(sampleM(1,:)==1), ebno_vec, N, K);
% 
% sampleM = GeneticA_pool_expand(N, K, 1.5,3,1000,6000,0.01,0.01);
% 
% [bler4, ber4] = Simulation_given_construction(max_iter, max_err, max_runs, find(sampleM(1,:)==1), ebno_vec, N, K);
% 
% sampleM = GeneticA_pool_expand(N, K, 1.5,3,1000,6000,0.03,0.02);
% 
% [bler5, ber5] = Simulation_given_construction(max_iter, max_err, max_runs, find(sampleM(1,:)==1), ebno_vec, N, K);
% 
% sampleM = GeneticA_pool_expand(N, K, 1.5,3,1000,6000,0.03,0.005);
% 
% [bler6, ber6] = Simulation_given_construction(max_iter, max_err, max_runs, find(sampleM(1,:)==1), ebno_vec, N, K);
% 
% sampleM = GeneticA_pool_expand(N, K, 1.5,3,1000,6000,0.03,0.05);
% 
% [bler7, ber7] = Simulation_given_construction(max_iter, max_err, max_runs, find(sampleM(1,:)==1), ebno_vec, N, K);

% sampleM = GeneticA(N, K, 1,3,2000,1000,0.03,0.01);

% [bler1, ber1] = Simulation_given_construction(max_iter, max_err, max_runs, find(sampleM(1,:)==1), ebno_vec, N, K);
% 
% disp('finish 1');
% 
% sampleM = GeneticA(N, K, 1,3,2000,2000,0.03,0.01);
% 
% [bler2, ber2] = Simulation_given_construction(max_iter, max_err, max_runs, find(sampleM(1,:)==1), ebno_vec, N, K);
% 
% disp('finish 2');
% 
% sampleM = GeneticA(N, K, 1,3,3000,1000,0.03,0.01);
% 
% [bler3, ber3] = Simulation_given_construction(max_iter, max_err, max_runs, find(sampleM(1,:)==1), ebno_vec, N, K);
% 
% disp('finish 3');
% 
% sampleM = GeneticA(N, K, 1,3,1000,3000,0.03,0.01);
% 
% [bler4, ber4] = Simulation_given_construction(max_iter, max_err, max_runs, find(sampleM(1,:)==1), ebno_vec, N, K);
% 
% disp('finish 4');
% 
% sampleM = GeneticA(N, K, 1,3,4000,3000,0.03,0.01);
% 
% [bler5, ber5] = Simulation_given_construction(max_iter, max_err, max_runs, find(sampleM(1,:)==1), ebno_vec, N, K);
% 
% disp('finish 5');



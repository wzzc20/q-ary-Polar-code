function [bler, ber] = Simulation_for2to4(max_iter, max_err, max_runs, resolution, ebno_vec, N, K,q,alpha)
R = K/N;
n = log2(N);
num_block_err_bp = zeros(length(ebno_vec), 3);
num_bit_err_bp = zeros(length(ebno_vec), 3);
num_runs = zeros(length(ebno_vec), 1);

load(['q',num2str(2),'.mat']);
addchart2 = addchart;
mulchart2 = mulchart;
divchart2 = divchart;
load(['q',num2str(4),'.mat']);
addchart4 = addchart;
mulchart4 = mulchart;
divchart4 = divchart;

%Indices for enc/decoding
[M_up, M_down] = index_Matrix(N);
lambda_offset = 2.^(0 : n);
llr_layer_vec = get_llr_layer(N);
bit_layer_vec = get_bit_layer(N);

lambda_offset4 = 2.^(0 : n-1);
llr_layer_vec4 = get_llr_layer(N/2);
bit_layer_vec4 = get_bit_layer(N/2);

%code construction
design_snr = 2.5;
sigma_cc = 1/sqrt(2 * R) * 10^(-design_snr/20);
[channels, ~] = GA(sigma_cc, N);
[~, channel_ordered] = sort(channels, 'descend');
info_bits = sort(channel_ordered(1 : K), 'ascend');
frozen_bits = ones(N , 1);
frozen_bits(info_bits) = 0;

tic
for i_run = 1 : max_runs 
    if mod(i_run, ceil(max_runs/resolution)) == 1
        disp(['Sim iteration running = ', num2str(i_run)]);
        disp(['N = ' num2str(N) ' K = ' num2str(K) ' GA construction SNR = ' num2str(design_snr) 'dB' ' Max Iter Number = ' num2str(max_iter)])
        disp('BP BLER')
        disp(num2str([ebno_vec' num_block_err_bp./num_runs]));
        disp('BP BER')
        disp(num2str([ebno_vec' num_bit_err_bp./num_runs/K]));
        disp(' ')
    end
    u = zeros(N, 1);
    info = floor(q*rand(K, 1)) ;
%     u = zeros(N, 1);
%     info = floor(2*rand(K, 1));
    u(info_bits) = info;
    x = polar_encoder(u, q,lambda_offset, llr_layer_vec,addchart2,mulchart2,divchart2,alpha);
    %bpsk = 1 - 2 * x;
    bpsk = 1-2*(double(dec2bin(x))-48);
    noise = randn(N, log2(q)); 
    for i_ebno = 1 : length(ebno_vec) 
        if all(num_block_err_bp(i_ebno,:) > max_err)
            continue;
        end
        num_runs(i_ebno) = num_runs(i_ebno) + 1;
        sigma = 1 / sqrt(2 * R) * 10^(-ebno_vec(i_ebno)/20);
        y = bpsk + sigma * noise;
        pro = zeros(size(y,1),q);
        for pro_iter = 1:q
            one_count = dec2bin(pro_iter-1,log2(q))=='1';
            one_count = 1-2*one_count;
            one_count = repmat(one_count,size(y,1),1);
            pro(:,pro_iter) = power(1/(sigma*sqrt(2*pi)),log2(q))*exp(-sum((y-one_count).^2,2))/(2*sigma^2);
        end
        %llr = 2/sigma^2 * y;
        %[info_esti_bp, ~, ~, ~] = BP_Decoder_LLR(info_bits, frozen_bits, llr, max_iter, M_up, M_down);
        info_esti_SC = SC_decoder(q,pro, K, frozen_bits, lambda_offset, llr_layer_vec, bit_layer_vec,addchart2,mulchart2,divchart2,alpha);
        info_esti_SC1 = SC_2to4_decoder(q,pro, K, frozen_bits, lambda_offset4, llr_layer_vec4, bit_layer_vec4,addchart4,mulchart4,divchart4,alpha);
        info_esti_SC2 = ML(q,y, K, frozen_bits, lambda_offset, llr_layer_vec, bit_layer_vec,addchart2,mulchart2,divchart2,alpha);
        if any(info_esti_SC ~= info)
            num_block_err_bp(i_ebno,1) =  num_block_err_bp(i_ebno,1) + 1;
            num_bit_err_bp(i_ebno,1) = num_bit_err_bp(i_ebno,1) + sum(info ~= info_esti_SC);
        end
        if any(info_esti_SC1 ~= info)
            num_block_err_bp(i_ebno,2) =  num_block_err_bp(i_ebno,2) + 1;
            num_bit_err_bp(i_ebno,2) = num_bit_err_bp(i_ebno,2) + sum(info ~= info_esti_SC1);
        end
        if any(info_esti_SC2 ~= info)
            num_block_err_bp(i_ebno,3) =  num_block_err_bp(i_ebno,3) + 1;
            num_bit_err_bp(i_ebno,3) = num_bit_err_bp(i_ebno,3) + sum(info ~= info_esti_SC2);
        end

    end
end
toc
bler = num_block_err_bp./num_runs;
ber = num_bit_err_bp./num_runs/K;
end

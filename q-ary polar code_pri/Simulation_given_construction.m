function [bler, ber] = Simulation_given_construction(max_iter, max_err, max_runs, info_bits, ebno_vec, N, K,q,addchart,mulchart,divchart,alpha)
R = K/N;
n = log2(N);
num_block_err_bp = zeros(length(ebno_vec), max_runs);
num_bit_err_bp = zeros(length(ebno_vec), max_runs);
%num_runs = zeros(length(ebno_vec), 1);

%Indices for enc/decoding
[M_up, M_down] = index_Matrix(N);
lambda_offset = 2.^(0 : n);
llr_layer_vec = get_llr_layer(N);
bit_layer_vec = get_bit_layer(N);

decoder_config.partially_frozen = false;
decoder_config.is_qary = true;
decoder_config.is_LLR = false;
decoder_config.update_decoder = true;
decoder_config.is_list = true;
decoder_config.L = 16;

%code construction
frozen_bits = ones(N , 1);
frozen_bits(info_bits) = 0;

for i_ebno = 1 : length(ebno_vec)
    parfor i_run = 1 : max_runs 
    u = zeros(N, 1);
    info = floor(q*rand(K, 1)) ;
%     u = zeros(N, 1);
%     info = floor(2*rand(K, 1));
    u(info_bits) = info;
    x = polar_encoder(u, q,lambda_offset, llr_layer_vec,addchart,mulchart,divchart,alpha);
    %bpsk = 1 - 2 * x;
    bpsk = 1-2*(double(dec2bin(x))-48);
    noise = randn(N, log2(q)); 
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
        pro = pro ./ sum(pro,2);
        pro = pro.';
        info_esti_SC = Qary_SC_Decoder(pro, N, log2(q), frozen_bits.', alpha, decoder_config);
        t = info_esti_SC.';
        info_esti_SC = endian_reverse(t,q);
        % info [3 1]' -> info_split [1 1 1 0]'  reversed by zja
        info_split = qary2binary(info,q);
        if any(info_esti_SC ~= info_split)
            num_block_err_bp(i_ebno,i_run) =  num_block_err_bp(i_ebno,i_run) + 1;
            num_bit_err_bp(i_ebno,i_run) = num_bit_err_bp(i_ebno,i_run) + sum(info_split ~= info_esti_SC);
        end

    end
end
bler = sum(num_block_err_bp,2)/max_runs;
ber = sum(num_bit_err_bp,2)/max_runs/K;
end
function info_esti = info_chosen_by_GenieSC(max_iter, max_err, max_runs,max_info_num, resolution, ebno_vec, N, K,q,alpha)
R = K/N;
n = log2(N);
num_bit_err_bp = zeros(length(ebno_vec), N);
num_runs = zeros(length(ebno_vec), 1);

load(['q',num2str(q),'.mat']);

%Indices for enc/decoding
[M_up, M_down] = index_Matrix(N);
lambda_offset = 2.^(0 : n);
llr_layer_vec = get_llr_layer(N);
bit_layer_vec = get_bit_layer(N);

tic
for info_num = 1:max_info_num
    if(mod(info_num,200)==0)
        disp(info_num);
    end
    info_bits = randperm(N);
    info_bits = info_bits(1:K);
    frozen_bits = ones(N , 1);
    frozen_bits(info_bits) = 0;

    for i_run = 1 : max_runs 
%         if mod(i_run, ceil(max_runs/resolution)) == 1
%             disp(['Sim iteration running = ', num2str(i_run)]);
%             disp(['N = ' num2str(N) ' K = ' num2str(K) ' GA construction SNR = ' num2str(design_snr) 'dB' ' Max Iter Number = ' num2str(max_iter)])
%             disp('BP BLER')
%             disp(num2str([ebno_vec' num_block_err_bp./num_runs]));
%             disp(' ')
%         end
        u = zeros(N, 1);
        info = floor(q*rand(K, 1)) ;
    %     u = zeros(N, 1);
    %     info = floor(2*rand(K, 1));
        u(info_bits) = info;
        x = polar_encoder(u, q,lambda_offset, llr_layer_vec,addchart,mulchart,divchart,alpha);
        %bpsk = 1 - 2 * x;
        bpsk = 1-2*(double(dec2bin(x))-48);
        noise = randn(N, log2(q)); 
        for i_ebno = 1 : length(ebno_vec) 
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
            info_esti_SC = SC_decoder(1,u,q,pro,K, frozen_bits, lambda_offset, llr_layer_vec, bit_layer_vec,addchart,mulchart,divchart,alpha);
            for u_iter = 1:K
                if info_esti_SC(u_iter) ~= info(u_iter)
                    num_bit_err_bp(i_ebno,info_bits(u_iter)) = num_bit_err_bp(i_ebno,info_bits(u_iter)) + 1;
                end
            end

        end
    end
end
toc
ber = prod(num_bit_err_bp);
[~,index] = sort(ber);
info_esti = index(1:K);
info_esti = sort(info_esti);
end

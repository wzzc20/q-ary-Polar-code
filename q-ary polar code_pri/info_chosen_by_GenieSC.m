function info_esti = info_chosen_by_GenieSC(max_iter, max_err, max_runs,max_info_num, num_of_CPUs, ebno_vec, N, K,q,alpha)
R = K/N;
n = log2(N);
num_bit_err_bp = cell(num_of_CPUs,1);
for num = 1:num_of_CPUs
    num_bit_err_bp{num} = zeros(length(ebno_vec), N);
end
%num_bit_err_bp = zeros(num_of_CPUs,length(ebno_vec), N);

load(['q',num2str(q),'.mat']);
% donot delete this code!
temp = addchart;
addchart = temp;
temp = mulchart;
mulchart = temp;
temp = divchart;
divchart = temp;

%Indices for enc/decoding
[M_up, M_down] = index_Matrix(N);
lambda_offset = 2.^(0 : n);
llr_layer_vec = get_llr_layer(N);
bit_layer_vec = get_bit_layer(N);

% Decoder configuration.
decoder_config.partially_frozen = false;
decoder_config.is_qary = true;
decoder_config.is_LLR = false;
decoder_config.update_decoder = true;
decoder_config.is_Genie = true;



tic
parfor cpu_num = 1:num_of_CPUs
    for info_num = 1:floor(max_info_num/num_of_CPUs)
    %     if(mod(info_num,200)==0)
    %         disp(info_num);
    %     end
        info_bits = randperm(N);
        info_bits = info_bits(1:K);
        info_bits = sort(info_bits);
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
                %num_runs(i_ebno) = num_runs(i_ebno) + 1;
                sigma = 1 / sqrt(2 * R) * 10^(-ebno_vec(i_ebno)/20);
                y = bpsk + sigma * noise;
                pro = zeros(size(y,1),q);
                for pro_iter = 1:q
                    one_count = dec2bin(pro_iter-1,log2(q))=='1';
                    one_count = 1-2*one_count;
                    one_count = repmat(one_count,size(y,1),1);
                    pro(:,pro_iter) = power(1/(sigma*sqrt(2*pi)),log2(q))*exp(-sum((y-one_count).^2,2))/(2*sigma^2);
                end
                pro = pro ./ sum(pro,2);
                pro = pro.';
                %llr = 2/sigma^2 * y;
                %[info_esti_bp, ~, ~, ~] = BP_Decoder_LLR(info_bits, frozen_bits, llr, max_iter, M_up, M_down);
                temp_decoder_config = decoder_config;
                temp_decoder_config.update_decoder = true;
                info_esti_SC = Qary_SC_Decoder(pro, N, log2(q), frozen_bits.', alpha, temp_decoder_config,u);
                t = info_esti_SC.';
                info_esti_SC = endian_reverse(t,q);
                % info [3 1]' -> info_split [1 1 1 0]'  reversed by zja
                info_split = qary2binary(info,q);
                for u_iter = 1:K
                    temp1 = (u_iter-1)*log2(q)+1;
                    temp2 = u_iter*log2(q);
                    num_bit_err_bp{cpu_num}(i_ebno,info_bits(u_iter)) = num_bit_err_bp{cpu_num}(i_ebno,info_bits(u_iter)) + sum(xor(info_esti_SC(temp1:temp2),info_split(temp1:temp2)));
                end

            end
        end
    end
end
toc

num_bit_err_bp_sum = zeros(length(ebno_vec), N);
for num = 1:num_of_CPUs
    num_bit_err_bp_sum = num_bit_err_bp_sum+num_bit_err_bp{num};
end
ber = prod(num_bit_err_bp_sum);
ber = ber(end:-1:1);
[~,index] = sort(ber);
index = N+1-index;
if(N==K)
info_esti = index(1:K/2);
else
    info_esti = index(1:K);
end
info_esti = sort(info_esti);
end

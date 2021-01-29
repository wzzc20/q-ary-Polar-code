function polar_info_esti = ML(q,y,K, frozen_bits, lambda_offset, llr_layer_vec, bit_layer_vec,addchart,mulchart,divchart,alpha)
N = size(y,1);%llr refers to channel LLR.
n = log2(N);
info = dec2bin(0:(2^K-1), K);
best_save = zeros(K,1);
distance = 1e10;
for  i = 1:2^K
    info1 = double(info(i,:))-48;
    u = zeros(N,1);
    u(frozen_bits==0) = info1;
    x = polar_encoder(u, q,lambda_offset, llr_layer_vec,addchart,mulchart,divchart,alpha);
    distance1 = sum(abs(1-2*(double(dec2bin(x))-48)-y).^2);
    if(distance1<distance)
        best_save = info1';
        distance = distance1;
    end
end
polar_info_esti = best_save;


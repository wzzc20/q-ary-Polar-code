function polar_info_esti = SC_2to4_decoder(q,pro2, K, frozen_bits, lambda_offset, llr_layer_vec, bit_layer_vec,addchart,mulchart,divchart,alpha)
q = 4;
N = size(pro2,1);%llr refers to channel LLR.
N_4 = N/2;  %N->N/2
n = log2(N_4);
P = zeros(N_4 - 1, q);%channel pro
C = zeros(N_4 - 1, 2);%C stores internal bit values
polar_info_esti = zeros(K, 1);
polar_u_F4_esti = zeros(N_4,1);
polar_u_F2_esti = zeros(N,1);
cnt_K = 1;

% pro2 from F2 to F4 (length from N to N/2)
pro = zeros(N_4,4);
for i = 1:N_4
    pro(i,1) = pro2(2*i-1,1)*pro2(2*i,1);  %00
    pro(i,2) = pro2(2*i-1,1)*pro2(2*i,2);  %01
    pro(i,3) = pro2(2*i-1,2)*pro2(2*i,1);  %10
    pro(i,4) = pro2(2*i-1,2)*pro2(2*i,2);  %11
end

for phi = 0 : N_4 - 1
    switch phi
        case 0%for decoding u_1
            index_1 = lambda_offset(n);
            for beta = 0 : index_1 - 1%use llr vector
                %P(beta + index_1) =  sign(llr(beta + 1)) * sign(llr(beta + 1 + index_1)) * min(abs(llr(beta + 1)), abs(llr(beta + 1 + index_1)));
                P(beta + index_1,:) = upcalculate(pro(beta+1,:),pro(beta+1+index_1,:),q,addchart,mulchart,alpha);
            end
            for i_layer = n - 2 : -1 : 0%use P vector
                index_1 = lambda_offset(i_layer + 1);
                index_2 = lambda_offset(i_layer + 2);
                for beta = index_1 : index_2 - 1
                    %P(beta) =  sign(P(beta + index_1)) * sign(P(beta + index_2)) * min(abs(P(beta + index_1)), abs(P(beta + index_2)));
                    P(beta,:) = upcalculate(P(beta + index_1,:),P(beta + index_2,:),q,addchart,mulchart,alpha);
                end
            end
        case N_4/2%for deocding u_{N/2 + 1}
            index_1 = lambda_offset(n);
            for beta = 0 : index_1 - 1%use llr vector. g function.
                %P(beta + index_1) = (1 - 2 * C(beta + index_1, 1)) * llr(beta + 1) + llr(beta + 1 + index_1);
                P(beta + index_1,:) = downcalculate(pro(beta + 1,:),pro(beta + 1 + index_1,:),C(beta + index_1, 1),q,addchart,mulchart,alpha);
            end
            for i_layer = n - 2 : -1 : 0%use P vector. f function
                index_1 = lambda_offset(i_layer + 1);
                index_2 = lambda_offset(i_layer + 2);
                for beta = index_1 : index_2 - 1
                    %P(beta) =  sign(P(beta + index_1)) * sign(P(beta + index_2)) * min(abs(P(beta + index_1)), abs(P(beta + index_2)));
                    P(beta,:) = upcalculate(P(beta + index_1,:),P(beta + index_2,:),q,addchart,mulchart,alpha);
                end
            end
        otherwise
            llr_layer = llr_layer_vec(phi + 1);
            index_1 = lambda_offset(llr_layer + 1);
            index_2 = lambda_offset(llr_layer + 2);
            for beta = index_1 : index_2 - 1%g function is first implemented.
                %P(beta) = (1 - 2 * C(beta, 1)) * P(beta + index_1) + P(beta + index_2);
                P(beta ,:) = downcalculate(P(beta + index_1,:),P(beta + index_2,:),C(beta, 1),q,addchart,mulchart,alpha);
            end
            for i_layer = llr_layer - 1 : -1 : 0%then f function is implemented.
                index_1 = lambda_offset(i_layer + 1);
                index_2 = lambda_offset(i_layer + 2);
                for beta = index_1 : index_2 - 1
                    %P(beta) =  sign(P(beta + index_1)) * sign(P(beta + index_2)) * min(abs(P(beta + index_1)), abs(P(beta + index_2)));
                    P(beta,:) = upcalculate(P(beta + index_1,:),P(beta + index_2,:),q,addchart,mulchart,alpha);
                end
            end
    end
    phi_mod_2 = mod(phi, 2);
%     if frozen_bits(phi + 1) == 1%frozen bit
%         C(1, 1 + phi_mod_2) = 0;
%     else%information bit
%         [~,index] = max(P(1,:));
%         C(1, 1 + phi_mod_2) = index-1;%store internal bit values
%         polar_info_esti(cnt_K) = index-1;
%         cnt_K = cnt_K + 1;
%     end
    
    % the traspose at the first stage
    % 00->00 01->11 10->10 11->01
    temp = P(1,2);
    P(1,2) = P(1,4);
    P(1,4) = temp;
    if frozen_bits(2*phi+1)==1
        P(1,3) = 0;
        P(1,4) = 0;
    end
    if frozen_bits(2*phi+2)==1
        P(1,2) = 0;
        P(1,4) = 0;        
    end

    [~,index] = max(P(1,:));
    if(index-1==0)
        C(1, 1 + phi_mod_2) = 0;%store internal bit values
    elseif(index-1==1)
        C(1, 1 + phi_mod_2) = 3;
    elseif(index-1==2)
        C(1, 1 + phi_mod_2) = 2;
    elseif(index-1==3)
        C(1, 1 + phi_mod_2) = 1;
    end
    polar_u_F4_esti(cnt_K) = index-1;
    if index-1 ==0
        polar_u_F2_esti(2*cnt_K-1) = 0;
        polar_u_F2_esti(2*cnt_K) = 0;
    elseif index-1 ==1
        polar_u_F2_esti(2*cnt_K-1) = 0;
        polar_u_F2_esti(2*cnt_K) = 1;
    elseif index-1 ==2
        polar_u_F2_esti(2*cnt_K-1) = 1;
        polar_u_F2_esti(2*cnt_K) = 0;
    else
        polar_u_F2_esti(2*cnt_K-1) = 1;
        polar_u_F2_esti(2*cnt_K) = 1;
    end
    cnt_K = cnt_K + 1;    
    
    if phi_mod_2  == 1 && phi ~= N_4 - 1
        bit_layer = bit_layer_vec(phi + 1);
        for i_layer = 0 : bit_layer - 1%give values to the 2nd column of C
            index_1 = lambda_offset(i_layer + 1);
            index_2 = lambda_offset(i_layer + 2);
            for beta = index_1 : index_2 - 1
                C(beta + index_1, 2) = addchart(C(beta, 1) +1, mulchart(C(beta, 2)+1,alpha+1)+1);
                C(beta + index_2, 2) = C(beta, 2);
            end
        end
        index_1 = lambda_offset(bit_layer + 1);
        index_2 = lambda_offset(bit_layer + 2);
        for beta = index_1 : index_2 - 1%give values to the 1st column of C
            C(beta + index_1, 1) = addchart(C(beta, 1)+1, mulchart(C(beta, 2)+1,alpha+1)+1);
            C(beta + index_2, 1) = C(beta, 2);
        end
    end
end
polar_info_esti = polar_u_F2_esti(frozen_bits==0);
end

function pro = upcalculate(pro1,pro2,q,addchart,mulchart,alpha)
    pro = zeros(1,q);
    %alpha*(0:q-1)
    alpha_mul = mulchart(alpha+1,:); 
    for u1_id = 1:q
        qid = zeros(1,q);
        for i=1:q
            qid(i) = addchart(u1_id,alpha_mul(i)+1);
        end
        index = qid+1;
        pro(1,u1_id) = sum(pro1(index).*pro2)/q;
    end
end

function pro = downcalculate(pro1,pro2,u1,q,addchart,mulchart,alpha)
    pro = zeros(1,q);
    for u2_id = 1:q
        u1lambda = addchart(u1+1,mulchart(u2_id,alpha+1)+1);
        u1lambdaid = u1lambda+1;
        pro(1,u2_id) = pro1(u1lambdaid)*pro2(u2_id)/q;
    end
end
function polar_info_esti = SC_decoder(q,llr, K, frozen_bits, lambda_offset, llr_layer_vec, bit_layer_vec)
N = size(llr,1);%llr refers to channel LLR.
n = log2(N);
P = zeros(N - 1, q);%channel llr is not include in P.
C = gf(zeros(N - 1, 2),log2(q));%C stores internal bit values
polar_info_esti = gf(zeros(K, 1),log2(q));
cnt_K = 1;
for phi = 0 : N - 1
    switch phi
        case 0%for decoding u_1
            index_1 = lambda_offset(n);
            for beta = 0 : index_1 - 1%use llr vector
                %P(beta + index_1) =  sign(llr(beta + 1)) * sign(llr(beta + 1 + index_1)) * min(abs(llr(beta + 1)), abs(llr(beta + 1 + index_1)));
                P(beta + index_1,:) = upcalculate(llr(beta+1,:),llr(beta+1+index_1,:),q);
            end
            for i_layer = n - 2 : -1 : 0%use P vector
                index_1 = lambda_offset(i_layer + 1);
                index_2 = lambda_offset(i_layer + 2);
                for beta = index_1 : index_2 - 1
                    %P(beta) =  sign(P(beta + index_1)) * sign(P(beta + index_2)) * min(abs(P(beta + index_1)), abs(P(beta + index_2)));
                    P(beta,:) = upcalculate(P(beta + index_1,:),P(beta + index_2,:),q);
                end
            end
        case N/2%for deocding u_{N/2 + 1}
            index_1 = lambda_offset(n);
            for beta = 0 : index_1 - 1%use llr vector. g function.
                %P(beta + index_1) = (1 - 2 * C(beta + index_1, 1)) * llr(beta + 1) + llr(beta + 1 + index_1);
                P(beta + index_1,:) = downcalculate(llr(beta + 1,:),llr(beta + 1 + index_1,:),C(beta + index_1, 1),q);
            end
            for i_layer = n - 2 : -1 : 0%use P vector. f function
                index_1 = lambda_offset(i_layer + 1);
                index_2 = lambda_offset(i_layer + 2);
                for beta = index_1 : index_2 - 1
                    %P(beta) =  sign(P(beta + index_1)) * sign(P(beta + index_2)) * min(abs(P(beta + index_1)), abs(P(beta + index_2)));
                    P(beta,:) = upcalculate(P(beta + index_1,:),P(beta + index_2,:),q);
                end
            end
        otherwise
            llr_layer = llr_layer_vec(phi + 1);
            index_1 = lambda_offset(llr_layer + 1);
            index_2 = lambda_offset(llr_layer + 2);
            for beta = index_1 : index_2 - 1%g function is first implemented.
                %P(beta) = (1 - 2 * C(beta, 1)) * P(beta + index_1) + P(beta + index_2);
                P(beta ,:) = downcalculate(P(beta + index_1,:),P(beta + index_2,:),C(beta, 1),q);
            end
            for i_layer = llr_layer - 1 : -1 : 0%then f function is implemented.
                index_1 = lambda_offset(i_layer + 1);
                index_2 = lambda_offset(i_layer + 2);
                for beta = index_1 : index_2 - 1
                    %P(beta) =  sign(P(beta + index_1)) * sign(P(beta + index_2)) * min(abs(P(beta + index_1)), abs(P(beta + index_2)));
                    P(beta,:) = upcalculate(P(beta + index_1,:),P(beta + index_2,:),q);
                end
            end
    end
    phi_mod_2 = mod(phi, 2);
    if frozen_bits(phi + 1) == 1%frozen bit
        C(1, 1 + phi_mod_2) = gf(0,log2(q));
    else%information bit
        [~,index] = min(P(1,:));
        C(1, 1 + phi_mod_2) = gf(index-1,log2(q));%store internal bit values
        polar_info_esti(cnt_K) = gf(index-1,log2(q));
        cnt_K = cnt_K + 1;
    end
    if phi_mod_2  == 1 && phi ~= N - 1
        bit_layer = bit_layer_vec(phi + 1);
        for i_layer = 0 : bit_layer - 1%give values to the 2nd column of C
            index_1 = lambda_offset(i_layer + 1);
            index_2 = lambda_offset(i_layer + 2);
            for beta = index_1 : index_2 - 1
                C(beta + index_1, 2) = C(beta, 1) + C(beta, 2);
                C(beta + index_2, 2) = C(beta, 2);
            end
        end
        index_1 = lambda_offset(bit_layer + 1);
        index_2 = lambda_offset(bit_layer + 2);
        for beta = index_1 : index_2 - 1%give values to the 1st column of C
            C(beta + index_1, 1) = C(beta, 1) + C(beta, 2);
            C(beta + index_2, 1) = C(beta, 2);
        end
    end
end
end

function llr = upcalculate(llr1,llr2,q)
    llr = zeros(1,q);
    for symid = 1:q
        numerator = sum(exp(-llr1-llr2));
        qid = gf(0:q-1,log2(q))+gf(symid-1,log2(q));
        index = double(qid.x)+1;
        denominator = sum(exp(-llr1(index)-llr2));
        llr(1,symid) = log(numerator/denominator);
    end
end

function llr = downcalculate(llr1,llr2,u1,q)
    llr = zeros(1,q);
    for symid = 1:q
        u1id = double(u1.x)+1;
        u1lambda = u1+(symid-1);
        u1lambdaid = double(u1lambda.x)+1;
        llr(1,symid) = -llr1(u1id)+llr1(u1lambdaid)+llr2(symid);
    end
end
q  =4;
llr1 = [0,0.5,-1,2];
llr2 = [0,-4,3,4];
    llr = zeros(1,q);
    for symid = 1:q
        numerator = sum(exp(-llr1-llr2));
        qid = gf(0:q-1,log2(q))+gf(symid-1,log2(q));
        index = double(qid.x)+1;
        denominator = sum(exp(-llr1(index)-llr2));
        llr(1,symid) = log(numerator/denominator);
    end
    llr
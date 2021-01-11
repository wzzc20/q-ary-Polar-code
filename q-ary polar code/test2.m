q  =4;
llr1 = [0,0.5,-1,2];
llr2 = [0,-4,3,4];
u1 = gf(0,2);
llr = zeros(1,q);
    for symid = 1:q
        u1id = double(u1.x)+1;
        u1lambda = u1+(symid-1);
        u1lambdaid = double(u1lambda.x)+1;
        llr(1,symid) = -llr1(u1id)+llr1(u1lambdaid)+llr2(symid);
    end
    llr
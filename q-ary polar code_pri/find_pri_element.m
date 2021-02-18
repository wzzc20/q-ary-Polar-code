addpath('chartsave/')

q = 16;
load(['q',num2str(q),'.mat']);
for i = 2:q-1
    mul = i;
    for iter = 1:q
        if(mul==1)
            disp(['i=' num2str(i) '  rank=' num2str(iter)]);
            break;
        end
        mul = mulchart(i+1,mul+1);
    end
end
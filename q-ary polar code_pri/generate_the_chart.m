q = 16;
%%¼Ó·¨
addchart = zeros(q,q);
for i = 1:q
    for j = 1:q
        temp = gf(i-1,log2(q))+gf(j-1,log2(q));
        addchart(i,j) = double(temp.x);
    end
end

mulchart = zeros(q,q);
for i = 1:q
    for j = 1:q
        temp = gf(i-1,log2(q))*gf(j-1,log2(q));
        mulchart(i,j) = double(temp.x);
    end
end

divchart = zeros(q,q);
for i = 1:q
    for j = 1:q
        if(j~=1)
        temp = gf(i-1,log2(q))/gf(j-1,log2(q));
        divchart(i,j) = double(temp.x);
        end
    end
end

save('chartsave\q16','addchart','mulchart','divchart');
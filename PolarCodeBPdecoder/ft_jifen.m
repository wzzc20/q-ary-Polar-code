x = normrnd(0,1,1,100000);
y = normrnd(0,1,100,100000);

z = x+y;
c = 0*x;

for  i =1:100000
    temp = z(:,i);
    temp1 = abs(temp);
    index = find(temp1==min(temp1));
    c(i) = temp(index);
end
varc = var(c);
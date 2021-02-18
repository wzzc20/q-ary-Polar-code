function [Fitness] = FitMeasure(P)
PopNum=size(P,1);
ChomLength=size(P,2);
Fitness=zeros(PopNum,1);
Incode=zeros(PopNum,1);
for i=1:1:PopNum
    sum=0;
    for j=1:1:ChomLength
        sum=sum+P(i,j)*2^(4-j);       
    end
    Incode(i)=sum;
end
%º∆À„  ”¶∂»
	for i = 1:1:PopNum
        Fitness(i) = exp(-0.001*Incode(i))*(cos(0.8*Incode(i)))^2;
	end
end


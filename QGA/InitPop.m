function InitQ = InitPop(PopNum,ChomLength )
InitQ=zeros(PopNum,ChomLength,2);  
for i=1:1:PopNum
    for j=1:1:ChomLength
        InitQ(i,j,1)=-1/sqrt(2);
        InitQ(i,j,2)=-1/sqrt(2);
    end
end
end


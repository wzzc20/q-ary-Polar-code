function [ P ] = Watch(Q,K)
PopNum=size(Q,1);
ChomLength=size(Q,2);
P=zeros(PopNum,ChomLength);
for i=1:1:PopNum
    for j=1:1:ChomLength
        if Q(i,j,1)^2<rand()
           P(i,j)=1;
        else
           P(i,j)=0;
        end
    end
    if(sum(P(i,:))<K)
        index = find(P(i,:)==0);
        alpha_in_index = Q(i,index,1);
        [~,index_alphasmall] = sort(alpha_in_index);
        index = index(index_alphasmall);
        P(i,index(1:K-sum(P(i,:))))=1;
    elseif(sum(P(i,:))>K)
        index = find(P(i,:)==1);
        beta_in_index = Q(i,index,2);
        [~,index_betasmall] = sort(beta_in_index);
        index = index(index_betasmall);
        P(i,index(1:sum(P(i,:))-K))=0;        
    end
end
end


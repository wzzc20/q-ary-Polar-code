function updated_Q = QUpdate(Fitness, BestFit,P, BestP, Q)

PopNum=size(Fitness,1);   %种群规模
ChomLength=size(Q,2);  %染色体位长度
updated_Q=zeros(PopNum,ChomLength,2);
delta=0.05*pi;

for i=1:1:PopNum
   for j=1:1:ChomLength
       alpha= Q(i,j,1);
       beta= Q(i,j,2);
       s=0;
       if (~Fitness(i)>BestFit && P(i,j)==0 && BestP(j)==1 && alpha*beta>0 ) || ...
          (Fitness(i)>BestFit && P(i,j)==1 && BestP(j)==0 && alpha*beta>0 )  || ...
          (Fitness(i)>BestFit && P(i,j)==0 && BestP(j)==1 && alpha*beta<0 )  || ...
          (~Fitness(i)>BestFit && P(i,j)==1 && BestP(j)==0 && alpha*beta<0 )
           s=1;    
       end    
       if (~Fitness(i)>BestFit && P(i,j)==0 && BestP(j)==1 && alpha*beta<0 ) || ...
          (Fitness(i)>BestFit && P(i,j)==0 && BestP(j)==1 && alpha*beta>0 )  || ...
          (~Fitness(i)>BestFit && P(i,j)==1 && BestP(j)==0 && alpha*beta>0 )  || ...
          (Fitness(i)>BestFit && P(i,j)==1 && BestP(j)==0 && alpha*beta<0 )
           s=-1;
       end 
%        if (~Fitness(i)>BestFit && P(i,j)==0 && BestP(j)==1 && beta==0 ) || ...
%           (Fitness(i)>BestFit && P(i,j)==0 && BestP(j)==1 && alpha==0 )  || ...
%           (~Fitness(i)>BestFit && P(i,j)==1 && BestP(j)==0 && alpha==0 )  || ...
%           (Fitness(i)>BestFit && P(i,j)==1 && BestP(j)==0 && beta==0 )
%            s=round(rand)*2-1;
%        end

       thete=delta*s;
       temp_alpha=cos(thete)*alpha-sin(thete)*beta;
       temp_beta=sin(thete)*alpha+cos(thete)*beta;
       updated_Q(i,j,1)=temp_alpha;
       updated_Q(i,j,2)=temp_beta;       
   end
end
end


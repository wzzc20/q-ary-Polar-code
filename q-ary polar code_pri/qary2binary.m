function [info_split] = qary2binary(info,q)
    len = length(info);
    info_split = zeros(len*log2(q),1);
    for i = 1:log2(q)
        info_split(log2(q)-i+1:log2(q):end-i+1,1)=mod(floor(info/2^(i-1)),2);
    end
end
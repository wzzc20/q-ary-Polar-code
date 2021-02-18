function [info_rev] = endian_reverse(info,q)
    info_rev = info*0;
    for i = 1:log2(q)
        info_rev(log2(q)-i+1:log2(q):end-i+1,1) = info(i:log2(q):end-log2(q)+i);
    end
end
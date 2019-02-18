function [BW_ind, peak_ind, peak_val, left_ind, right_ind] = BandWidth(spec,nsc)

[ind,val] = pickpeak(spec,nsc);
len = length(spec);
left_val = val;
left_ind = ind;
while left_val > val/sqrt(2)
    left_ind = left_ind-1;
    if left_ind == 0
        left_val=-inf;
    else
        left_val = spec(left_ind);
    end    
end
right_val = val;
right_ind = ind;
while right_val > val/sqrt(2)
    right_ind = right_ind+1;
    if right_ind == len+1;
        right_val = -inf;
    else
        right_val = spec(right_ind);
    end    
end

BW_ind = right_ind - left_ind;
peak_ind = ind;
peak_val = val;
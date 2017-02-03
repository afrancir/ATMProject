function topKey = UTopKey(U)
% function topKey = UTopKey(U)
    if isempty(U)
        topKey = [Inf,Inf];
    else 
        topKey = U(1,2:3);
    end 
end 
function topNodePriority = UTop(U)
% function topKey = UTopKey(U)
    if isempty(U)
        display('UTopKey error: U is empty')
        topNodePriority = [];
    else 
        topNodePriority = U(1,1);
    end 
end 
function U = UUpdate(U,node,newKey)
% function U = UUpdate(U,node,newKey)
index          = U(:,1)==node;
U(index,2:end) = newKey;
U              = sortrows(U,[2 3]);

end 
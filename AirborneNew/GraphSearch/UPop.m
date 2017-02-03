function [U,node] = UPop(U)
node = U(1,1);  % node with minimum priority

if size(U,1)==1 
    U = [];
else
    U = U(2:end,:); % queue with node removed
end 
% note: not really necessary but more elegant. Otherwise if size(U,1)==1
% doing U = U(2:end,:); would return a 3x1 empty matrix, which also works.

end 
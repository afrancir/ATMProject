function U = URemove(U,node)
% function U = URemove(U,node)
% description: removes node with ID #node from the priority queue U
index = find(U(:,1)==node);
% a check just in case misfunctioning were happening. U(:,1) is the list of
% nodes ID's present in the queue. Every node may contain as much as a
% single entrance so
if ~isscalar(index)
    display('Error: U has multiple entrances for the same node')
end

switch index 
    case 1
        U = U(2:end,:);
    case size(U,1)
        U = U(1:end-1,:);
    otherwise
        U = U([1:(index-1);(index+1):end],:);
end 


end %function
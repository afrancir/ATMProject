function U = UInsert(U,node,key)
% function U = UInsert(U,node,key)
% description: inserts 'node' to the vertices priority queue U with priority
% established by 'key'
U = [U;[node key]];
U = sortrows(U,[2 3]);
end 
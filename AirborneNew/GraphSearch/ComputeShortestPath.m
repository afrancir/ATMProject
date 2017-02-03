function [U,g,rhs] = ComputeShortestPath(G, U, edges, g, rhs, costs, sStart, sGoal,km)
% function [U,g,rhs] = ComputeShortestPath(G, U, edges, g, rhs, costs, sStart, sGoal)
%
% notice that G and "edges" contain the same info (just different class of
% data) so in a future one or the other should be chosen.

while (isKey1GreaterKey2(CalculateKey(sStart,g,rhs,sGoal,km),UTopKey(U))...
        || rhs(sStart)~=g(sStart) )
    kOld  = UTopKey(U);
    [U,u] = UPop(U);
    
    if isKey1GreaterKey2(CalculateKey(u,g,rhs,sGoal,km),kOld)
        U = UInsert(U,u,CalculateKey(u,g,rhs,sGoal,km)); %UInsert(U,node,key)
    elseif g(u)>rhs(u)
        g(u) = rhs(u);
        % for all "s" in Pred(u) --> UpdateVertex(s)
        preIDs = predecessors(G,u);
        for i = 1:numel(preIDs)
            s = preIDs(i);
            [U,rhs] = UpdateVertex(edges,U,rhs,g,costs,sGoal,s,km);
        end 
    else
        g(u) = Inf;
        % for all "s" in Pred(u) & "u" itself --> UpdateVertex(s)
        preIDsAndU = [u;predecessors(G,u)];
        for i = 1:numel(preIDs)
            s = preIDsAndU(i);
            [U,rhs] = UpdateVertex(edges,U,rhs,g,costs,sGoal,s,km);
        end 
    end 
    
end 

end %function



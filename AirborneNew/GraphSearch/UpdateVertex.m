function [U,rhs] = UpdateVertex(edges,U,rhs,g,costs,sGoal,u,km)
% function [U,rhs] = UpdateVertex(edges,U,rhs,g,costs,sGoal,u,km)

if u~=sGoal
    sucIndex = (edges(:,1) == u); %successors
    sucIDs   = edges(sucIndex,2);
    
    rhs(u)   = min(costs(sucIndex) + g(sucIDs));
end 

if any(u == U(:,1)) %if node u is in the prioirity queue
    U = URemove(U,u);
end 

if g(u)~=rhs(u)
    U = UInsert(U,u,CalculateKey(u,g,rhs,sGoal,km));
end 

end 
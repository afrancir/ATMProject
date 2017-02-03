%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%$$               Arnau Francí-Rodon  8-20-2016                      $$                                                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function shortestPath = Main_D_Lite_asexternalfunction(nodes, edges, ...
%     costs, sGoal, sStart)
%
% Description: This program computes the time-dependent 
% shortest path detour using D-Lite. 
% 
% INPUTS: 
%       - nodes: nx1 vector. 
%       - edges: mx2 matrix with all the directed connections.
%       - costs: mx1 vector with the corresponding costs of every.
%                connection.
%       - sGoal, sStart: scalar representing initial and final node ID.
% OUTPUTS:
%       - shortestPath: 1XN vector with the nodes' ID constituting the 
%                       shortest path.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function shortestPath = Main_D_Lite_asexternalfunction(nodes, edges, ...
    costs, sGoal, sStart)

Nnodes = size(nodes,1); 
Nedges = size(edges,1); 

if Nedges~= size(costs,1) %throw error
    display('Error: costs and edges dimensional mismatch')
end 
% it is confortable to use Matlab's graph package, so we will be working
% with that 
G = digraph(edges(:,1),edges(:,2));

%%% INITIALIZING STUFF $$(procedure Initialize{})$$.
U = [];
km = 0;
g   = Inf(Nnodes,1);
rhs = Inf(Nnodes,1);
rhs(sGoal) = 0;
U = UInsert(U,sGoal,CalculateKey(sGoal,g,rhs,sGoal,km));

shortestPath = sStart; %node ID's strain to follow for the shortest path 
%%% MAIN
sLast = sStart;
% {22} Initialize()
[U,g,rhs] = ComputeShortestPath(G, U, edges, g, rhs, costs, sStart, sGoal,km);

while(sStart ~= sGoal)
    %%% {25} if g(sStart) = Inf then there is no knwon path
    sucIndex = (edges(:,1) == sStart); %successors connections
    sucIDs   = edges(sucIndex,2); %sStart successors ID
    
    [~,sucIDsMin]  = min(costs(sucIndex) + g(sucIDs)); %moving sStart to next cheapest 
    sStart         = sucIDs(sucIDsMin);
    %%% move to sStart... how to interpret that
    shortestPath = [shortestPath, sStart]; %#ok<AGROW>
    %%%
    %%% SCAN GRAPH FOR CHANGED EDGE COSTS... FOR THE MOMENT WE TRY THE
    %%% ALGORITHM FOR STATIC COSTS BUT  AMETHOD TO DETECT THAT NEEDS TO BE
    %%% DEVELOPED WHENEVER COSTS HAVE CHANGED. 
    %%%
    costsChangedMarker = 0;
    if costsChangedMarker
        km = km + h(sLast, sStart);       %#ok<UNRCH>
        sLast = sStart;
        
        % for all changed costs 
        %   update costs 
        %   UpdateVertex(u)
        % end 
        [U,g,rhs] = ComputeShortestPath(G, U, edges, g, rhs, costs, sStart, sGoal,km);
    end 
end 

display(shortestPath)

end %function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%$$               Arnau Francí-Rodon  8-14-2016                        $$%                                                                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% description: This program computes the time-dependent 
% shortest path detour using D-Lite. 
% 
% INPUTS: Spatial-Temporal Graph (nodes and edge costs). Goal nodes are 
% every node after the avoidance polygon.  
% 
% Prioirity queue.  nx2 cell-matrix with the ID of the node in the first
% column and the pair of keys k = [k_1,k_2] in the second column. 

%%%DEFINING THE GRAPH
nodes = (1:10)'; 
edges = [1 2; ...
             1 3; ...
             2 3; ...
             3 4; ...
             3 5; ...
             3 6; ...
             6 7; ...
             6 8; ...
             7 9; ...
             8 7; ...
             8 10] ;
costs     = [1; ...
             1; ...
             4; ...
             8; ...
             1; ...
             5; ...
             2; ...
             3; ...
             1; ...
             1; ...
             1] ;
Nnodes = size(nodes,1); 
Nedges = size(edges,1); 
% the initial start node is #1 and the goal node is #9
sGoal  = 9;
sStart = 1;

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


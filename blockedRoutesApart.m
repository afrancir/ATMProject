%%% SCRIPT TO SELECT THOSE ROUTES WHICH ARE BLOCKED BY AVOIDANCE ZONES 

% the avoidance zones are allLatlonCH. Need to check if any of the routes
% is blocked and which are the avoidance zones (registrer them) causing
% such blocking. Therefore we need a loop for every

%%% INPUTS
    %%%allLatLonCH;
    %%%network;
%%% Variables for trying script. MADEUP VALUES.
load('blockedRoutesInput')
% loads ''LatLOnCh and 'network' variables


%%%%% function [blockedRoutes] = blockedRoutesApart(allLatlonCH,network)
%
%

nTotalRoutes = network(end,1); %blotched job to count it this way
affectedCounter = 0; %this way we can cut down the matrix size 
% (not all routes will be affected). 
%%% allLatlonCH;

field1 = 'ID';
field2 = 'timeSteps';
field3 = 'nodesBlocked';
field4 = 'polygonID'; %'polygonsIDForThisTimeStep';
value1 = (1:nTotalRoutes)';
value2 = cell(nTotalRoutes,1);
value3 = cell(nTotalRoutes,1);
value4 = cell(nTotalRoutes,1);  

blockedRoutes = struct(field1,value1,field2,value2,field3,value3,field4,value4);

% % %%% for example, it can be that 5 nodes are affected for route 15, which is
% % %%% the first to be afected
% % blockedRoutes(1).ID = 15;
% % blockedRoutes(1).timeSteps    = [3 4 5 12 13];
% % blockedRoutes(1).nodesBlocked = [14 15 16 22 23];
% % % as we can see there are 2 problematic zones, blocked by different
% % % avoidance zones
% % blockedRoutes(1).polygonID = [7 7 7 2 2];
% % %%%

% only considering first TimeStep for the moment. Probably we need to
% change that in a future
%%% LatLonCH = allLatlonCH{1};

for routeN=1:nTotalRoutes
    
    %%% determine which are affected
    routePoints = network(network(:,1)==routeN,2:3);
    xq = routePoints(:,2);
    yq = routePoints(:,1);
    
    %%%% read*, important, at the bottom
    IsFirstTime = true;
    for PolIter = 1:size(LatLonCH,1)        
        dummyValues = LatLonCH{PolIter};
        xv = dummyValues(:,2) ;
        yv = dummyValues(:,1);
        in = inpolygon(xq,yq,xv,yv);
        
        if numel(xq(in))~=0 % some points included; then note down those
            if IsFirstTime 
               affectedCounter = affectedCounter+1; 
            end 
            blockedRoutes(affectedCounter).ID = routeN;
            blockedRoutes(affectedCounter).timeSteps    = [blockedRoutes(affectedCounter).timeSteps; find(in)];
            blockedRoutes(affectedCounter).nodesBlocked = [blockedRoutes(affectedCounter).nodesBlocked; find(in)];
            blockedRoutes(affectedCounter).polygonID = [blockedRoutes(affectedCounter).polygonID; PolIter];            
% %             break; %if well implemented, there cant be more than 1 polygon blocking (otherwise they would have merged when taking CH)            
            IsFirstTime = false;
        end 
    end 
         
end
blockedRoutes = blockedRoutes(1:affectedCounter);


%%% I would say it works but I need to generate some CH such that they
%%% encompass some of the routes 

%%% Remember always lat - y and lon - x !!

%%%%* I am not considering timeSteps since I don't really know when they
%%%% start or they stop, I believe the only solution would be to
%%%% agglomerate all the different timeSteps polygons into a final one more
%%%% conservative map of the avoidance zones. 
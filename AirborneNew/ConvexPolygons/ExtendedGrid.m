
fakeAvoidanceZoneLat = [34.5 ; 38 ; 37;  35];
fakeAvoidanceZoneLon = [-116; -115.5; -112; -114];
fHi                  = convhull(fakeAvoidanceZoneLat,fakeAvoidanceZoneLon);
fakeHullLatlon       = [fakeAvoidanceZoneLat(fHi) fakeAvoidanceZoneLon(fHi)];

%Extend the network. The ONLY PUNISHMENT WILL BE THE NUMBER OF STEPS TO
%COME BACK TO A predesigned route
LatMax = max(initial_nodes(:,2)); LatMin = min(initial_nodes(:,2));
LonMax = max(initial_nodes(:,3)); LonMin = min(initial_nodes(:,3));
% edges also 
initialEdgesLatlon = NaN(size(initial_edges,1),2*size(initial_edges,2));
for i = 1:size(initial_edges,1)    
    firstPair  = find(initial_nodes(:,1) == initial_edges(i,1));
    initialEdgesLatlon(i,1:2) = initial_nodes(firstPair,2:3);
    secondPair = find(initial_nodes(:,1) == initial_edges(i,2));
    initialEdgesLatlon(i,3:4) = initial_nodes(secondPair,2:3);
end

% to change in a future, we take the Median of lat/lon separation for every
% two pair of points 
LatSpace = median(abs(initial_nodes(2:end,2) - initial_nodes(1:end-1,2)));
LonSpace = median(abs(initial_nodes(2:end,3) - initial_nodes(1:end-1,3)));

genNodesGridLat = (LatMax:-LatSpace:LatMin-LatSpace)'; % LatMin-LatSpace to go a little further
genNodesGridLon = (LonMin:LonSpace:LonMax+LonSpace)';  %''
[NodesLon,NodesLat] = meshgrid(genNodesGridLon,genNodesGridLat);
NodesMatrix = cell(size(NodesLat));

%HIGHLY INEEFFICIENT PROGRAMMING 
NodeID        = initial_nodes(end,1);
dummy         = initial_nodes(end,1);
IndicesMatrix = NaN(numel(NodesMatrix),3);

for kk = 1:size(NodesMatrix,1)
    for hh = 1:size(NodesMatrix,2)
        NodeID = NodeID + 1;
        IndicesMatrix(NodeID-dummy,:) = [NodeID, kk, hh];
        NodesMatrix{kk,hh} = [NodeID, NodesLat(kk,hh), NodesLon(kk,hh)];
    end 
end
NewNodes    = [reshape(NodesLat,[],1),reshape(NodesLon,[],1)];

%% Suppres nodes inside avoidance polygons. 
%%%Both from original and extended grid nodes
%       ,
XV  = fakeAvoidanceZoneLon(fHi); YV = fakeAvoidanceZoneLat(fHi);
X   = initial_nodes(:,3)       ; Y  = initial_nodes(:,2);
toSuppres = inpolygon(X,Y,XV,YV);
initial_nodes(toSuppres,:) = [];
X   = NewNodes(:,2)            ; Y  = NewNodes(:,1);
toSuppres = inpolygon(X,Y,XV,YV);
NewNodes(toSuppres,:) = [];
%%


NewEdges = [];
%Construct edges
for kk = 1:size(NodesMatrix,1)
    for hh = 1:size(NodesMatrix,2)
        originNodeID = NodesMatrix{kk,hh}(1);
% %         if kk~=1
% %             NewEdges = [NewEdges ; [originNodeID NodesMatrix{kk-1,hh}(1)]]; %#ok<AGROW>          
% %         end
        if hh~=1
            NewEdges = [NewEdges ; [originNodeID NodesMatrix{kk,hh-1}(1)]]; %#ok<AGROW>
        end
        if kk~=size(NodesMatrix,1)
            NewEdges = [NewEdges ; [originNodeID NodesMatrix{kk+1,hh}(1)]]; %#ok<AGROW>  
        end 
% %         if hh~=size(NodesMatrix,2)
% %             NewEdges = [NewEdges ; [originNodeID NodesMatrix{kk,hh+1}(1)]]; %#ok<AGROW>
% %         end 
    end 
end 
% Reexpress edges from ID to latlon
NewEdgesLatlon = NaN(size(NewEdges,1),2*size(NewEdges,2)); %latlon pair for one ID
for ii = 1:size(NewEdgesLatlon,1)
    idStart = NewEdges(ii,1);  idFinal = NewEdges(ii,2);
    locStart = find(IndicesMatrix(:,1)==idStart);
    locEnd   = find(IndicesMatrix(:,1)==idFinal);
    rowColStart = IndicesMatrix(locStart,2:end);
    rowColEnd   = IndicesMatrix(locEnd,2:end);
    latlonStart = NodesMatrix{rowColStart(1),rowColStart(2)};
    latlonFinal = NodesMatrix{rowColEnd(1),rowColEnd(2)};
    NewEdgesLatlon(ii,:) = [latlonStart(2:end) latlonFinal(2:end)];
    display(NewEdgesLatlon(ii,:))
end 
%% Suppres Edges inside avoidance polygons. VAYA CHAPUZA
%%%Both from original and extended grid nodes

XV  = fakeAvoidanceZoneLon(fHi); YV = fakeAvoidanceZoneLat(fHi);
X   = NewEdgesLatlon(:,2)       ; Y  = NewEdgesLatlon(:,1);
toSuppres1 = inpolygon(X,Y,XV,YV);
X   = NewEdgesLatlon(:,4)       ; Y  = NewEdgesLatlon(:,3);
toSuppres2 = inpolygon(X,Y,XV,YV);
toSuppres  = toSuppres1 | toSuppres2;
NewEdgesLatlon(toSuppres,:) = [];


X   = initialEdgesLatlon(:,2)       ; Y  = initialEdgesLatlon(:,1);
toSuppres1 = inpolygon(X,Y,XV,YV);
X   = initialEdgesLatlon(:,4)       ; Y  = initialEdgesLatlon(:,3);
toSuppres2 = inpolygon(X,Y,XV,YV);
toSuppres  = toSuppres1 | toSuppres2;
initialEdgesLatlon(toSuppres,:) = [];


%% GRAPHICS
p = gcf;
figure(p.Number+1)
clf
% plot original route nodes 
plot(initial_nodes(:,3),initial_nodes(:,2),'*','Color',[0 0 .8])
hold on 
% plot Avoidance Zones
plot(fakeAvoidanceZoneLon(fHi),fakeAvoidanceZoneLat(fHi),'Color',[0 0 0],'linewidth',3)
plot(fakeAvoidanceZoneLon(fHi),fakeAvoidanceZoneLat(fHi),'k--', 'Color',[1 1 0],'linewidth',2)
hold on
% plot New Generated Nodes
plot(NewNodes(:,2),NewNodes(:,1),'*','Color',[.8 .8 1])
hold on
NewEdgesLatlonDIff  = NewEdgesLatlon(:,3:4)-NewEdgesLatlon(:,1:2);
quiver(NewEdgesLatlon(:,2),NewEdgesLatlon(:,1),NewEdgesLatlonDIff(:,2),NewEdgesLatlonDIff(:,1),...
    0,'Color',[.9 .8 .8])
hold on
InitialEdgesLatlonDIff  = initialEdgesLatlon(:,3:4)-initialEdgesLatlon(:,1:2);
quiver(initialEdgesLatlon(:,2),initialEdgesLatlon(:,1),InitialEdgesLatlonDIff(:,2),InitialEdgesLatlonDIff(:,1),...
    0,'Color',[0.1 0.1 0.1])
hold on
axis([ -118.5300 -111.5600 34.1970   40.8190])
% axis tight
%%% only valid for oure example of SLC 
plot(initial_nodes(1,3),initial_nodes(1,2),'.','markersize',40,'Color',[1 0 0])
text(initial_nodes(1,3),initial_nodes(1,2),'OR')
plot(-118.3800,34.1980,'.','markersize',40,'Color',[1 0 0])
text(-118.3800,34.1980,'DEST')


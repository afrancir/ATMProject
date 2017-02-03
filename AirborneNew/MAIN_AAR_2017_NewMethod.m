clc
clear all %#ok<CLALL>
close all
%%% Computing gridded approximation of planning domain

variable_of_interest = 'Geopotential_height_convective_cloud_top';
LAT_lim              = [16.281 58.3654];
LON_lim              = [-139.8561 -57.3811];
delta_LAT            = 1;
delta_LON            = 1;
files_analyzed       = {'ruc2_130_20110701_0000_000.grb2','ruc2_130_20110701_0100_000.grb2','ruc2_130_20110701_0200_000.grb2',...
                        'ruc2_130_20110701_0300_000.grb2','ruc2_130_20110701_0400_000.grb2','ruc2_130_20110701_0500_000.grb2'};
                    
N_hours          = 5;
t_in             = 1; % [h] gap between the current time and the beginning of the planning horizon
dt               = 5; % [min] time-step of our model
timespan         = 3; % [h] time-span of planning horizon
planning_horizon = 60*t_in:dt:60*(t_in+timespan);
weather_horizon  = 60*t_in:dt:60*N_hours;

[GRID,LONLAT_vec,WEATHER_interp] = planning_domain(LAT_lim,LON_lim,delta_LAT,delta_LON,files_analyzed,variable_of_interest,N_hours,t_in,dt);
%% Generation of the avoidance zones from the Weather product
%####### normalization of the weather weights  ARNAU 7/31/2016
    looper = numel(size(WEATHER_interp));
    max_dummy = WEATHER_interp;
    for i=1:looper
       max_dummy = max(max_dummy,[],'omitnan');
    end 
    max_element  = max_dummy;
    min_element  = 0;
    WEATHER_interp = WEATHER_interp./(max_element-min_element);
%####### the threshold method is quite simple, but enough for the moment
    thresh = 0.9;
    weatherPolygons = WEATHER_interp > thresh;
%% ####### Generating the Convex Clusters enclosing avoidance zones 
cd ConvexPolygons
allLatlonCH = PolygonsProcessorMainAsAFunctionGRID(weatherPolygons, ...
    GRID);
cd ..
%%

%%% Loading a Network component, adjusting the Network component so that is
%%% is compatible with the time-step selected

AIRPORT_ID = 15;  % in this case, SLC  ##Arnau: (We are only studying Salt Lake City)
ROUTE_ID   = 1;   % internal flights   ##Arnau: Only internal flights? Why?
PATH       = pwd; % path to current folder
% storing initial network (i.e., all internal flights within  SLC center)
network    = load(strcat(PATH,'/Cells_timestep/',num2str(AIRPORT_ID),'_',num2str(ROUTE_ID),'.txt')); 

% the output is a .txt file written in the directory
get_aggregate_routes_timestep(AIRPORT_ID,ROUTE_ID,dt,PATH)

%%% Defining all O/D pairs for SLC and using the same routine applied to
%%% each specific O/D pair:
% %         routes_ID     = all_OD_pairs{1,OD_pair};
% %% original routes connecting same O/D. For example we can have an 
% %         alt_routes_ID = all_OD_pairs{2,OD_pair};
% %% that means, all the alternative routes available for our initial routes
% %         ref_airport   = all_OD_pairs{3,OD_pair};
% %% Departure Airport ID latLon (always 2?)
all_OD_pairs       = readAllODPairs();

%%% Focusing on a specific O/D pair. Creating local network, solving
%%% minimum cost path problem and storing new edges/nodes that have been
%%% created

ALL_NEW_CONN  = [];
ALL_NEW_NODES = [];

% Initialization of the output containing all controls for the network
% component considered. Each column defines a specific destination airport,
% while the first row defines ground holding, the second row pre-departure
% rerouting and the third airborne rerouting

% CONTROLS{1,x} = gh_ODPair or ground holding Control.
% CONTROLS{2,x} = gh_ODPair or ground holding Control.
%
CONTROLS = cell(3,size(all_OD_pairs,2)); 

% BlockedRoutesApart
for OD_pairID=1:length(all_OD_pairs(1,:))
    
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp(['Focusing on airport ',num2str(OD_pairID)])
        
    if isnan(all_OD_pairs{1,OD_pairID}) == 0                    
        
        CONTROLS{1,OD_pairID} = [];
        CONTROLS{2,OD_pairID} = [];
        CONTROLS{3,OD_pairID} = [];
        
        
    else
        routes_ID     = all_OD_pairs{1,OD_pairID};
        alt_routes_ID = all_OD_pairs{2,OD_pairID};
        ref_airport   = all_OD_pairs{3,OD_pairID};
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Storing ground holding controls %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Creating a cell array with as many elements as the number of
        % routes going to this particular destination
        gh_OD_pair = cell(numel(routes_ID),1);
        for i=1:numel(routes_ID)
            nom_path            = find(network(:,1)==routes_ID(i));
            alt_path            = [nom_path(1);nom_path]; %?
            gh                  = cell(3,1);
            gh{1,1}             = nom_path(1); %
            gh{2,1}             = nom_path; %whats the point of that??
            gh{3,1}             = alt_path;
            gh_OD_pair{i,1}     = gh;
        end
        
        CONTROLS{1,OD_pairID}     = gh_OD_pair;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Pre-departure and Airborne Rerouting controls %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%% FIRST WE NEED TO IDENTIFY THOSE FLIGHTS WHICH ARE UNDER CONFLICT
        %%% AND JUST IN THAT CASE CREATE THE AUGMENTED NETWORK
                
    end
    
end

%%% Now we have all the new connections created when introducing airborne
%%% rerouting. We need to assign a node ID to each new node and introduce
%%% all the new edges. Additionally, each connection defines a specific
%%% control and we need to map this control for the optimization problem

cont    = length(network(:,1))+1;

NEW_EDGES_ALL       = [];
NEW_NODES_ALL       = [];
NODES_FOR_EACH_CONN = cell(length(ALL_NEW_CONN(:,1)),1);

for i=1:length(ALL_NEW_CONN(:,1))
    
    % The number of new edges we are introducing is the number of new nodes
    % we are introducing plus one
    new_edges        = ALL_NEW_CONN(i,5)+1;
    % Initialization of matrix where we will store all new edges associated
    % with this particular connection
    new_edges_matrix = zeros(new_edges,2);
    
    if new_edges==1 % in this case, no additional node has been created
        new_edges_matrix(1,:) = [ALL_NEW_CONN(i,2) ALL_NEW_CONN(i,4)];
        
        NEW_EDGES_ALL          = [NEW_EDGES_ALL;new_edges_matrix];
        NODES_FOR_EACH_CONN{i} = NaN;
        
    else
        new_nodes              = ALL_NEW_CONN(i,5);       % number of new nodes
        new_nodes_matrix       = zeros(new_nodes,5);
        new_nodes_ID           = cont:cont+new_nodes-1;
        NODES_FOR_EACH_CONN{i} = transpose(new_nodes_ID);
        for j=1:new_edges
            % first edge connects node where we leave original route with
            % first newly created node
            if j==1
                new_edges_matrix(j,:) = [ALL_NEW_CONN(i,2) new_nodes_ID(j)];
                % last edge connects last newly created node with node of the
                % alternate route where we are merging
            elseif j==new_edges
                new_edges_matrix(j,:) = [new_nodes_ID(end) ALL_NEW_CONN(i,4)];
                % otherwise, we have an edge connecting 2 newly created nodes
            else
                new_edges_matrix(j,:) = [new_nodes_ID(j-1) new_nodes_ID(j)];
            end
        end
        
        for k=1:new_nodes
            new_nodes_matrix(k,:) = [new_nodes_ID(k) ALL_NEW_NODES{i}(k,2:5)];
        end
        
        NEW_EDGES_ALL = [NEW_EDGES_ALL;new_edges_matrix];
        NEW_NODES_ALL = [NEW_NODES_ALL;new_nodes_matrix];
        
        % Updating value for new node ID
        cont = new_nodes_ID(end)+1;
        
    end
end

%%% We can now define airborne rerouting controls as well, using the same
%%% policy. Every control will be associated with a cell with 3 elements.
%%% The first one represents the ID of the node where control is carried
%%% out, the second one the path the aircraft would follow with no control,
%%% the third the path the aircraft would follow if control is activated

for OD_pairID=1:length(all_OD_pairs(1,:))
    % Determining if we have any rerouting connection for this specific
    % destination
    idx_controls_for_this_ODpair = find(ALL_NEW_CONN(:,6)==OD_pairID);
    N_controls_for_this_ODpair   = numel(idx_controls_for_this_ODpair);
    % If yes, store all the controls
    if ~isempty(idx_controls_for_this_ODpair)
        arr = cell(N_controls_for_this_ODpair,1);
        for i=1:N_controls_for_this_ODpair
            arr_i     = cell(3,1);                                                               % cell that will contain information for this specific control
            nom_route = network(ALL_NEW_CONN(idx_controls_for_this_ODpair(i),2),1);              % ID of nominal route
            alt_route = network(ALL_NEW_CONN(idx_controls_for_this_ODpair(i),4),1);              % ID of alternate route
            idx_nom   = find(network(:,1)==nom_route);                                           % find all nodes of nominal route
            idx_alt   = find(network(:,1)==alt_route);                                           % find all nodes of alternate route
            nom_path  = transpose(ALL_NEW_CONN(idx_controls_for_this_ODpair(i),2):idx_nom(end)); % path without control
            
            % we have some intermediate nodes
            if ~isnan(NODES_FOR_EACH_CONN{idx_controls_for_this_ODpair(i)})
                alt_path = [ALL_NEW_CONN(idx_controls_for_this_ODpair(i),2);...
                            NODES_FOR_EACH_CONN{idx_controls_for_this_ODpair(i)};...
                            transpose(ALL_NEW_CONN(idx_controls_for_this_ODpair(i),4):idx_alt(end))];
            
            % we go directly from the nominal to the alternate route
            else
                alt_path = [ALL_NEW_CONN(idx_controls_for_this_ODpair(i),2);...
                            transpose(ALL_NEW_CONN(idx_controls_for_this_ODpair(i),4):idx_alt(end))];
            end
            
            arr_i{1,1} = ALL_NEW_CONN(idx_controls_for_this_ODpair(i),2);
            arr_i{2,1} = nom_path;
            arr_i{3,1} = alt_path;
            arr{i,1}   = arr_i;
            
        end
        
        CONTROLS{3,OD_pairID} = arr;
        
    % Otherwise, do nothing
    else
        CONTROLS{3,OD_pairID} = [];
    end
    
    
    
end

%%% Creating a matrix that has the same size of the cell array CONTROLS.
%%% All elements characterized by a 0 are associated with a
%%% control/destination pair that has at least one control. Elements with a
%%% 1 are associated with a control/destination pair that shows no controls 

check_controls = cellfun('isempty',CONTROLS);

% Element map_controls(i,j) defines the number of controls of type i (1 =
% ground holding, 2 = pre-departure rerouting, 3 = airborne rerouting),
% associated with the destination j. At least one column will show three
% zeros (i.e., the column associated with the origin airport defining this
% specific network component)
map_controls   = zeros(length(CONTROLS(:,1)),length(CONTROLS(1,:)));  

for i=1:length(CONTROLS(:,1))
    for j=1:length(CONTROLS(1,:))
        % in this case, the associated cell array element was empty, thus
        % we don't have any control associated
        if check_controls(i,j)
        else
            map_controls(i,j) = length(CONTROLS{i,j});
        end
    end
end

% Getting overall number of controls for the network component
N_controls        = sum(sum(map_controls));
% Defining cell array that will store all the controls
Controls_sequence = cell(N_controls,1);
% Initializing indexing to store controls in the cell array
cont_controls     = 1;

for i=1:length(CONTROLS(:,1))
    for j=1:length(CONTROLS(1,:))
        if map_controls(i,j)~=0
            idx = cont_controls:cont_controls+length(CONTROLS{i,j})-1;
            for k=1:length(CONTROLS{i,j})
                dummy                       = cell(4,1);
                dummy{1}                    = i;
                dummy{2}                    = CONTROLS{i,j}{k}{1};
                dummy{3}                    = CONTROLS{i,j}{k}{2};
                dummy{4}                    = CONTROLS{i,j}{k}{3};
                Controls_sequence{idx(k),1} = dummy;
            end
            cont_controls = cont_controls+length(CONTROLS{i,j});
        else
        end
    end
end


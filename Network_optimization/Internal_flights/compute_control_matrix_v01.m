% Written by Alessandro Bombelli, 27th Sept. 2016
% Function that computes matrices Da and Db for a specific network
% component
% Note: in v01 only ground holding and pre-departure rerouting are
% considered. The extension to airborne rerouting follows the same routine

function [gh_ID,prr_mat] = compute_control_matrix_v01(ID_AIRP)

airports     = load(strcat(pwd,'/NETWORK/AIRPORTS/AIRPORTS_coord.txt'));
orig_airport = airports(airports(:,1)==ID_AIRP,:);

AGG_ROUTES =  load(strcat(pwd,'/NETWORK/INT_AggRoutes_timestep/',num2str(orig_airport(1)),'/',num2str(orig_airport(1)),'.txt'));
N_r        = AGG_ROUTES(end,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 1: define the overall number of controls %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ground holding: as many controls as the number of aggregate routes

Nc_gh = N_r;          % number of ground holding controls
gh_ID = zeros(N_r,1);

for i=1:N_r
    idx          = find(AGG_ROUTES(:,1)==i);
    gh_ID(i)     = idx(1); % first node of the current aggregate route
end

% find all possible destinations
possible_destinations = unique(AGG_ROUTES(:,3));

prr_mat = [];

for i=1:numel(possible_destinations)
    idx_this_dest     = find(AGG_ROUTES(:,3)==possible_destinations(i));
    ar_to_this_dest   = unique(AGG_ROUTES(idx_this_dest,1));
    N_ar_to_this_dest = numel(ar_to_this_dest);
    % Only one aggregate route to this destination, no pre-departure
    % rerouting
    if N_ar_to_this_dest == 1
    % At least two routes. We have pre-departure rerouting options
    else
        idx_origin_node = zeros(N_ar_to_this_dest,1);
        % vector where every element is the ID of the first node of the
        % associated route
        for j=1:N_ar_to_this_dest
            idx_this_route     = find(AGG_ROUTES(:,1)==ar_to_this_dest(j));
            idx_origin_node(j) = idx_this_route(1);
        end
        for j=1:N_ar_to_this_dest
            original_node = idx_origin_node(j)*ones(N_ar_to_this_dest-1,1);
            new_node      = setdiff(idx_origin_node(j),idx_origin_node);
            prr_mat       = vertcat(prr_mat,horzcat(original_node,new_node));
        end
    end
end






return
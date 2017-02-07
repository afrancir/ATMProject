% Written by Alessandro Bombelli, 6 February 2017
%_________________________________________________________________________%
% Simple example to show how to build all the matrices of interest (both
% for the dynamics and for the constraints) for a single network component.
% The idea can be easily extended to more network components using
% horizontal concatenation for the constraints that couple different
% network components, and a block diagonal structure for teh dynamics and
% the other constraints
%
% Reference: Strategic Air Traffic Planning with Fréchet Distance 
% Aggregation and Rerouting, Alessandro Bombelli, Lluis Soler, 
% Eric Trumbauer, Kenneth D. Mease, Journal of Guidance, Control, and
% Dynamics

clc
clear all
close all

%%## Arnau 2/3/2017 - adaptation for diff. OS //Start
OSTypeString = {'win' ,'lin'};
winOrLin = input('Please press "1" for a Windows OS or "2" for a Linux-based OS \n');
if (winOrLin==1 || winOrLin==2)
    OS = OSTypeString(winOrLin);
else 
    msg ='This is not a valid option. 1 for Win or 2 for Linux';
    error(msg)
end 


if strcmp(winOrLin,'win')
    AIRPORTS    = load(strcat(pwd,'\AIRPORTS\AIRPORTS_coord.txt'));
else
    AIRPORTS    = load(strcat(pwd,'/AIRPORTS/AIRPORTS_coord.txt'));
end 
%%## Arnau 2/3/2017 - adaptation for diff. OS //End

airport_ID = AIRPORTS(:,1);

% Focusing on internal flow from LAX
this_airport = 34;

% Loading Network
if strcmp(winOrLin,'win')
    Network = load(strcat(pwd,'\Internal_flights\NETWORK\INT_AggRoutes_timestep\'...
        ,num2str(this_airport),'\',num2str(this_airport),'.txt'));
else
    Network = load(strcat(pwd,'/Internal_flights/NETWORK/INT_AggRoutes_timestep/'...
        ,num2str(this_airport),'/',num2str(this_airport),'.txt'));
end 

%_________________________________________________________________________%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Definition of matrices A and C %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_routes = Network(end,1);
A_temp   = [];
C_temp   = [];
for i=1:N_routes
    nodes_this_route = numel(find(Network(:,1)==i));
    A_temp           = blkdiag(A_temp,diag(ones(nodes_this_route-1,1),-1));
    C_temp           = blkdiag(C_temp,vertcat(1,zeros(nodes_this_route-1,1)));
end

A = sparse(A_temp);
C = sparse(C_temp);

clear A_temp C_temp

sparsity_A = nnz(A)/length(Network(:,1))^2*100;
sparsity_C = nnz(C)/(length(Network(:,1))*N_routes)*100;

%_________________________________________________________________________%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Determination of controls:   %%%
%%% GH = ground holding          %%%
%%% PR = pre-departure rerouting %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

other_airports = setdiff(airport_ID,this_airport);

GH = [];
PR = [];

for j=1:numel(other_airports)
    idx_all_entries_to_this_destination = find(Network(:,3)==other_airports(j));
    if ~isempty(idx_all_entries_to_this_destination)
        ID_routes = unique(Network(idx_all_entries_to_this_destination,1));
        disp(['Number of routes between airport ',num2str(this_airport),' and airport ',num2str(other_airports(j)),': ',num2str(numel(ID_routes))]);
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        %%%%%%%%%%%%%%%%%%%%%%
        %%% Ground holding %%%
        %%%%%%%%%%%%%%%%%%%%%%
        for i=1:numel(ID_routes)
            all_idx         = find(Network(:,1)==ID_routes(i));
            idx_origin_node = all_idx(1);
            GH              = vertcat(GH,[1 idx_origin_node idx_origin_node]);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Pre-departure rerouting %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if numel(ID_routes)>1
            origin_node_ID = zeros(numel(ID_routes),1);
            for ii=1:numel(ID_routes)
                all_idx_this_route = find(Network(:,1)==ID_routes(ii));
                origin_node_ID(ii) = all_idx_this_route(1);
            end
            for ii=1:numel(origin_node_ID)
                other_routes = setdiff(origin_node_ID,origin_node_ID(ii));
                this_PR      = zeros(numel(other_routes),3);
                for jj=1:numel(other_routes)
                    this_PR(jj,:) = [2 origin_node_ID(ii) other_routes(jj)];
                end
                PR = vertcat(PR,this_PR);
            end
        end
    else
        disp(['Number of routes between airport ',num2str(this_airport),' and airport ',num2str(other_airports(j)),': ',num2str(0)]);
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    end
end

%_________________________________________________________________________%

% Note: this code works if both ground holding and pre-departure rerouting
% are considered. If only ground holding is considered, the code will
% output an error. Please modify accordingly in that case

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Determination of control matrices Da and Db %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Da_temp_GH = zeros(length(Network(:,1)),length(GH(:,1)));
Da_temp_PR = zeros(length(Network(:,1)),length(PR(:,1)));
Db_temp_GH = zeros(length(Network(:,1)),length(GH(:,1)));
Db_temp_PR = zeros(length(Network(:,1)),length(PR(:,1)));

for i=1:length(GH(:,1))
    Da_temp_GH(GH(i,2),i) = -1;
    Db_temp_GH(GH(i,2),i) = 1;
end

for i=1:length(PR(:,1))
    Da_temp_PR(PR(i,2),i) = -1;
    Da_temp_PR(PR(i,3),i) = 1;
end

Da = sparse(horzcat(Da_temp_GH,Da_temp_PR));
Db = sparse(horzcat(Db_temp_GH,Db_temp_PR));

clear Da_temp_GH Da_temp_PR Db_temp_GH Db_temp_PR;

sparsity_Da = nnz(Da)/(length(Network(:,1))*(length(GH(:,1))+length(PR(:,1))))*100;
sparsity_Db = nnz(Db)/(length(Network(:,1))*(length(GH(:,1))+length(PR(:,1))))*100;

%_________________________________________________________________________%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Determination of full matrices A_nc, C_nc, B_nc %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The network is associated to a 5 minute time-step. With 48 time-steps, we
% are considering a planning horizon of 4 hours
N_t = 48;

A_nc_temp  = [];
C_nc_temp  = [];
B1_nc_temp = [];
B2_nc_temp = [];

%%%%%%%%%%%%%%%%%%%%%%%
%%% Assembling A_nc %%%
%%%%%%%%%%%%%%%%%%%%%%%
for i=1:N_t
    
    if i == 1
        new_block = A;
    else
        new_block = A*new_block;
    end
    
    A_nc_temp = blkdiag(A_nc_temp,new_block);
end

A_nc = sparse(A_nc_temp);

clear A_nc_temp;

%%%%%%%%%%%%%%%%%%%%%%%
%%% Assembling C_nc %%%
%%%%%%%%%%%%%%%%%%%%%%%
for i=1:N_t
    row = [];
   for j=1:N_t
       if j<=i
           if j==1
               new_block = C;
           else
               new_block = A*new_block;
           end
           row = horzcat(row,new_block);
       else
           row = horzcat(row,zeros(size(C)));
       end
   end
   C_nc_temp = vertcat(C_nc_temp,row);
end

C_nc = sparse(C_nc_temp);

clear C_nc_temp;


%%%%%%%%%%%%%%%%%%%%%%%
%%% Assembling B_nc %%%
%%%%%%%%%%%%%%%%%%%%%%%

B_nc_temp = [];

for i=1:N_t
    row = [];
    for j=1:N_t
        if j<i
            block = A^(i-j)*Da+A^(i-j-1)*Db;
        elseif j==i
            block = Da;
        else
            block = zeros(size(Da));
        end
        row = horzcat(row,block);
        row = sparse(row);
    end
    B_nc_temp = vertcat(B_nc_temp,row);
end

B_nc = sparse(B_nc_temp);

clear B_nc_temp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initial condition vector %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X0         = zeros(length(Network(:,1)),1);
X0_nc_temp = repmat(X0,N_t,1);
X0_nc      = sparse(X0_nc_temp);

clear X0_nc_temp;

% Note: the full dynamics are described with the following matrix equation
%
% X = A_nc*X0_nc+C_nc*b+B_nc*u
%
% where:
%
% X: state variables for all different time-steps, vertically stacked 
%    X = [X(t1) X(t2) ... X(t_N_t)]'
%
% A_nc*X0: term mapping the initial condition into the dynamics
%
% C_nc*b: term mapping scheduled departures into the dynamics
%
% B_nc*u: term mapping controls into the dynamics


%_________________________________________________________________________%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Inequality Constraints %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Scheduled departures: we schedule one aircraft to depart from each route
% at each time-step
b = 1*ones(length(GH(:,1))*N_t,1); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Model Consistency Constraints %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ALL_C   = vertcat(GH,PR);
N_r     = length(GH(:,1));
N_c     = length(ALL_C(:,1));
MC_temp = zeros(N_r*N_t,N_c*N_t);

for i=1:N_t
    for j=1:N_r
        idx_origin   = GH(j,2);
        idx_controls = find(ALL_C(:,2)==idx_origin);
        idx_GH       = find(ALL_C(:,2)==idx_origin & ALL_C(:,3)==idx_origin);
        
        MC_temp(N_r*(i-1)+j,N_c*(i-1)+idx_controls) = 1;
        % For all time-steps greater than t1, we need to account the effect
        % of flights held on the ground from the previous time-step
        if i>1
            MC_temp(N_r*(i-1)+j,N_c*(i-2)+idx_GH) = -1;
        else
        end
    end
end

MC = sparse(MC_temp);

clear MC_temp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Sector Capacity Constraint %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note: in this example we set a very high threshold for Sector Capacity.
% This constraints is not active as a consequence (i.e., it is never
% violated)

% Obtain IDs of all sectors
all_sectors_IDs = unique(Network(:,end));
N_s             = numel(all_sectors_IDs);

% Node-to-sector matrix
N2S_mat_temp = zeros(N_s,length(Network(:,1)));
for i=1:N_s
    % Determine which sector we are considering
    this_sector     = all_sectors_IDs(i);
    % Determine all nodes belonging to the current sector
    idx_this_sector = find(Network(:,end)==this_sector);
    % Modify elements of N2S_mat accordingly
    N2S_mat_temp(i,idx_this_sector) = 1;
end

N2S_mat = sparse(N2S_mat_temp);

clear N2S_mat_temp;

% Now, build a block-diagonal structure to replicate the N2S_mat matrix to
% account for all the different time-steps we are accounting for

MS_temp = [];

for i=1:N_t
    MS_temp = blkdiag(MS_temp,N2S_mat);
end

MS = sparse(MS_temp);

clear MS_temp;

b_S = 200*ones(N_s*N_t,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Departure Capacity constraints %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note: in this example we set a very high threshold for Departure 
% Capacity. This constraints is not active as a consequence (i.e., 
% it is never violated)

idx_departures                 = GH(:,2);
dep_vec_temp                   = zeros(1,length(Network(:,1)));
dep_vec_temp(1,idx_departures) = 1;
MD_temp                   = [];

n_h                            = 4; % [h] planning horizon is 4 hours

for i=1:4
    % we repeat the vector 12 times since we have 12 time-steps in 1 hour
    % with a 5-minute time-step. With a different time-step, the following
    % code line should be changed accordingly
    MD_temp = blkdiag(MD_temp,repmat(dep_vec_temp,1,12));
end

MD = sparse(MD_temp);

clear MD_temp;

b_D = 500*ones(n_h,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% No-fly-zones constraints %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We block three routes for a determined time interval

idx_nfz        = [70;160;280];
t_nfz          = zeros(numel(idx_nfz),N_t);
t_nfz(1,13:24) = 1;
t_nfz(2,13:24) = 1;
t_nfz(3,13:24) = 1;
nfz_mat_temp   = [];

for i=1:numel(idx_nfz)
    this_blocked_node           = idx_nfz(i);
    timesteps_this_blocked_node = find(t_nfz(i,:));
    for j=1:numel(timesteps_this_blocked_node)
        nfz_mat_temp_new_row = zeros(1,N_t*length(Network(:,1)));
        nfz_mat_temp_new_row(1,(timesteps_this_blocked_node(j)-1)*length(Network(:,1))+this_blocked_node) = 1;
        nfz_mat_temp = vertcat(nfz_mat_temp,nfz_mat_temp_new_row);
    end
end

M_NFZ = sparse(nfz_mat_temp);

clear nfz_mat_temp;

b_NFZ = zeros(length(M_NFZ(:,1)),1);

%_________________________________________________________________________%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Definition of inequality matrices and vectors %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For model consistency constraints, we directly work on the constraint
% vector. The other constraints are manipulations of the state variables
% (i.e., of the dynamics), and thus are obtained pre-multiplying the
% dynamics with the matrices we defined above

A_MC     = MC;
b_MC     = b;

A_SecCap = MS*B_nc;
b_SecCap = b_S-MS*(A_nc*X0_nc+C_nc*b);

A_DepCap = MD*B_nc;
b_DepCap = b_D-MD*(A_nc*X0_nc+C_nc*b);

A_NoFlyZones = M_NFZ*B_nc;
b_NoFlyZones = b_NFZ-M_NFZ*(A_nc*X0_nc+C_nc*b);

%_________________________________________________________________________%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Lower boundary for controls %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set a lower boundary for controls, i.e., zero. Note tha an upper boundary
% is not necessary, since it is implicitely applied via model consistency
% constraints

LB = zeros(length(A_MC(1,:)),1);

%_________________________________________________________________________%

%%%%%%%%%%%%%%%%%%%%%
%%% Cost function %%%
%%%%%%%%%%%%%%%%%%%%%

w_length_GH = ones(length(GH(:,1)),1);
w_length_PR = zeros(length(PR(:,1)),1);

% Scan all the different pre-departure rerouting constraints
for i=1:length(PR(:,1))
    idx_original_route    = find(GH(:,2)==PR(i,2));
    idx_reroute           = find(GH(:,2)==PR(i,3));
    length_original_route = numel(find(Network(:,1)==idx_original_route));
    length_reroute        = numel(find(Network(:,1)==idx_reroute));
    w_length_PR(i)          = max(1,length_reroute-length_original_route);
    
end

w_time = vertcat(w_length_GH,w_length_PR);

w = [];

for i=1:N_t
    w_this_time = w_time*(N_t+1-i)/N_t;
    w           = vertcat(w,w_this_time);
end

%_________________________________________________________________________%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Solving ILP problem %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

A_ineq = vertcat(A_MC,A_SecCap,A_DepCap,A_NoFlyZones);
b_ineq = vertcat(b_MC,b_SecCap,b_DepCap,b_NoFlyZones);

%intcon = ones(numel(w),1);
%X      = intlinprog(w,intcon,A_ineq,b_ineq,[],[],LB,[]);


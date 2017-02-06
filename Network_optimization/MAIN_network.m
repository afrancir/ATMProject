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

%%## Arnau 2/4/2017 - comment
%For this example only focusing on departures from this_airport:
this_airport = 34;

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

% this is only a chech of how sparse they are, not used in practice
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

% this is only a chech of how sparse they are, not used in practice
sparsity_Da = nnz(Da)/(length(Network(:,1))*(length(GH(:,1))+length(PR(:,1))))*100;
sparsity_Db = nnz(Db)/(length(Network(:,1))*(length(GH(:,1))+length(PR(:,1))))*100;

%_________________________________________________________________________%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Determination of full matrices A_nc, C_nc, B_nc %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        new_block = new_block*A;
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
for i=1:N_t
    row_B1 = [];
   for j=1:N_t
       if j<=i
           if j==1
               new_block = Da;
           else
               new_block = A*new_block;
           end
           row_B1 = horzcat(row_B1,new_block);
       else
           row_B1 = horzcat(row_B1,zeros(size(Da)));
       end
   end
   B1_nc_temp = vertcat(B1_nc_temp,row_B1);
end

for i=1:N_t
    if i==1
        row_B2 = zeros(length(Network(:,1)),length(Da(1,:))*N_t);
    else
        row_B2 = [];
        for j=1:N_t
            if j<i
                new_block = A^(i-j-1)*Db;
            else
                new_block = zeros(size(Db));
            end
            
            row_B2 = horzcat(row_B2,new_block);
        end
    end
    B2_nc_temp = vertcat(B2_nc_temp,row_B2);
end

B_nc = sparse(B1_nc_temp+B2_nc_temp);

clear B1_nc_temp B2_nc_temp;

%%## Arnau 2/5/2017 - Whole section INTEGER LINEAR PROGRAM PART //Start
%% INTEGER LINEAR PROGRAM PART
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N_s : number of sectors 
% cMax: (N_t*N_s)x1 vector defining sector capacity at each time. 
% dMax: (N_t*N_a)x1 vector defining departure capacity at each time. 
% aMax: (N_t*N_a)x1 vector defining arrival capacity at each time. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% temporarily setting some variables to easy values for testing:
N_s  = 1;
cMax = inf(N_s*N_t,1); %no capacity limit 

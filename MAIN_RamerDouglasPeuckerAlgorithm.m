% Written by Alessandro Bombelli, 17th January 2017
% Ramer-Douglas-Peucker algorithm to undersample a trajectory defined by a
% sequence of great circle arcs. In the original framework, the algorithm
% is applied to the planar case (thus, the Euclidean distance is used). In
% this formulation, the Haversine distance is used instead, to account for
% the fact that great circle arcs are positioned onto a sphere
%
% In the example provided here, a recorded flight trajectory from Las Vegas
% McCarran International Airport to Denver International Airport is loaded.
% Then, the function RamerDouglasPeucker is used to obtain a reduced
% trajectory that is a close approximation of the original trajectory. The
% quality of the approximation is controlled by the parameter epsilon
%
% References:
% Urs Ramer, "An iterative procedure for the polygonal approximation of 
% plane curves", Computer Graphics and Image Processing, 1(3), 
%244–256 (1972)
% David Douglas & Thomas Peucker, "Algorithms for the reduction of the 
% number of points required to represent a digitized line or its 
% caricature", The Canadian Cartographer 10(2), 112–122 (1973) 

clc
clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loading the Centers that form the planning domain %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid     = fopen(strcat(pwd,'/AUX_DATA/PLANNING_DOMAIN.txt'), 'rt');
CENTERS = textscan(fid,'%s');
fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loading the airports that form the planning domain %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
airports   = load(strcat(pwd,'/NETWORK/AIRPORTS/AIRPORTS_coord.txt'));

origin_airport      = 33; % Las Vegas McCarran International Airport
destination_airport = 18; % Denver International Airport

% Loading all trajectories from Las Vegas McCarran International Airport to
% Denver International Airport
ALL_TRAJ       = load(strcat(pwd,'/NETWORK/INT_traj_filt_info/',num2str(origin_airport),'/',num2str(origin_airport),'_',num2str(destination_airport),'_int.txt'));
original_IDs   = unique(ALL_TRAJ(:,1));

% Selecting the first trajectory of the set
Trajectory = ALL_TRAJ(ALL_TRAJ(:,1)==original_IDs(1),4:5);

% Definition of the distance threshold (in radians)
epsilon = 0.005/6378;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The Ramer-Douglas-Peucker algorithm is used to reduce the %%%
%%% initial size of the trajectory                            %%%
%%% idx_to_keep = indices of the original trajectory that     %%%
%%% define the reduced trajectory                             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idx_to_keep = RamerDouglasPeucker(Trajectory,epsilon);

% Obtain the reduced trajectory by only selecting the indices that were
% computed by the Ramer-Douglas-Peucker algorithm 
Trajectory_reduced = Trajectory(find(idx_to_keep),:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot original trajectory vs. reduced trajectory %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
hold on
for ii=1:length(Trajectory(:,1))-1
    [lat,lon] = gcwaypts(Trajectory(ii,1),Trajectory(ii,2),Trajectory(ii+1,1),Trajectory(ii+1,2));
    plot(lon,lat,'Color','k','Linestyle','-','Linewidth',2,'Marker','*','Markersize',4)
end
for ii=1:length(Trajectory_reduced(:,1))-1
    [lat,lon] = gcwaypts(Trajectory_reduced(ii,1),Trajectory_reduced(ii,2),Trajectory_reduced(ii+1,1),Trajectory_reduced(ii+1,2));
    plot(lon,lat,'Color','b','Linestyle','-','Linewidth',2,'Marker','none')
end
plot(airports(find(airports(:,1)==origin_airport),3),airports(find(airports(:,1)==origin_airport),2),'Color','r','Linewidth',2,'Marker','s','Markersize',10)
plot(airports(find(airports(:,1)==destination_airport),3),airports(find(airports(:,1)==destination_airport),2),'Color','r','Linewidth',2,'Marker','s','Markersize',10) 
grid on
axis equal
xlabel('Longitude [deg]','Fontname','Avantgarde','Fontsize',14)
ylabel('Latiitude [deg]','Fontname','Avantgarde','Fontsize',14)
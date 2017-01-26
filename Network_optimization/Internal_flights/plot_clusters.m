function plot_clusters(ID1,ID2)


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


origin_airport = airports(airports(:,1)==ID1,:);
dest_airport   = airports(airports(:,1)==ID2,:);

ALL_TRAJ = load(strcat(pwd,'/NETWORK/INT_traj_filt_info/',num2str(origin_airport(1)),...
                '/',num2str(origin_airport(1)),'_',num2str(dest_airport(1)),'_int.txt'));
CLUSTERS = load(strcat(pwd,'/NETWORK/INT_Clusters/',num2str(origin_airport(1)),'/',num2str(origin_airport(1)),'_',num2str(dest_airport(1)),'.txt'));


figure()
hold on
for i=1:numel(CLUSTERS(:,1))
    this_cluster = CLUSTERS(i,1:numel(find(CLUSTERS(i,:)~=0)));
    RGB          = [rand(1) rand(1) rand(1)];
    for j=1:numel(this_cluster)
        this_ID = this_cluster(j);
        idx     = find(ALL_TRAJ(:,1)==this_ID);
        plot(ALL_TRAJ(idx,5),ALL_TRAJ(idx,4),'Color',RGB,'Linewidth',1,'Linestyle','-','Marker','none','Markersize',4)
    end
end
plot(airports(:,3),airports(:,2),'Color','r','Linewidth',2,'Linestyle','none','Marker','s','Markersize',10)
plot(origin_airport(3),origin_airport(2),'Color','b','Linewidth',2,'Linestyle','none','Marker','s','Markersize',10)
plot(dest_airport(3),dest_airport(2),'Color','k','Linewidth',2,'Linestyle','none','Marker','s','Markersize',10)
for i=1:length(CENTERS{1,1})
    current_Center_latlon = load(strcat(pwd,'/NETWORK/Center_boundaries/',char(CENTERS{1,1}(i)),'.txt'));
    plot(current_Center_latlon(:,2),current_Center_latlon(:,1),'Color','k','Linewidth',2,'Marker','none')
end
grid on
axis equal
xlabel('Longitude [deg]','Fontname','Avantgarde','Fontsize',14)
ylabel('Latiitude [deg]','Fontname','Avantgarde','Fontsize',14)
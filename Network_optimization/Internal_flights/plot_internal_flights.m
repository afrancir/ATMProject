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
airport_ID = airports(:,1);

%%%%%%%%%%%%%%
%%% Colors %%%
%%%%%%%%%%%%%%

RGB1 = 0.6*ones(1,3);
RGB2 = 0.1*ones(1,3);

figure()
hold on
%for i=1:numel(airport_ID)
for i=3    
    for j=1:numel(airport_ID)
        if exist(strcat(pwd,'/NETWORK/INT_traj_filt_info/',num2str(airport_ID(i)),...
                '/',num2str(airport_ID(i)),'_',num2str(airport_ID(j)),'_int.txt'), 'file')
            this_OD = load(strcat(pwd,'/NETWORK/INT_traj_filt_info/',num2str(airport_ID(i)),...
                '/',num2str(airport_ID(i)),'_',num2str(airport_ID(j)),'_int.txt'));
            plot(this_OD(:,5),this_OD(:,4),'Color',RGB1,'Linestyle','none','Linewidth',.5,'Marker','*','Markersize',.5)
        else
        end
    end
end
plot(airports(:,3),airports(:,2),'Color',RGB1,'Linewidth',1.5,'Linestyle','none','Marker','s','Markersize',8)
for i=1:length(CENTERS{1,1})
    current_Center_latlon = load(strcat(pwd,'/NETWORK/Center_boundaries/',char(CENTERS{1,1}(i)),'.txt'));
    plot(current_Center_latlon(:,2),current_Center_latlon(:,1),'Color','k','Linewidth',2,'Marker','none')
end
grid on
axis equal
xlabel('Longitude [deg]','Fontname','Avantgarde','Fontsize',14)
ylabel('Latiitude [deg]','Fontname','Avantgarde','Fontsize',14)

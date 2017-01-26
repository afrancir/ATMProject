clc
clear all
close all
tic

% Checking if the folder exists already or not
if isequal(exist(strcat(pwd,'/NETWORK/INT_FDMatrix'),'dir'),7) % 7 = directory
else
    mkdir(strcat(pwd,'/NETWORK/INT_FDMatrix'));
end

AIRPORTS    = load(strcat(pwd,'/NETWORK/AIRPORTS/AIRPORTS_coord.txt'));
airport_ID = AIRPORTS(:,1);

for i=1:numel(airport_ID)
    % Checking if the folder exists already or not
    if isequal(exist(strcat(pwd,'/NETWORK/INT_FDMatrix/',num2str(airport_ID(i))),'dir'),7) % 7 = directory
    else
        mkdir(strcat(pwd,'/NETWORK/INT_FDMatrix/',num2str(airport_ID(i))));
    end
    other_airports = setdiff(airport_ID,airport_ID(i));
    for j=1:numel(other_airports)
        disp(['Trajectories from airport ',num2str(airport_ID(i)),' to airport ',num2str(other_airports(j))]);
        if exist(strcat(pwd,'/NETWORK/INT_traj_filt_info/',num2str(airport_ID(i)),...
                '/',num2str(airport_ID(i)),'_',num2str(other_airports(j)),'_int.txt'), 'file')
            this_OD = load(strcat(pwd,'/NETWORK/INT_traj_filt_info/',num2str(airport_ID(i)),...
                '/',num2str(airport_ID(i)),'_',num2str(other_airports(j)),'_int.txt'));
            traj_ID  = unique(this_OD(:,1));
            n_traj   = numel(traj_ID);
            if n_traj>=2
                FM       = zeros(n_traj,n_traj);
                for k=1:n_traj-1
                    P = this_OD(this_OD(:,1)==traj_ID(k),4:5);
                    for kk=k+1:n_traj
                        Q        = this_OD(this_OD(:,1)==traj_ID(kk),4:5);
                        FM(k,kk) = func_discrete_Frechet_lat_lon_NM_v02(P,Q);
                    end
                end
                dlmwrite(strcat(pwd,'/NETWORK/INT_FDMatrix/',num2str(airport_ID(i)),'/',num2str(airport_ID(i)),'_',num2str(other_airports(j)),'.txt'),FM,' ');
            else
            end
        else
        end
        disp('%%%%%%%%%%%%%%%%%%%%%');
    end
end
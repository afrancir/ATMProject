clc
clear all
close all


AIRPORTS    = load(strcat(pwd,'/NETWORK/AIRPORTS/AIRPORTS_coord.txt'));
airport_ID = AIRPORTS(:,1);

for i=1:numel(airport_ID)
    
    Network_component = [];
    cont              = 1;
    
    
    % Checking if the folder exists already or not
    if isequal(exist(strcat(pwd,'/NETWORK/INT_AggRoutes/',num2str(airport_ID(i))),'dir'),7) % 7 = directory
    else
        mkdir(strcat(pwd,'/NETWORK/INT_AggRoutes/',num2str(airport_ID(i))));
    end
    other_airports = setdiff(airport_ID,airport_ID(i));
    for j=1:numel(other_airports)
        
        % Check that this OD pair is characterized by some aggregate routes
        if exist(strcat(pwd,'/NETWORK/INT_AggRoutes/',num2str(airport_ID(i)),...
                '/',num2str(airport_ID(i)),'_',num2str(other_airports(j)),'.txt'), 'file')
            
            % Loading all aggregate routes
            ALL_TRAJ = load(strcat(pwd,'/NETWORK/INT_AggRoutes/',num2str(airport_ID(i)),...
                '/',num2str(airport_ID(i)),'_',num2str(other_airports(j)),'.txt'));
            n_ar_this_OD = ALL_TRAJ(end,1); % number of aggregate routes
            
            for k=1:n_ar_this_OD
                idx               = find(ALL_TRAJ(:,1)==k);
                this_agg_route    = horzcat(cont*ones(numel(idx),1),ALL_TRAJ(idx,2:end));
                Network_component = vertcat(Network_component,this_agg_route);
                cont              = cont+1;
            end
            
            
            
        else
        end
        
    end
    
    dlmwrite(strcat(pwd,'/NETWORK/INT_AggRoutes/',num2str(airport_ID(i)),'/',num2str(airport_ID(i)),'.txt'),Network_component,' ');
    
end

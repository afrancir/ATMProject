function heading = traj_heading(trajectory)

trajectory = trajectory*pi/180; % from [deg] to [rad]
heading    = zeros(numel(trajectory(:,1))-1,1);

for i=1:numel(trajectory(:,1))-1
    lat_a = trajectory(i,1);    % [rad] lat of first point
    lon_a = trajectory(i,2);    % [rad] lon of first point
    lat_b = trajectory(i+1,1);  % [rad] lat of second point
    lon_b = trajectory(i+1,2);  % [rad] lon of second point
    
    X          = cos(lat_b)*sin(lon_b-lon_a);
    Y          = cos(lat_a)*sin(lat_b)-sin(lat_a)*cos(lat_b)*cos(lon_b-lon_a);
    heading(i) = atan2(X,Y)*180/pi;
    
    % Normalizing heading to be in the [0,360 deg] interval
    if heading(i)>=360
        heading(i) = heading(i)-360*floor(heading(i)/360);
    else
    end
    % Normalizing heading to be in the [0,360 deg] interval
    if heading(i)<0
        heading(i)=heading(i)+360;
    else
    end
    
end

return
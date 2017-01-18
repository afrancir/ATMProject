function delta_sigma = Haversine(POS1,POS2,formula)

lat1 = POS1(1)*pi/180; % [rad]
lon1 = POS1(2)*pi/180; % [rad]
lat2 = POS2(1)*pi/180; % [rad]
lon2 = POS2(2)*pi/180; % [rad]


delta_lon = lon2-lon1; % [rad]
delta_lat = lat2-lat1; % [rad]

if formula == 1
    
    % Haversine
    delta_sigma  = 2*asin(sqrt(sin(delta_lat/2)*sin(delta_lat/2)+cos(lat1)*cos(lat2)*sin(delta_lon/2)*sin(delta_lon/2)));
    
else
    
    % Vincenty
    num          = sqrt((cos(lat2)*sin(delta_lon))^2+(cos(lat1)*sin(lat2)-sin(lat1)*cos(lat2)*cos(delta_lon))^2);
    den          = sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(delta_lon);
    delta_sigma  = atan2(num,den);
    
end

return
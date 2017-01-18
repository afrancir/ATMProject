function cart = latlon2cart(POS,r)

lat = POS(1)*pi/180; % [rad]
lon = POS(2)*pi/180; % [rad]

cart = [r*cos(lat)*cos(lon); r*cos(lat)*sin(lon); r*sin(lat)];
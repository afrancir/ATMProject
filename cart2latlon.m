% Input:  cartesian coordinates [x;y;z] of a point on a sphere
% Output: lat-lon pair in degrees

function latlon = cart2latlon(cart)

R = norm(cart);
x = cart(1);
y = cart(2);
z = cart(3);

lat = asin(z/R);
lon = atan2(y,x);

latlon(1) = lat*180/pi;
latlon(2) = lon*180/pi;

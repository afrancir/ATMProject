function [T] = gc_int(Pc,P1,P2,Re,alpha)

formula = 1;                        % formula used to compute great circle
d       = Haversine(P1,P2,formula); % Great circle segment from P1 to P2
Pc_cart = latlon2cart(Pc,Re);       % Converting Pc into cartesian coordinates
P1_cart = latlon2cart(P1,Re);       % Converting P1 into cartesian coordinates
P2_cart = latlon2cart(P2,Re);       % Converting P2 into cartesian coordiantes
xc      = Pc_cart(1);
yc      = Pc_cart(2);
zc      = Pc_cart(3);
x1      = P1_cart(1);
y1      = P1_cart(2);
z1      = P1_cart(3);
x2      = P2_cart(1);
y2      = P2_cart(2);
z2      = P2_cart(3);

% The great circle segment from P1 to P2 can be expressed in parametrized
% form as (1) [x;y;z] = sin((1-f)*d)/sin(d)*[x1;y1;z1]+sin(f*d)/sin(d)*[x2;y2;z2]
% where f ranges from zero to one. At the same time, the locus of points
% that belong to the sphere's surface and that are distant alpha radians
% from Pc along the surface (i.e., the locus of all great circle segments
% that pass through Pc and display a length of 2*alpha) can be expressed as
% (2) 2*(x*xc+y*yc+z*zc)=2*Re^2*cos(alpha). If we combine the two equations,
% we get the intersection between (1) and (2). Whether there are two
% intersections (secant case), there's one intersection (tangent case) or
% no intersection, depends if the following trigonometric equation is
% solvable
% A*cos(f*d)+B*sin(f*d) = C = R*cos(fd-gamma)
A     = x1*xc+y1*yc+z1*zc;
B     = (x2*xc+y2*yc+z2*zc-cos(d)*(x1*xc+y1*yc+z1*zc))/sin(d);
C     = Re^2*(1-2*(sin(alpha/2))^2);
R     = sqrt(A^2+B^2);
gamma = atan2(B,A);

if C == R % tangency
    
    phi = gamma - 2*pi*floor(gamma/(2*pi));
    f   = phi/d-2*pi/d*floor((phi/d)/(2*pi/d));
    LB  = f;
    UB  = f;
    
elseif C > R % no intersection
    
    LB  = -1;
    UB  = -1;
    
    
else % great circle is secant
    
    phi1 = (acos(C/R)+gamma) - 2*pi*floor((acos(C/R)+gamma)/(2*pi));
    phi2 = (2*pi-acos(C/R)+gamma) - 2*pi*floor((2*pi-acos(C/R)+gamma)/(2*pi));
    f1   = phi1/d-2*pi/d*floor((phi1/d)/(2*pi/d));
    f2   = phi2/d-2*pi/d*floor((phi2/d)/(2*pi/d));
    
    if f1 <= 1 && f2 <= 1 % 2 feasible intersections
        
        LB = min(f1,f2);
        UB = max(f1,f2);
        
    elseif (f1 <= 1 && f2 >= 1) || (f1 >= 1 && f2 <= 1) % one feasible intersection
        
        d_P1 = Haversine(P1,Pc,formula);
        
        if d_P1 < alpha
            LB = 0;
            UB = min(f1,f2);
        else
            LB = min(f1,f2);
            UB = 1;
        end
        
    else % no feasible intersections. Need to verify whether great circle segment in internal or external
        
        d_P1 = Haversine(P1,Pc,formula);
        d_P2 = Haversine(P2,Pc,formula);
        
        if d_P1 < alpha && d_P2 < alpha
            
            LB = 0;
            UB = 1;
            
        else
            
            LB = -1;
            UB = -1;
            
        end
        
    end
    
end

T = [LB;UB];

return
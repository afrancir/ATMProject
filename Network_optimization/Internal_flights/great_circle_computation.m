% Written by Alessandro Bombelli, June 21, 2016
% Code that implements formulas that compute intersections between grat
% circle segments

clc
clear all
close all

%P = [zeros(9,1) (-4:4)'];
%P = gcwaypts(-2,-3,2,5,8);
%Q = gcwaypts(-3,-4,3,4,8);

P = [0.2 -4;0.5 -3;1 -2;0.2 -0.5;0.4 1;0.4 1.5;0.1 1.7;0.5 2;0.05 2.3];
Q = [-0.2 -4;-0.4 -3;-0.7 -2;-0.2 -0.5;-0.4 1;-0.4 1.5;-0.1 1.7;-0.5 2;-0.05 2.3;-0.08 2.4];

p = numel(P(:,1))-1;
q = numel(Q(:,1))-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Critical distances: case (a) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd1       = Haversine(P(1,:),Q(1,:),1);
cd2       = Haversine(P(end,:),Q(end,:),1);
cd_case_a = vertcat(cd1,cd2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Critical distances: case (b) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd_case_b1 = zeros((p+1)*q,1);
cont       = 1;

for i=1:p+1
    for j=1:q
        A = Q(j,:);
        B = Q(j+1,:);
        C = P(i,:);
        
        Acart    = latlon2cart(A,1);
        Bcart    = latlon2cart(B,1);
        Ccart    = latlon2cart(C,1);
        G        = cross(Acart,Bcart)/norm(cross(Acart,Bcart));
        F        = cross(Ccart,G)/norm(cross(Ccart,G));
        P1cart   = cross(F,G)/norm(cross(F,G));
        P2cart   = -P1cart;
        P1       = cart2latlon(P1cart);
        P2       = cart2latlon(P2cart);
        
        gc1      = Haversine(C,P1,1);
        gc2      = Haversine(C,P2,1);
        
        epsilon          = min([gc1 gc2]);
        cd_case_b1(cont) = epsilon;
        cont             = cont+1;
    end
end

cd_case_b2 = zeros((q+1)*p,1);
cont       = 1;

for i=1:q+1
    for j=1:p
        A = P(j,:);
        B = P(j+1,:);
        C = Q(i,:);
        
        Acart    = latlon2cart(A,1);
        Bcart    = latlon2cart(B,1);
        Ccart    = latlon2cart(C,1);
        G        = cross(Acart,Bcart)/norm(cross(Acart,Bcart));
        F        = cross(Ccart,G)/norm(cross(Ccart,G));
        P1cart   = cross(F,G)/norm(cross(F,G));
        P2cart   = -P1cart;
        P1       = cart2latlon(P1cart);
        P2       = cart2latlon(P2cart);
        
        gc1      = Haversine(C,P1,1);
        gc2      = Haversine(C,P2,1);
        
        epsilon          = min([gc1 gc2]);
        cd_case_b2(cont) = epsilon;
        cont             = cont+1;
    end
end

cd_case_b = vertcat(cd_case_b1,cd_case_b2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Critical distances: case (c) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd_case_c1 = zeros(p*q,1);
cont       = 1;

for i=1:p
    for j=1:q
        
        a0        = P(i,:);
        a1        = P(i+1,:);
        [lat,lon] = gcwaypts(a0(1),a0(2),a1(1),a1(2),2);
        a_med     = [lat(2) lon(2)];
        
        a0cart    = latlon2cart(a0,1);
        a1cart    = latlon2cart(a1,1);
        a_medcart = latlon2cart(a_med,1);
        
        nA     = cross(a0cart,a1cart)/norm(cross(a0cart,a1cart));
        e      = a_medcart/norm(a_medcart);
        theta  = pi/2;
        nB     = cos(theta)*nA+sin(theta)*(cross(e,nA))+(1-cos(theta))*(dot(e,nA))*e;
        
        c0     = Q(j,:);
        c1     = Q(j+1,:);
        
        c0cart = latlon2cart(c0,1);
        c1cart = latlon2cart(c1,1);
        nC     = cross(c0cart,c1cart)/norm(cross(c0cart,c1cart));
        
        t      = cross(nB,nC)/norm(cross(nB,nC));
        
        s1     = dot(cross(c0cart,nC),t);
        s2     = dot(cross(c1cart,nC),t);
        
        % if s1>=0 && s2<=0, the great circle from the midpoint of the great circle
        % arc intersects the other great circle arc, and the intersection point is
        % along versor t
        if s1<=0 && s2>=0
            latlon       = cart2latlon(t);
            intersection = latlon;
            epsilon      = Haversine(a_med,intersection,1);
            % if s1<=0 && s2>=0, the great circle from the midpoint of the great circle
            % arc intersects the other great circle arc, and the intersection point is
            % along versor -t
        elseif s1>=0 && s2<=0
            latlon = cart2latlon(-t);
            intersection = latlon;
            epsilon      = Haversine(a_med,intersection,1);
        else
            epsilon      = NaN;
        end
        cd_case_c1(cont) = epsilon;
        cont             = cont+1;
    end
end

cd_case_c2 = zeros(q*p,1);
cont       = 1;

for i=1:q
    for j=1:p
        
        a0        = Q(i,:);
        a1        = Q(i+1,:);
        [lat,lon] = gcwaypts(a0(1),a0(2),a1(1),a1(2),2);
        a_med     = [lat(2) lon(2)];
        
        a0cart    = latlon2cart(a0,1);
        a1cart    = latlon2cart(a1,1);
        a_medcart = latlon2cart(a_med,1);
        
        nA     = cross(a0cart,a1cart)/norm(cross(a0cart,a1cart));
        e      = a_medcart/norm(a_medcart);
        theta  = pi/2;
        nB     = cos(theta)*nA+sin(theta)*(cross(e,nA))+(1-cos(theta))*(dot(e,nA))*e;
        
        c0     = P(j,:);
        c1     = P(j+1,:);
        
        c0cart = latlon2cart(c0,1);
        c1cart = latlon2cart(c1,1);
        nC     = cross(c0cart,c1cart)/norm(cross(c0cart,c1cart));
        
        t      = cross(nB,nC)/norm(cross(nB,nC));
        
        s1     = dot(cross(c0cart,nC),t);
        s2     = dot(cross(c1cart,nC),t);
        
        % if s1>=0 && s2<=0, the great circle from the midpoint of the great circle
        % arc intersects the other great circle arc, and the intersection point is
        % along versor t
        if s1<=0 && s2>=0
            latlon       = cart2latlon(t);
            intersection = latlon;
            epsilon      = Haversine(a_med,intersection,1);
            % if s1<=0 && s2>=0, the great circle from the midpoint of the great circle
            % arc intersects the other great circle arc, and the intersection point is
            % along versor -t
        elseif s1>=0 && s2<=0
            latlon = cart2latlon(-t);
            intersection = latlon;
            epsilon      = Haversine(a_med,intersection,1);
        else
            epsilon      = NaN;
        end
        cd_case_c2(cont) = epsilon;
        cont             = cont+1;
    end
end

cd_case_c = vertcat(cd_case_c1,cd_case_c2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Assembling all critical distances %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd_all = vertcat(cd_case_a,cd_case_b,cd_case_c);
cd_all = cd_all(~isnan(cd_all));
size(cd_all)
cd_all = unique(cd_all);
size(cd_all)
[cd_all,~] = mmunique(cd_all,0.0000000001);


alpha0 = cd_all(round(numel(cd_all)/2));
tic
[FD,no_sol,sol]           = FD_computation_v02(P,Q,1,alpha0,0.0000001);
time=toc;
tic
[no_sol2,sol2,cd_all,FD2] = FD_computation_v03(P,Q,1);
time2=toc;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% a0        = [0 -2];
% a1        = [0 2];
% [lat,lon] = gcwaypts(a0(1),a0(2),a1(1),a1(2),2);
% a_med     = [lat(2) lon(2)];
% 
% a0cart    = latlon2cart(a0,1);
% a1cart    = latlon2cart(a1,1);
% a_medcart = latlon2cart(a_med,1);
% 
% nA     = cross(a0cart,a1cart)/norm(cross(a0cart,a1cart));
% e      = a_medcart/norm(a_medcart);
% theta  = pi/2;
% nB     = cos(theta)*nA+sin(theta)*(cross(e,nA))+(1-cos(theta))*(dot(e,nA))*e;
% 
% c0     = [2 -4];
% c1     = [2  0.1];
% 
% c0cart = latlon2cart(c0,1);
% c1cart = latlon2cart(c1,1);
% nC     = cross(c0cart,c1cart)/norm(cross(c0cart,c1cart));
% 
% t      = cross(nB,nC)/norm(cross(nB,nC));
% 
% s1     = dot(cross(c0cart,nC),t);
% s2     = dot(cross(c1cart,nC),t);
% 
% % if s1>=0 && s2<=0, the great circle from the midpoint of the great circle
% % arc intersects the other great circle arc, and the intersection point is
% % along versor t
% if s1<=0 && s2>=0
%     latlon       = cart2latlon(t);
%     intersection = latlon;
%     epsilon      = Haversine(a_med,intersection,1);
% % if s1<=0 && s2>=0, the great circle from the midpoint of the great circle
% % arc intersects the other great circle arc, and the intersection point is
% % along versor -t
% elseif s1>=0 && s2<=0
%     latlon = cart2latlon(-t);
%     intersection = latlon;
%     epsilon      = Haversine(a_med,intersection,1);
% else
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% A = [0 -2];
% B = [0 2];
% C = [5 0];
% 
% Acart    = latlon2cart(A,1);
% Bcart    = latlon2cart(B,1);
% Ccart    = latlon2cart(C,1);
% G        = cross(Acart,Bcart)/norm(cross(Acart,Bcart));
% F        = cross(Ccart,G)/norm(cross(Ccart,G));
% P1cart   = cross(F,G)/norm(cross(F,G));
% P2cart   = -P1cart;
% P1       = cart2latlon(P1cart);
% P2       = cart2latlon(P2cart);
% 
% gc1      = Haversine(C,P1,1);
% gc2      = Haversine(C,P2,1);
% 
% epsilon = min([gc1 gc2]);





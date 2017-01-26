function [cd_all,FD,sol,no_sol,Lmin,Lmax,Bmin,Bmax,A,sequence,Node_ID] = FD_computation_v03(P,Q,r)

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
cd_all = unique(cd_all);
[cd_all,~] = mmunique(cd_all,0.0000001);

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Computation of critical distances completed')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

found_FD    = 1;
no_sol      = [];
sol         = [];
LB_idx      = 1;
UB_idx      = numel(cd_all);

% Initial guess is the central value among all the critical distances 
alpha0_idx = round(numel(cd_all)/2);
alpha0     = cd_all(alpha0_idx);

while found_FD % keep iterating until a solution is found
    [~,~,~,~,A_list,~] = free_space_shape_v02(P,Q,r,alpha0);
    startpoint         = 1;
    endpoint           = length(A_list);
    [is_solution,~]    = DFS(startpoint,endpoint,A_list);
    
    % no solution found
    if is_solution == 0
        LB_idx     = alpha0_idx;
        no_sol     = vertcat(no_sol,[alpha0_idx alpha0]);
        range      = alpha0_idx+1:UB_idx;
        if numel(range)==1
            found_FD = 0;
            [~,~,~,~,A_list] = free_space_shape_v02(P,Q,r,cd_all(alpha0_idx+1));
            startpoint       = 1;
            endpoint         = length(A_list);
            is_solution      = DFS(startpoint,endpoint,A_list);
            if is_solution==0
                sol    = vertcat(sol,[alpha0_idx+2 cd_all(alpha0_idx+2)]);
            else
                sol    = vertcat(sol,[alpha0_idx+1 cd_all(alpha0_idx+1)]);
            end
        else
            alpha0_idx = range(round(numel(range)/2));
            alpha0     = cd_all(alpha0_idx);
        end
        % solution found
    else
        UB_idx     = alpha0_idx;
        sol        = vertcat(sol,[alpha0_idx alpha0]);
        range      = LB_idx:alpha0_idx-1;
        if numel(range)==1
            found_FD = 0;
            [~,~,~,~,A_list] = free_space_shape_v02(P,Q,r,cd_all(alpha0_idx-1));
            startpoint       = 1;
            endpoint         = length(A_list);
            is_solution      = DFS(startpoint,endpoint,A_list);
            if is_solution==0
                sol    = vertcat(sol,[alpha0_idx cd_all(alpha0_idx)]);
            else
                sol    = vertcat(sol,[alpha0_idx-1 cd_all(alpha0_idx-1)]);
            end
        else
            alpha0_idx = range(round(numel(range)/2));
            alpha0     = cd_all(alpha0_idx);
        end
    end
end

FD = cd_all(min(sol(:,1)));
[Lmin,Lmax,Bmin,Bmax,A_list,Node_ID] = free_space_shape_v02(P,Q,r,FD);
startpoint                           = 1;
endpoint                             = length(A_list);
[~,A]                                = DFS(startpoint,endpoint,A_list);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idx           = find(A(1,:)==endpoint);
N_steps       = A(3,idx(1));
sequence      = zeros(N_steps,1);
sequence(end) = endpoint;

for i=1:N_steps-1
    idx = find(A(1,:)==sequence(N_steps-i+1));
    for j=1:numel(idx)
        if A(3,idx(j)) == N_steps-i+1
            sequence(N_steps-i) = A(2,idx(j));
            break;
        else
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return
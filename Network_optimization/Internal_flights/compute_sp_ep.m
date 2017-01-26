function [Overlap_points,P_length,Q_length] = compute_sp_ep(Merged_Sequences,P,Q,alpha,FD,Bmin,Bmax,Lmin,Lmax)

Overlap_points = zeros(length(Merged_Sequences(:,1)),16);
N_f             = 20;

p = length(P(:,1))-1;
q = length(Q(:,1))-1;

P_length = zeros(p,1);
Q_length = zeros(q,1);

% Determine length of each edge for polygonal curve P
for i=1:p
    P_length(i,1) = Haversine(P(i,:),P(i+1,:),1);
end

% Determine length of each edge for polygonal curve Q
for i=1:q
    Q_length(i,1) = Haversine(Q(i,:),Q(i+1,:),1);
end

for i=1:length(Merged_Sequences(:,1))
    
    Type_sp = Merged_Sequences(i,1);
    r_sp    = Merged_Sequences(i,2);
    c_sp    = Merged_Sequences(i,3);
    Type_ep = Merged_Sequences(i,5);
    r_ep    = Merged_Sequences(i,6);
    c_ep    = Merged_Sequences(i,7);
    
    %%%%%%%%%%%%%%%%%%%%%%
    %%% Starting point %%%
    %%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%
        %%% Corner node %%%
        %%%%%%%%%%%%%%%%%%%
        if Type_sp == 1
            % First row of free space
            if r_sp == 1
                if c_sp == 1
                    sp_Q             = Q(r_sp,:);
                    sp_P             = P(c_sp,:);
                    sp_best_f_Q      = 0;
                    sp_best_f_P      = 0;
                    sp_unit_square_Q = r_sp;
                    sp_unit_square_P = c_sp;
                elseif c_sp == p+1
                    sp_Q             = Q(r_sp,:);
                    sp_P             = P(c_sp,:);
                    sp_best_f_Q      = 0;
                    sp_best_f_P      = 1;
                    sp_unit_square_Q = r_sp;
                    sp_unit_square_P = p;
                else
                    sp_Q             = Q(r_sp,:);
                    sp_P             = P(c_sp,:);
                    sp_best_f_Q      = 0;
                    sp_best_f_P      = 0;
                    sp_unit_square_Q = r_sp;
                    sp_unit_square_P = c_sp;
                end
            % Last row of free space    
            elseif r_sp == q+1
                if c_sp == 1
                    sp_Q             = Q(r_sp,:);
                    sp_P             = P(c_sp,:);
                    sp_best_f_Q      = 1;
                    sp_best_f_P      = 0;
                    sp_unit_square_Q = q;
                    sp_unit_square_P = c_sp;
                elseif c_sp == p+1
                    sp_Q             = Q(r_sp,:);
                    sp_P             = P(c_sp,:);
                    sp_best_f_Q      = 1;
                    sp_best_f_P      = 1;
                    sp_unit_square_Q = q;
                    sp_unit_square_P = p;
                else
                    sp_Q             = Q(r_sp,:);
                    sp_P             = P(c_sp,:);
                    sp_best_f_Q      = 1;
                    sp_best_f_P      = 0;
                    sp_unit_square_Q = q;
                    sp_unit_square_P = c_sp;
                end
            % Internal row of the free space diagram    
            else
                if c_sp == 1
                    sp_Q             = Q(r_sp,:);
                    sp_P             = P(c_sp,:);
                    sp_best_f_Q      = 0;
                    sp_best_f_P      = 0;
                    sp_unit_square_Q = r_sp;
                    sp_unit_square_P = c_sp;
                elseif c_sp == p+1
                    sp_Q             = Q(r_sp,:);
                    sp_P             = P(c_sp,:);
                    sp_best_f_Q      = 0;
                    sp_best_f_P      = 1;
                    sp_unit_square_Q = r_sp;
                    sp_unit_square_P = p;
                else
                    sp_Q             = Q(r_sp,:);
                    sp_P             = P(c_sp,:);
                    sp_best_f_Q      = 0;
                    sp_best_f_P      = 0;
                    sp_unit_square_Q = r_sp;
                    sp_unit_square_P = c_sp;
                end
            end
        %%%%%%%%%%%%%%%%%%%%%%%
        %%% Horizontal node %%%
        %%%%%%%%%%%%%%%%%%%%%%%
        elseif Type_sp == 2
            % Upper boundary of the free space diagram. The starting point
            % is on the upper boundary itself
            if r_sp == 1
                sp_Q             = Q(r_sp,:);
                P_tail           = P(c_sp,:);
                P_head           = P(c_sp+1,:);
                P_tail_cart      = latlon2cart(P_tail,1);
                P_head_cart      = latlon2cart(P_head,1);
                sp_best_f_Q      = 0;
                sp_best_f_P      = Bmin(r_sp,c_sp);
                DP               = P_length(c_sp);
                P_cart_sp        = P_tail_cart*sin((1-sp_best_f_P)*DP)/sin(DP)+P_head_cart*sin(sp_best_f_P*DP)/sin(DP);
                sp_P             = cart2latlon(P_cart_sp);
                sp_unit_square_Q = r_sp;
                sp_unit_square_P = c_sp;
            % Otherwise, the unitary square above needs to be inspected    
            else
                Q_tail           = Q(r_sp-1,:);
                Q_head           = Q(r_sp,:);
                P_tail           = P(c_sp,:);
                P_head           = P(c_sp+1,:);
                Q_tail_cart      = latlon2cart(Q_tail,1);
                Q_head_cart      = latlon2cart(Q_head,1);
                P_tail_cart      = latlon2cart(P_tail,1);
                P_head_cart      = latlon2cart(P_head,1);
                DQ               = Q_length(r_sp-1);
                DP               = P_length(c_sp);
                f_vec_P          = linspace(0,Bmax(r_sp,c_sp),N_f);
                f_vec_Q          = linspace(0,1,N_f);
                dist_mat         = zeros(numel(f_vec_Q),numel(f_vec_P));
                dist_head        = zeros(size(dist_mat));
                for ii=1:numel(f_vec_Q)
                    Q_cart   = Q_tail_cart*sin((1-f_vec_Q(ii))*DQ)/sin(DQ)+Q_head_cart*sin(f_vec_Q(ii)*DQ)/sin(DQ);
                    Q_ll     = cart2latlon(Q_cart);
                    for jj=1:numel(f_vec_P)
                        P_cart           = P_tail_cart*sin((1-f_vec_P(jj))*DP)/sin(DP)+P_head_cart*sin(f_vec_P(jj)*DP)/sin(DP);
                        P_ll             = cart2latlon(P_cart);
                        dist_mat(ii,jj)  = Haversine(P_ll,Q_ll,1);
                        dist_head(ii,jj) = Haversine(P_ll,P_head,1)+Haversine(Q_ll,Q_head,1);
                    end
                end
                [r,c]    = find(dist_mat<=alpha*FD);
                dist_head_vec = zeros(numel(r),1);
                for ii=1:numel(r)
                    dist_head_vec(ii) = dist_head(r(ii),c(ii));
                end
                [~,idx]          = max(dist_head_vec);
                sp_best_f_Q      = f_vec_Q(r(idx));
                sp_best_f_P      = f_vec_P(c(idx));
                P_cart_sp        = P_tail_cart*sin((1-sp_best_f_P)*DP)/sin(DP)+P_head_cart*sin(sp_best_f_P*DP)/sin(DP);
                Q_cart_sp        = Q_tail_cart*sin((1-sp_best_f_Q)*DQ)/sin(DQ)+Q_head_cart*sin(sp_best_f_Q*DQ)/sin(DQ);
                sp_P             = cart2latlon(P_cart_sp);
                sp_Q             = cart2latlon(Q_cart_sp);
                sp_unit_square_Q = r_sp-1;
                sp_unit_square_P = c_sp;
            end
        %%%%%%%%%%%%%%%%%%%%%
        %%% Vertical node %%%
        %%%%%%%%%%%%%%%%%%%%%    
        else
            % Left boundary of the free space diagram. The starting point
            % is on the left boundary itself
            if c_sp == 1
                sp_P             = P(c_sp,:);
                Q_tail           = Q(r_sp,:);
                Q_head           = Q(r_sp+1,:);
                Q_tail_cart      = latlon2cart(Q_tail,1);
                Q_head_cart      = latlon2cart(Q_head,1);
                sp_best_f_P      = 0;
                sp_best_f_Q      = Lmin(r_sp,c_sp);
                DQ               = Q_length(r_sp);
                Q_cart_sp        = Q_tail_cart*sin((1-sp_best_f_Q)*DQ)/sin(DQ)+Q_head_cart*sin(sp_best_f_Q*DQ)/sin(DQ);
                sp_Q             = cart2latlon(Q_cart_sp);
                sp_unit_square_Q = r_sp;
                sp_unit_square_P = c_sp;
            else
                Q_tail           = Q(r_sp,:);
                Q_head           = Q(r_sp+1,:);
                P_tail           = P(c_sp-1,:);
                P_head           = P(c_sp,:);
                Q_tail_cart      = latlon2cart(Q_tail,1);
                Q_head_cart      = latlon2cart(Q_head,1);
                P_tail_cart      = latlon2cart(P_tail,1);
                P_head_cart      = latlon2cart(P_head,1);
                DQ               = Q_length(r_sp);
                DP               = P_length(c_sp-1);
                f_vec_P          = linspace(0,1,N_f);
                f_vec_Q          = linspace(0,Lmax(r_sp,c_sp),N_f);
                dist_mat         = zeros(numel(f_vec_Q),numel(f_vec_P));
                dist_head        = zeros(size(dist_mat));
                for ii=1:numel(f_vec_Q)
                    Q_cart   = Q_tail_cart*sin((1-f_vec_Q(ii))*DQ)/sin(DQ)+Q_head_cart*sin(f_vec_Q(ii)*DQ)/sin(DQ);
                    Q_ll     = cart2latlon(Q_cart);
                    for jj=1:numel(f_vec_P)
                        P_cart           = P_tail_cart*sin((1-f_vec_P(jj))*DP)/sin(DP)+P_head_cart*sin(f_vec_P(jj)*DP)/sin(DP);
                        P_ll             = cart2latlon(P_cart);
                        dist_mat(ii,jj)  = Haversine(P_ll,Q_ll,1);
                        dist_head(ii,jj) = Haversine(P_ll,P_head,1)+Haversine(Q_ll,Q_head,1);
                    end
                end
                [r,c]    = find(dist_mat<=alpha*FD);
                dist_head_vec = zeros(numel(r),1);
                for ii=1:numel(r)
                    dist_head_vec(ii) = dist_head(r(ii),c(ii));
                end
                [~,idx]          = max(dist_head_vec);
                sp_best_f_Q      = f_vec_Q(r(idx));
                sp_best_f_P      = f_vec_P(c(idx));
                P_cart_sp        = P_tail_cart*sin((1-sp_best_f_P)*DP)/sin(DP)+P_head_cart*sin(sp_best_f_P*DP)/sin(DP);
                Q_cart_sp        = Q_tail_cart*sin((1-sp_best_f_Q)*DQ)/sin(DQ)+Q_head_cart*sin(sp_best_f_Q*DQ)/sin(DQ);
                sp_P             = cart2latlon(P_cart_sp);
                sp_Q             = cart2latlon(Q_cart_sp);
                sp_unit_square_Q = r_sp;
                sp_unit_square_P = c_sp-1;
            end
        end
    
    
    %%%%%%%%%%%%%%%%
    %%% Endpoint %%%
    %%%%%%%%%%%%%%%%
    
        %%%%%%%%%%%%%%%%%%%
        %%% Corner node %%%
        %%%%%%%%%%%%%%%%%%%
        if Type_ep == 1
            % First row of free space
            if r_ep == 1
                if c_ep == 1
                    ep_Q             = Q(1,:);
                    ep_P             = P(c_ep,:);
                    ep_best_f_Q      = 0;
                    ep_best_f_P      = 0;
                    ep_unit_square_Q = 1;
                    ep_unit_square_P = c_ep;
                elseif c_ep == p+1
                    ep_Q             = Q(1,:);
                    ep_P             = P(c_ep,:);
                    ep_best_f_Q      = 0;
                    ep_best_f_P      = 1;
                    ep_unit_square_Q = 1;
                    ep_unit_square_P = p;
                else
                    ep_Q             = Q(1,:);
                    ep_P             = P(c_ep,:);
                    ep_best_f_Q      = 0;
                    ep_best_f_P      = 0;
                    ep_unit_square_Q = 1;
                    ep_unit_square_P = c_ep;
                end
            % Last row of free space    
            elseif r_ep == q+1
                if c_ep == 1
                    ep_Q             = Q(end,:);
                    ep_P             = P(c_ep,:);
                    ep_best_f_Q      = 1;
                    ep_best_f_P      = 0;
                    ep_unit_square_Q = q;
                    ep_unit_square_P = c_ep;
                elseif c_ep == p+1
                    ep_Q             = Q(end,:);
                    ep_P             = P(c_ep,:);
                    ep_best_f_Q      = 1;
                    ep_best_f_P      = 1;
                    ep_unit_square_Q = q;
                    ep_unit_square_P = p;
                else
                    ep_Q             = Q(end,:);
                    ep_P             = P(c_ep,:);
                    ep_best_f_Q      = 1;
                    ep_best_f_P      = 0;
                    ep_unit_square_Q = q;
                    ep_unit_square_P = c_ep;
                end
            % Internal row of the free space diagram    
            else
                if c_ep == 1
                    ep_Q             = Q(r_ep,:);
                    ep_P             = P(c_ep,:);
                    ep_best_f_Q      = 0;
                    ep_best_f_P      = 0;
                    ep_unit_square_Q = r_ep;
                    ep_unit_square_P = c_ep;
                elseif c_ep == p+1
                    ep_Q             = Q(r_ep,:);
                    ep_P             = P(c_ep,:);
                    ep_best_f_Q      = 0;
                    ep_best_f_P      = 1;
                    ep_unit_square_Q = r_ep;
                    ep_unit_square_P = p;
                else
                    ep_Q             = Q(r_ep,:);
                    ep_P             = P(c_ep,:);
                    ep_best_f_Q      = 0;
                    ep_best_f_P      = 0;
                    ep_unit_square_Q = r_ep;
                    ep_unit_square_P = c_ep;
                end
            end
        %%%%%%%%%%%%%%%%%%%%%%%
        %%% Horizontal node %%%
        %%%%%%%%%%%%%%%%%%%%%%%
        elseif Type_ep == 2
            % Lower boundary of the free space diagram. The  endpoint
            % is on the lower boundary itself
            if r_ep == q+1
                ep_Q             = Q(end,:);
                P_tail           = P(c_ep,:);
                P_head           = P(c_ep+1,:);
                P_tail_cart      = latlon2cart(P_tail,1);
                P_head_cart      = latlon2cart(P_head,1);
                ep_best_f_Q      = 1;
                ep_best_f_P      = Bmax(r_ep,c_ep);
                DP               = P_length(c_ep);
                P_cart_ep        = P_tail_cart*sin((1-ep_best_f_P)*DP)/sin(DP)+P_head_cart*sin(ep_best_f_P*DP)/sin(DP);
                ep_P             = cart2latlon(P_cart_ep);
                ep_unit_square_Q = q;
                ep_unit_square_P = c_sp;
            % Otherwise, the unitary square below needs to be inspected    
            else
                Q_tail           = Q(r_ep,:);
                Q_head           = Q(r_ep+1,:);
                P_tail           = P(c_ep,:);
                P_head           = P(c_ep+1,:);
                Q_tail_cart      = latlon2cart(Q_tail,1);
                Q_head_cart      = latlon2cart(Q_head,1);
                P_tail_cart      = latlon2cart(P_tail,1);
                P_head_cart      = latlon2cart(P_head,1);
                DQ               = Q_length(r_ep);
                DP               = P_length(c_ep);
                f_vec_P          = linspace(Bmin(r_ep,c_ep),1,N_f);
                f_vec_Q          = linspace(0,1,N_f);
                dist_mat         = zeros(numel(f_vec_Q),numel(f_vec_P));
                dist_tail        = zeros(size(dist_mat));
                for ii=1:numel(f_vec_Q)
                    Q_cart   = Q_tail_cart*sin((1-f_vec_Q(ii))*DQ)/sin(DQ)+Q_head_cart*sin(f_vec_Q(ii)*DQ)/sin(DQ);
                    Q_ll     = cart2latlon(Q_cart);
                    for jj=1:numel(f_vec_P)
                        P_cart           = P_tail_cart*sin((1-f_vec_P(jj))*DP)/sin(DP)+P_head_cart*sin(f_vec_P(jj)*DP)/sin(DP);
                        P_ll             = cart2latlon(P_cart);
                        dist_mat(ii,jj)  = Haversine(P_ll,Q_ll,1);
                        dist_tail(ii,jj) = Haversine(P_tail,P_ll,1)+Haversine(Q_tail,Q_ll,1);
                    end
                end
                [r,c]    = find(dist_mat<=alpha*FD);
                dist_tail_vec = zeros(numel(r),1);
                for ii=1:numel(r)
                    dist_tail_vec(ii) = dist_tail(r(ii),c(ii));
                end
                [~,idx]          = max(dist_tail_vec);
                ep_best_f_Q      = f_vec_Q(r(idx));
                ep_best_f_P      = f_vec_P(c(idx));
                P_cart_ep        = P_tail_cart*sin((1-ep_best_f_P)*DP)/sin(DP)+P_head_cart*sin(ep_best_f_P*DP)/sin(DP);
                Q_cart_ep        = Q_tail_cart*sin((1-ep_best_f_Q)*DQ)/sin(DQ)+Q_head_cart*sin(ep_best_f_Q*DQ)/sin(DQ);
                ep_P             = cart2latlon(P_cart_ep);
                ep_Q             = cart2latlon(Q_cart_ep);
                ep_unit_square_Q = r_ep;
                ep_unit_square_P = c_ep;
            end
        else
        %%%%%%%%%%%%%%%%%%%%%
        %%% Vertical node %%%
        %%%%%%%%%%%%%%%%%%%%%
        % Right boundary of the free space diagram. The endpoint
        % is on the right boundary itself
            if c_ep == p+1
                ep_P             = P(c_ep,:);
                Q_tail           = Q(r_ep,:);
                Q_head           = Q(r_ep+1,:);
                Q_tail_cart      = latlon2cart(Q_tail,1);
                Q_head_cart      = latlon2cart(Q_head,1);
                ep_best_f_P      = 1;
                ep_best_f_Q      = Lmax(r_ep,c_ep);
                DQ               = Q_length(r_ep);
                Q_cart_ep        = Q_tail_cart*sin((1-ep_best_f_Q)*DQ)/sin(DQ)+Q_head_cart*sin(ep_best_f_Q*DQ)/sin(DQ);
                ep_Q             = cart2latlon(Q_cart_ep);
                sp_unit_square_Q = r_ep;
                sp_unit_square_P = p;
            else
                Q_tail           = Q(r_ep,:);
                Q_head           = Q(r_ep+1,:);
                P_tail           = P(c_ep,:);
                P_head           = P(c_ep+1,:);
                Q_tail_cart      = latlon2cart(Q_tail,1);
                Q_head_cart      = latlon2cart(Q_head,1);
                P_tail_cart      = latlon2cart(P_tail,1);
                P_head_cart      = latlon2cart(P_head,1);
                DQ               = Q_length(r_ep);
                DP               = P_length(c_ep);
                f_vec_P          = linspace(0,1,N_f);
                f_vec_Q          = linspace(Lmin(r_ep,c_ep),1,N_f);
                dist_mat         = zeros(numel(f_vec_Q),numel(f_vec_P));
                dist_tail        = zeros(size(dist_mat));
                for ii=1:numel(f_vec_Q)
                    Q_cart   = Q_tail_cart*sin((1-f_vec_Q(ii))*DQ)/sin(DQ)+Q_head_cart*sin(f_vec_Q(ii)*DQ)/sin(DQ);
                    Q_ll     = cart2latlon(Q_cart);
                    for jj=1:numel(f_vec_P)
                        P_cart           = P_tail_cart*sin((1-f_vec_P(jj))*DP)/sin(DP)+P_head_cart*sin(f_vec_P(jj)*DP)/sin(DP);
                        P_ll             = cart2latlon(P_cart);
                        dist_mat(ii,jj)  = Haversine(P_ll,Q_ll,1);
                        dist_tail(ii,jj) = Haversine(P_tail,P_ll,1)+Haversine(Q_tail,Q_ll,1);
                    end
                end
                [r,c]    = find(dist_mat<=alpha*FD);
                dist_tail_vec = zeros(numel(r),1);
                for ii=1:numel(r)
                    dist_tail_vec(ii) = dist_tail(r(ii),c(ii));
                end
                [~,idx]          = max(dist_tail_vec);
                ep_best_f_Q      = f_vec_Q(r(idx));
                ep_best_f_P      = f_vec_P(c(idx));
                P_cart_ep        = P_tail_cart*sin((1-ep_best_f_P)*DP)/sin(DP)+P_head_cart*sin(ep_best_f_P*DP)/sin(DP);
                Q_cart_ep        = Q_tail_cart*sin((1-ep_best_f_Q)*DQ)/sin(DQ)+Q_head_cart*sin(ep_best_f_Q*DQ)/sin(DQ);
                ep_P             = cart2latlon(P_cart_ep);
                ep_Q             = cart2latlon(Q_cart_ep);
                ep_unit_square_Q = r_ep;
                ep_unit_square_P = c_ep;
            end
        end
        
        Overlap_points(i,:) = horzcat(sp_unit_square_P,sp_best_f_P,sp_P,sp_unit_square_Q,sp_best_f_Q,sp_Q,...
                                       ep_unit_square_P,ep_best_f_P,ep_P,ep_unit_square_Q,ep_best_f_Q,ep_Q);
        
end

    return
function idx_to_keep = RamerDouglasPeucker(p,eps)

% INPUTS:
% - p: Nx2 matrix that defines the original trajectory. Each row is a
% (lat,lon) pair
% - eps: distance threshold used to undersample the trajectory. See
% References to understand how epsilon is used in the algorithm
% OUTPUTS:
% idx_to_keep: Nx1 vector, where each element can be either 1 or 0. All the
% 1's define the indices of the original trajectory that should also be
% kept for the reduced trajectory. The initial and final element of
% idx_to_keep are 1 be definition

% idx_to_keep is a variable that is spans multiple functions
% (RamerDouglasPeucker and RDP)
idx_to_keep  = zeros(length(p(:,1)),1);
% Initial index
idx_in       = 1;
% Final index
idx_fin      = length(p(:,1));
% Initial trajectory
p_unfiltered = p;

RDP(p_unfiltered,eps,idx_in,idx_fin)

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Recursive function %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
function RDP(p_unfiltered,eps,idx_in,idx_fin)

if idx_fin-idx_in == 1 % no points in between
    idx_to_keep(idx_in)  = 1;
    idx_to_keep(idx_fin) = 1;
else
    distance = zeros(idx_fin-idx_in-1,1);
    A        = p_unfiltered(idx_in,:);
    B        = p_unfiltered(idx_fin,:);
    for j=1:idx_fin-idx_in-1
        C           = p_unfiltered(idx_in+j,:);
        Acart       = latlon2cart(A,1);
        Bcart       = latlon2cart(B,1);
        Ccart       = latlon2cart(C,1);
        G           = cross(Acart,Bcart)/norm(cross(Acart,Bcart));
        F        = cross(Ccart,G)/norm(cross(Ccart,G));
        P1cart   = cross(F,G)/norm(cross(F,G));
        P2cart   = -P1cart;
        P1       = cart2latlon(P1cart);
        P2       = cart2latlon(P2cart);
        
        gc1      = Haversine(C,P1,1);
        gc2      = Haversine(C,P2,1);
        distance(j) = min([gc1 gc2]);
    end
    [dmax,idx_dmax] = max(distance);
    idx_dmax        = idx_dmax+idx_in;
    if dmax <= eps
        idx_to_keep(idx_in)  = 1;
        idx_to_keep(idx_fin) = 1;
    else
        RDP(p_unfiltered,eps,idx_in,idx_dmax);
        RDP(p_unfiltered,eps,idx_dmax,idx_fin);
    end
end

end
end
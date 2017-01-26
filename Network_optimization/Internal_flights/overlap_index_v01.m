function [P_length,Q_length,dist_vector,NODES] = overlap_index_v01(P,Q,sequence,Node_ID,Lmin,Lmax,Bmin,Bmax)

p = length(P(:,1))-1; % number of edges of first polygonal curve
q = length(Q(:,1))-1; % number of edges of second polygonal curve
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

% For each node of the sequence, retrieve the two associated points
% (one on P, one on Q), and compute their distance
dist_vector = zeros(numel(sequence),1);
NODES       = zeros(numel(sequence),4);
N_p         = 50;

for i=1:numel(sequence)
    ID    = sequence(i);
    BLOCK = ceil(ID/(3*p+2));
    POS   = ID-(BLOCK-1)*(3*p+2);
    node_type = Node_ID(sequence(i),1);
    if node_type == 1 % corner node
        R = ceil(ID/(3*p+2));
        C = ceil(POS/2);
        point_on_P = P(C,:);
        point_on_Q = Q(R,:);
        d          = Haversine(point_on_P,point_on_Q,1);
    elseif node_type == 2 % horizontal edge node
        R = ceil(ID/(3*p+2));
        C = ceil(POS/2);
        f1 = Bmin(R,C);
        f2 = Bmax(R,C);
        tail            = latlon2cart(P(C,:),1);
        head            = latlon2cart(P(C+1,:),1);
        D               = P_length(C);
        point_on_P_cart1 = tail*sin((1-f1)*D)/sin(D)+head*sin(f1*D)/sin(D);
        point_on_P1      = cart2latlon(point_on_P_cart1);
        point_on_P_cart2 = tail*sin((1-f2)*D)/sin(D)+head*sin(f2*D)/sin(D);
        point_on_P2      = cart2latlon(point_on_P_cart2);
        [lat,lon] = gcwaypts(point_on_P1(1),point_on_P1(2),point_on_P2(1),point_on_P2(2),N_p);
        d_vec     = zeros(numel(lat),1);
        for j=1:numel(lat)
            d_vec(j) = Haversine([lat(j) lon(j)],point_on_Q,1);
        end
        d = min(d_vec);
    else % vertical edge node
        R = ceil(ID/(3*p+2));
        C = ID-(BLOCK-1)*(3*p+2)-(2*p+1);
        f1 = Lmin(R,C);
        f2 = Lmax(R,C);
        tail            = latlon2cart(Q(R,:),1);
        head            = latlon2cart(Q(R+1,:),1);
        D               = Q_length(R);
        point_on_Q_cart1 = tail*sin((1-f1)*D)/sin(D)+head*sin(f1*D)/sin(D);
        point_on_Q1      = cart2latlon(point_on_Q_cart1);
        point_on_Q_cart2 = tail*sin((1-f2)*D)/sin(D)+head*sin(f2*D)/sin(D);
        point_on_Q2      = cart2latlon(point_on_Q_cart2);
        [lat,lon] = gcwaypts(point_on_Q1(1),point_on_Q1(2),point_on_Q2(1),point_on_Q2(2),N_p);
        d_vec     = zeros(numel(lat),1);
        for j=1:numel(lat)
            d_vec(j) = Haversine([lat(j) lon(j)],point_on_P,1);
        end
        d = min(d_vec);
    end
    
    NODES(i,1)     = i;
    NODES(i,2)     = ID;
    NODES(i,3)     = node_type;
    NODES(i,4)     = R;
    NODES(i,5)     = C;
    dist_vector(i) = d;
    
end

return
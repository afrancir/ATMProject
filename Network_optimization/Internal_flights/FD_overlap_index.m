function Overlap_points = FD_overlap_index(P,Q,alpha)

angle     = 1;
good_idx  = heading_filtering(P,angle);
idx       = find(good_idx);
good_idx2 = heading_filtering(Q,angle);
idx2      = find(good_idx2);

P_red    = P(idx,:);
Q_red    = Q(idx2,:);

[FD,~,~] = FD_computation_v02(P_red,Q_red,1,0.05,1);

Re    = 1;
[Lmin,Lmax,Bmin,Bmax,A_list,Node_ID_eps] = free_space_shape_v02(P_red,Q_red,Re,alpha*FD);

p = length(P_red(:,1))-1;
q = length(Q_red(:,1))-1;
N = (2*p+1)*(q+1)+(p+1)*q;

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

valid_overlap_sequences_ID     = 0;
valid_overlap_sequences_int_ID = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 1: identify the active nodes in the free space diagram. With %%%
%%% active ndoes, we mean all the nodes that characterizes feasible   %%%
%%% regions within the free space diagram itself                      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

is_defined = zeros(N,1);

for i=1:N
    
    BLOCK = ceil(i/(3*p+2));
    POS   = i-(BLOCK-1)*(3*p+2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% First block of the free space diagram %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if BLOCK == 1
        % Either corner or horizontal node
        if POS<=2*p+1
            % Corner node
            if mod(POS,2)~=0
                
                % First corner node of the block
                if POS == 1
                    if Bmin(1,1) == 0 || Lmin(1,1) == 0
                        is_defined(i) = 1;
                    else
                    end
                    % Last corner node of the block
                elseif POS == 2*p+1
                    if Bmin(1,end) == 0 || Lmin(1,end) == 0
                        is_defined(i) = 1;
                    else
                    end
                    % Internal corner node
                else
                    R     = ceil(i/(3*p+2));
                    C     = ceil(POS/2);
                    if Bmax(R,C-1) == 1 || Bmin(R,C) == 0 || Lmin(R,C) == 0
                        is_defined(i) = 1;
                    else
                    end
                end
                % Horizontal node
            else
                R     = ceil(i/(3*p+2));
                C     = ceil(POS/2);
                if  Bmin(R,C) >=0
                    is_defined(i) = 1;
                else
                end
            end
            % Vertical node
        else
            R     = ceil(i/(3*p+2));
            C     = POS-(2*p+1);
            if  Lmin(R,C) >=0
                is_defined(i) = 1;
            else
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Last block of the free space diagram %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif BLOCK == q+1
        
        if mod(POS,2)~=0
            
            % First corner node of the block
            if POS == 1
                if Bmax(end,1) == 1 || Lmax(end,1) == 1
                    is_defined(i) = 1;
                else
                end
                % Last corner node of the block
            elseif POS == 2*p+1
                if Bmax(end,end) == 1 || Lmax(end,end) == 1
                    is_defined(i) = 1;
                else
                end
                % Internal corner node
            else
                R     = ceil(i/(3*p+2));
                C     = ceil(POS/2);
                if Bmax(R,C-1) == 1 || Bmin(R,C) == 0 || Lmax(R-1,C) == 1
                    is_defined(i) = 1;
                else
                end
            end
            % Horizontal node
        else
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Intermediate blocks of the free space diagram %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        % Either corner or horizontal node
        if POS<=2*p+1
            % Corner node
            if mod(POS,2)~=0
                
                R     = ceil(i/(3*p+2));
                C     = ceil(POS/2);
                
                % First corner node of the block
                if POS == 1
                    if Bmin(R,1) == 0 || Lmax(R-1,1) == 1 || Lmin(R,1) == 0
                        is_defined(i) = 1;
                    else
                    end
                    % Last corner node of the block
                elseif POS == 2*p+1
                    if Bmax(R,end) == 1 || Lmax(R-1,end) == 1 || Lmin(R,end) == 0
                        is_defined(i) = 1;
                    else
                    end
                    % Internal corner node
                else
                    if Bmax(R,C-1) == 1 || Bmin(R,C) == 0 || Lmax(R-1,C) == 1 || Lmin(R,C) == 0
                        is_defined(i) = 1;
                    else
                    end
                end
                % Horizontal node
            else
                R     = ceil(i/(3*p+2));
                C     = ceil(POS/2);
                if  Bmin(R,C) >=0
                    is_defined(i) = 1;
                else
                end
            end
            % Vertical node
        else
            R     = ceil(i/(3*p+2));
            C     = POS-(2*p+1);
            if  Lmin(R,C) >=0
                is_defined(i) = 1;
            else
            end
        end
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 2: once the active nodes have been identified, determine %%%
%%% sub-regions of the free space diagram that are likely to      %%%
%%% show an overlapping region. Note that this step is not        %%%
%%% necessary, since we could be using Depth First Search with    %%%
%%% the initial set of active nodes, but it helps reducing the    %%%
%%% computational time. It is based on a heuristic which can be   %%%
%%% improved in future versions of the algorithm                  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

good_nodes = find(is_defined);

idx_c    = unique(Node_ID_eps(good_nodes,3));
dummy    = idx_c(2:end)-idx_c(1:end-1);
find_gap = find(dummy~=1);

if isempty(find_gap)
else
    ID_mat = zeros(numel(find_gap)+1,2);
    for i=1:numel(find_gap)+1
        if i == 1
            if find_gap(i) == 1
                c_min = idx_c(1);
                c_max = idx_c(1);
            else
                c_min = idx_c(1);
                c_max = idx_c(find_gap(i));
            end
        elseif i == numel(find_gap)+1
            if find_gap(i-1) == numel(dummy)
                c_min = idx_c(end);
                c_max = idx_c(end);
            else
                c_min = idx_c(find_gap(i-1)+1);
                c_max = idx_c(end);
            end
        else
            c_min = idx_c(find_gap(i-1)+1);
            c_max = idx_c(find_gap(i));
        end
        % Indices of the Node_ID_eps(good_nodes,:) matrix that satisfy the
        % lower boundary constraint
        idx_c_min   = find(Node_ID_eps(good_nodes,3)>=c_min);
        % Indices of the Node_ID_eps(good_nodes,:) matrix that satisfy the
        % upper boundary constraint
        idx_c_max   = find(Node_ID_eps(good_nodes,3)<=c_max);
        % Get the common indices, and use them w.r.t. the first vector
        % [values,idx_a,idx_b] = intersect(a,b)
        [~,idx_all] = intersect(idx_c_min,idx_c_max);
        [~,idx_min] = min(Node_ID_eps(good_nodes(idx_c_min(idx_all)),4));
        [~,idx_max] = max(Node_ID_eps(good_nodes(idx_c_min(idx_all)),4));
        
        % User-defined parameter. Given the current sequence, we might be
        % interested in adding some previous starting indices, and some
        % following ending indices, to catch portions of the current
        % sequence we might lose otherwise. Note that this potential issue
        % is caused by the heuristic process mentioned in the "STEP 2"
        % description. If we get rid of that, we are sure 100% we are not
        % losing information. Experiments showed that in almost all cases,
        % information are not lost anyway, with a considerable decrease in
        % the computational time.
        eta = 10;
        
        if idx_c_min(idx_all(idx_min))-eta < 1
            ID_mat(i,1) = Node_ID_eps(good_nodes(1),4);
        else
            ID_mat(i,1) = Node_ID_eps(good_nodes(idx_c_min(idx_all(idx_min))-eta),4);
        end
        if idx_c_min(idx_all(idx_max))+eta > numel(good_nodes)
            ID_mat(i,2) = Node_ID_eps(good_nodes(end),4);
        else
            ID_mat(i,2) = Node_ID_eps(good_nodes(idx_c_min(idx_all(idx_max))+eta),4);
        end
        
        
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 3: use a Depth First Search algorithm to solve every sequence %%%
%%% and identify the longest sub-sequence with a solution. Start with  %%%
%%% the smallest initial point and the largest endpoint, and shrink    %%%
%%% the interval until a solution has been found. Note that we are     %%%
%%% using the property that the graph is directed to update the        %%%
%%% vectors defining the starting point and endpoint for each          %%%
%%% sub-sequence                                                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cont_overlap_sequences = 1;
% We initialize the matrix with a sufficiently high number of rows. Later
% on, we will only save the non-zero rows (if any)
overlap_sequences      = zeros(2*p+2*q+1,2);

for ii = 1:length(ID_mat(:,1))
    
    start_node_idx = find(good_nodes==ID_mat(ii,1));
    end_node_idx   = find(good_nodes==ID_mat(ii,2));
    sp             = sort(good_nodes(start_node_idx:end_node_idx));
    ep             = sort(sp,'descend');
    
    keep_searching         = 1;
    
    while keep_searching~=0
        
        for i=1:numel(sp)
            for j=1:numel(ep)
                
                if j==numel(ep)
                    overlap_sequences(cont_overlap_sequences,:) = [sp(i) ep(j)];
                    cont_overlap_sequences                      = cont_overlap_sequences+1;
                    idx = find(good_nodes==ep(j));
                    if idx == numel(good_nodes)
                        keep_searching = 0;
                    else
                        sp = good_nodes(idx+1:end);
                        ep = flipud(sp);
                    end
                    is_solution = 1;
                    break;
                else
                end
                
                
                [is_solution,~]    = DFS(sp(i),ep(j),A_list);
                
                if is_solution
                    overlap_sequences(cont_overlap_sequences,:) = [sp(i) ep(j)];
                    cont_overlap_sequences                      = cont_overlap_sequences+1;
                    idx = find(good_nodes==ep(j));
                    if idx == numel(good_nodes)
                        keep_searching = 0;
                    else
                        sp = good_nodes(idx+1:end);
                        ep = flipud(sp);
                    end
                    break;
                else
                end
            end
            
            if is_solution
                break;
            else
            end
            
        end
        
    end
    
end

% Initial number of valid overlap sequences. Note that this number can be
% lowered later, in case some sequences can be actually merged
n_valid_overlap_sequences = numel(find(overlap_sequences(:,1)~=0));

if n_valid_overlap_sequences~=0
    valid_overlap_sequences_ID = 1;
    overlap_sequences          = overlap_sequences(1:n_valid_overlap_sequences,:);
end

% If there's more than one path solving the problem, check if some paths
% can be merged
if valid_overlap_sequences_ID
    if length(overlap_sequences(:,1)) > 1
        [~,idx]                  = sort(overlap_sequences(:,1));
        overlap_sequences        = overlap_sequences(idx,:);
        merged_overlap_sequences = zeros(size(overlap_sequences));
        cont_os                  = 2;
        
        for i=1:length(overlap_sequences(:,1))-1
            if i == 1
                start_p = overlap_sequences(1,1);
                end_p   = max([overlap_sequences(1,2),overlap_sequences(2,2)]);
                
                [is_solution,~]    = DFS(start_p,end_p,A_list);
                
                if is_solution
                    merged_overlap_sequences(1,:) = [start_p end_p];
                else
                    merged_overlap_sequences(1:2,:) = overlap_sequences(1:2,:);
                end
                
            else
                idx     = find(merged_overlap_sequences(:,1)~=0);
                start_p = min([merged_overlap_sequences(idx(end),1),overlap_sequences(cont_os,1)]);
                end_p   = max([merged_overlap_sequences(idx(end),2),overlap_sequences(cont_os,2)]);
                
                [is_solution,~]    = DFS(start_p,end_p,A_list);
                
                if is_solution
                    merged_overlap_sequences(idx(end),:) = [start_p end_p];
                else
                    merged_overlap_sequences(idx(end)+1,:) = overlap_sequences(cont_os,:);
                end
                
            end
            
            cont_os = cont_os+1;
            
        end
        
        idx                       = find(merged_overlap_sequences(:,1)~=0);
        overlap_sequences         = merged_overlap_sequences(1:idx(end),:);
        n_valid_overlap_sequences = numel(find(overlap_sequences(:,1)~=0));
        
    else
    end
    
else
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 4: if we have identified sequences that show overlapping, %%%
%%% we compute the sub-blocks within the free-space diagram those  %%%
%%% sequences are associated with. Then, we compute if the two     %%%
%%% tracks are characterized by intersections, and we characterize %%%
%%% each intersection with the associated unitary square. We then  %%%
%%% check if these unitary squares are already inside the          %%%
%%% sub-blocks defined above or not. In the second case, we        %%%
%%% determine the extension of the overlapping regions within each %%%
%%% unitary square not contained in any sub-block, and store them  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if valid_overlap_sequences_ID
    unitary_square_matrix = zeros(n_valid_overlap_sequences,4);
    for i=1:n_valid_overlap_sequences
        Type_in  = Node_ID_eps(overlap_sequences(i,1),1);
        R_in     = Node_ID_eps(overlap_sequences(i,1),2);
        C_in     = Node_ID_eps(overlap_sequences(i,1),3);
        Type_fin = Node_ID_eps(overlap_sequences(i,2),1);
        R_fin    = Node_ID_eps(overlap_sequences(i,2),2);
        C_fin    = Node_ID_eps(overlap_sequences(i,2),3);
        
        %%%%%%%%%%%%%%%%%%%%%%
        %%% Starting point %%%
        %%%%%%%%%%%%%%%%%%%%%%
        
        % Corner node
        if Type_in == 1
            
            if R_in == q+1
                r_us_in = q;
                if C_in == p+1
                    c_us_in = p;
                else
                    c_us_in = C_in;
                end
            else
                r_us_in = R_in;
                if C_in == p+1
                    c_us_in = p;
                else
                    c_us_in = C_in;
                end
            end
            % Horizontal node
        elseif Type_in == 2
            
            if R_in == 1
                r_us_in = 1;
                c_us_in = C_in;
            else
                r_us_in = R_in-1;
                c_us_in = C_in;
            end
            % Vertical node
        else
            if C_in == 1
                r_us_in = R_in;
                c_us_in = 1;
            else
                r_us_in = R_in;
                c_us_in = C_in-1;
            end
        end
        
        %%%%%%%%%%%%%%%%
        %%% Endpoint %%%
        %%%%%%%%%%%%%%%%
        
        % Corner node
        if Type_fin == 1
            
            if R_fin == 1
                r_us_fin = 1;
                if C_fin == 1
                    c_us_fin = 1;
                else
                    c_us_fin = C_fin-1;
                end
            else
                r_us_fin = R_fin-1;
                if C_fin == 1
                    c_us_fin = 1;
                else
                    c_us_fin = C_fin-1;
                end
            end
            % Horizontal node
        elseif Type_fin == 2
            
            if R_fin == q+1
                r_us_fin = q;
                c_us_fin = C_fin;
            else
                r_us_fin = R_fin;
                c_us_fin = C_fin;
            end
            % Vertical node
        else
            if C_fin == p+1
                r_us_fin = R_fin;
                c_us_fin = p;
            else
                r_us_fin = R_fin;
                c_us_fin = C_fin;
            end
        end
        unitary_square_matrix(i,:) = [r_us_in r_us_fin c_us_in c_us_fin];
    end
else
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Here we are determining whether the two tracks show some            %%%
%%% intersections. To do so, we use triangular geometry as suggested in %%%
%%% http://enrico.spinielli.net/understanding-great-circle-arcs_57/     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int_matrix = zeros(q,p);

for i=1:q
    a0     = latlon2cart(Q_red(i,:),1);
    a1     = latlon2cart(Q_red(i+1,:),1);
    p_vers = cross(a0,a1)/norm(cross(a0,a1));
    for j=1:p
        b0     = latlon2cart(P_red(j,:),1);
        b1     = latlon2cart(P_red(j+1,:),1);
        q_vers = cross(b0,b1)/norm(cross(b0,b1));
        t_vers = cross(p_vers,q_vers)/norm(cross(p_vers,q_vers));
        s1     = dot(cross(a0,p_vers),t_vers);
        s2     = dot(cross(a1,p_vers),t_vers);
        s3     = dot(cross(b0,q_vers),t_vers);
        s4     = dot(cross(b1,q_vers),t_vers);
        if -s1>=0 && s2>=0 && -s3>=0 && s4>=0
            int_matrix(i,j) = 1;
        elseif -s1<=0 && s2<=0 && -s3<=0 && s4<=0
            int_matrix(i,j) = 1;
        else
        end
    end
end

[r,c] = find(int_matrix);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% If there's at least one intersection between routes, we check %%%
%%% if the intersections are already included in the sequences we %%%
%%% have computed or not                                          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If there's no intersection at all, just skip this check
if ~isempty(r)
    already_accounted_for = zeros(numel(r),1);
    % If we have identified overlapping regions, scan every unitary square
    % characterized by an intersection and check it it's already included
    % into one of the overlapping regions
    if valid_overlap_sequences_ID
        for i=1:numel(r)
            for j=1:length(unitary_square_matrix(:,1))
                % If the unitary square is inside one of the overlapping
                % regions, do not scan the overlapping regions left and
                % go to the next intersection
                if r(i)>=unitary_square_matrix(j,1) && r(i)<=unitary_square_matrix(j,2) && ...
                        c(i)>=unitary_square_matrix(j,3) && c(i)<=unitary_square_matrix(j,4)
                    already_accounted_for(i) = 1;
                    break;
                else
                end
                
            end
        end
        % If there are intersections, but no overlapping regions, we need to
        % consider all the intersections, thus we do not modify the
        % already_accounted_for vector
    else
    end
    
    idx = find(already_accounted_for==0);
    
    if ~isempty(idx)
        
        valid_overlap_sequences_int_ID = 1;
        
        Overlap_points_int = zeros(numel(idx),16);
        N_f                = 20;
        
        for i=1:numel(idx)
            % row of unitary square considered
            us_r = r(idx(i));
            % column of unitary square considered
            us_c = c(idx(i));
            
            Q_tail           = Q_red(us_r,:);
            Q_head           = Q_red(us_r+1,:);
            P_tail           = P_red(us_c,:);
            P_head           = P_red(us_c+1,:);
            Q_tail_cart      = latlon2cart(Q_tail,1);
            Q_head_cart      = latlon2cart(Q_head,1);
            P_tail_cart      = latlon2cart(P_tail,1);
            P_head_cart      = latlon2cart(P_head,1);
            DQ               = Q_length(us_r);
            DP               = P_length(us_c);
            f_vec_P          = linspace(0,1,N_f);
            f_vec_Q          = linspace(0,1,N_f);
            dist_mat         = zeros(numel(f_vec_Q),numel(f_vec_P));
            dist_tail        = zeros(size(dist_mat));
            dist_head        = zeros(size(dist_mat));
            for ii=1:numel(f_vec_Q)
                Q_cart   = Q_tail_cart*sin((1-f_vec_Q(ii))*DQ)/sin(DQ)+Q_head_cart*sin(f_vec_Q(ii)*DQ)/sin(DQ);
                Q_ll     = cart2latlon(Q_cart);
                for jj=1:numel(f_vec_P)
                    P_cart           = P_tail_cart*sin((1-f_vec_P(jj))*DP)/sin(DP)+P_head_cart*sin(f_vec_P(jj)*DP)/sin(DP);
                    P_ll             = cart2latlon(P_cart);
                    dist_mat(ii,jj)  = Haversine(P_ll,Q_ll,1);
                    dist_tail(ii,jj) = Haversine(P_tail,P_ll,1)+Haversine(Q_tail,Q_ll,1);
                    dist_head(ii,jj) = Haversine(P_ll,P_head,1)+Haversine(Q_ll,Q_head,1);
                end
            end
            [r,c]    = find(dist_mat<=alpha*FD);
            dist_tail_vec = zeros(numel(r),1);
            dist_head_vec = zeros(numel(r),1);
            for ii=1:numel(r)
                dist_tail_vec(ii) = dist_tail(r(ii),c(ii));
                dist_head_vec(ii) = dist_head(r(ii),c(ii));
            end
            [~,idx_sp]       = max(dist_head_vec);
            [~,idx_ep]       = max(dist_tail_vec);
            sp_best_f_Q      = f_vec_Q(r(idx_sp));
            sp_best_f_P      = f_vec_P(c(idx_sp));
            ep_best_f_Q      = f_vec_Q(r(idx_ep));
            ep_best_f_P      = f_vec_P(c(idx_ep));
            P_cart_sp        = P_tail_cart*sin((1-sp_best_f_P)*DP)/sin(DP)+P_head_cart*sin(sp_best_f_P*DP)/sin(DP);
            Q_cart_sp        = Q_tail_cart*sin((1-sp_best_f_Q)*DQ)/sin(DQ)+Q_head_cart*sin(sp_best_f_Q*DQ)/sin(DQ);
            P_cart_ep        = P_tail_cart*sin((1-ep_best_f_P)*DP)/sin(DP)+P_head_cart*sin(ep_best_f_P*DP)/sin(DP);
            Q_cart_ep        = Q_tail_cart*sin((1-ep_best_f_Q)*DQ)/sin(DQ)+Q_head_cart*sin(ep_best_f_Q*DQ)/sin(DQ);
            sp_P             = cart2latlon(P_cart_sp);
            sp_Q             = cart2latlon(Q_cart_sp);
            ep_P             = cart2latlon(P_cart_ep);
            ep_Q             = cart2latlon(Q_cart_ep);
            sp_unit_square_Q = us_r;
            sp_unit_square_P = us_c;
            ep_unit_square_Q = us_r;
            ep_unit_square_P = us_c;
            
            Overlap_points_int(i,:) = horzcat(sp_unit_square_P,sp_best_f_P,sp_P,sp_unit_square_Q,sp_best_f_Q,sp_Q,...
                ep_unit_square_P,ep_best_f_P,ep_P,ep_unit_square_Q,ep_best_f_Q,ep_Q);
            
        end
        % All intersection have already been accounted for
    else
    end
    
    
else
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 5: is some sequences have been identified in STEP 3, translate %%%
%%% each starting point/endpoint into a couple of lat/lon pairs, one    %%%
%%% on P and one on Q. For each couple of lat/lon pairs, we store the   %%%
%%% row and column of the unitary square they refer to, and the f value %%%
%%% (0<=f<=1) that defines where the lat/lon pair is located along the  %%%
%%% sides of the unitary square                                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if valid_overlap_sequences_ID
    
    Merged_Sequences = zeros(length(overlap_sequences(:,1)),8);
    for i=1:length(overlap_sequences(:,1))
        Merged_Sequences(i,:) = [Node_ID_eps(overlap_sequences(i,1),:) Node_ID_eps(overlap_sequences(i,2),:)];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% For each sequence, translate the starting and endpoint of the %%%
    %%% sequence into the corresponding points on curves P and Q, and %%%
    %%% store the values in a new matrix                              %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [Overlap_points,~,~] = compute_sp_ep(Merged_Sequences,P_red,Q_red,alpha,FD,Bmin,Bmax,Lmin,Lmax);
    
else
end



is_there_overlapping = 1;

if valid_overlap_sequences_ID && valid_overlap_sequences_int_ID
    disp('Sequences both from free space and from intersections')
    Overlap_points = vertcat(Overlap_points,Overlap_points_int);
elseif valid_overlap_sequences_ID
    disp('Sequences from free space')
elseif valid_overlap_sequences_int_ID
    disp('Sequences from intersections')
    Overlap_points = Overlap_points_int;
else
    is_there_overlapping = 0;
    disp('No sequences')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 6: for each overlapping region, defined as initial/final    %%%
%%% lat/lon and initial/final unitary square and f values, use these %%%
%%% information to compute the overall length of such overlapping    %%%
%%% region                                                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if is_there_overlapping
    overlapping_region_length = zeros(length(Overlap_points(:,1)),4);
    
    overall_length_P = sum(P_length);
    overall_length_Q = sum(Q_length);
    
    for i = 1:length(Overlap_points(:,1))
        in_us_P  = Overlap_points(i,1);
        fin_us_P = Overlap_points(i,9);
        in_us_Q  = Overlap_points(i,5);
        fin_us_Q = Overlap_points(i,13);
        %%%%%%%%%%%%%%%%%%%%%%
        %%% Flight track P %%%
        %%%%%%%%%%%%%%%%%%%%%%
        
        % Current overlapping region lies within a single unitary square
        if fin_us_P-in_us_P == 0
            overlapping_region_length(i,1) = Haversine(Overlap_points(i,3:4),Overlap_points(i,11:12),1);
            overlapping_region_length(i,2) = overlapping_region_length(i,1)/overall_length_P*100;
        % Current overlapping region lies within two consecutive 
        % unitary squares    
        elseif fin_us_P-in_us_P == 1
            % Compute first contribution
            P_tail                         = P_red(in_us_P,:);
            P_head                         = P_red(in_us_P+1,:);
            P_tail_cart                    = latlon2cart(P_tail,1);
            P_head_cart                    = latlon2cart(P_head,1);
            DP                             = P_length(in_us_P);
            f_P                            = Overlap_points(i,2);
            P_cart                         = P_tail_cart*sin((1-f_P)*DP)/sin(DP)+P_head_cart*sin(f_P*DP)/sin(DP);
            P                              = cart2latlon(P_cart);
            P_l1                           = Haversine(P,P_head,1);
           % Compute second contribution
            P_tail                         = P_red(in_us_P+1,:);
            P_head                         = P_red(in_us_P+2,:);
            P_tail_cart                    = latlon2cart(P_tail,1);
            P_head_cart                    = latlon2cart(P_head,1);
            DP                             = P_length(in_us_P+1);
            f_P                            = Overlap_points(i,10);
            P_cart                         = P_tail_cart*sin((1-f_P)*DP)/sin(DP)+P_head_cart*sin(f_P*DP)/sin(DP);
            P                              = cart2latlon(P_cart);
            P_l2                           = Haversine(P_tail,P,1);
            
            overlapping_region_length(i,1) = P_l1+P_l2;
            overlapping_region_length(i,2) = overlapping_region_length(i,1)/overall_length_P*100;
        % Current overlapping region lies within at least three consecutive 
        % unitary squares
        else
            % Compute first contribution
            P_tail                         = P_red(in_us_P,:);
            P_head                         = P_red(in_us_P+1,:);
            P_tail_cart                    = latlon2cart(P_tail,1);
            P_head_cart                    = latlon2cart(P_head,1);
            DP                             = P_length(in_us_P);
            f_P                            = Overlap_points(i,2);
            P_cart                         = P_tail_cart*sin((1-f_P)*DP)/sin(DP)+P_head_cart*sin(f_P*DP)/sin(DP);
            P                              = cart2latlon(P_cart);
            P_l1                           = Haversine(P,P_head,1);
            % Compute second contribution
            P_l2                           = sum(P_length(in_us_P+1:fin_us_P-1));
            % Compute third contribution
            P_tail                         = P_red(fin_us_P,:);
            P_head                         = P_red(fin_us_P+1,:);
            P_tail_cart                    = latlon2cart(P_tail,1);
            P_head_cart                    = latlon2cart(P_head,1);
            DP                             = P_length(in_us_P);
            f_P                            = Overlap_points(i,10);
            P_cart                         = P_tail_cart*sin((1-f_P)*DP)/sin(DP)+P_head_cart*sin(f_P*DP)/sin(DP);
            P                              = cart2latlon(P_cart);
            P_l3                           = Haversine(P_tail,P,1);
            
            overlapping_region_length(i,1) = P_l1+P_l2+P_l3;
            overlapping_region_length(i,2) = overlapping_region_length(i,1)/overall_length_P*100;
        end
        %%%%%%%%%%%%%%%%%%%%%%
        %%% Flight track Q %%%
        %%%%%%%%%%%%%%%%%%%%%%
        
        % Current overlapping region lies within a single unitary square
        if fin_us_Q-in_us_Q == 0
            overlapping_region_length(i,3) = Haversine(Overlap_points(i,7:8),Overlap_points(i,15:16),1);
            overlapping_region_length(i,4) = overlapping_region_length(i,3)/overall_length_Q*100;
        % Current overlapping region lies within two consecutive 
        % unitary squares
        elseif fin_us_Q-in_us_Q == 1
            % Compute first contribution
            Q_tail                         = Q_red(in_us_Q,:);
            Q_head                         = Q_red(in_us_Q+1,:);
            Q_tail_cart                    = latlon2cart(Q_tail,1);
            Q_head_cart                    = latlon2cart(Q_head,1);
            DQ                             = Q_length(in_us_Q);
            f_Q                            = Overlap_points(i,6);
            Q_cart                         = Q_tail_cart*sin((1-f_Q)*DQ)/sin(DQ)+Q_head_cart*sin(f_Q*DQ)/sin(DQ);
            Q                              = cart2latlon(Q_cart);
            Q_l1                           = Haversine(Q,Q_head,1);
           % Compute second contribution
            Q_tail                         = Q_red(in_us_P+1,:);
            Q_head                         = Q_red(in_us_P+2,:);
            Q_tail_cart                    = latlon2cart(Q_tail,1);
            Q_head_cart                    = latlon2cart(Q_head,1);
            DQ                             = Q_length(in_us_Q+1);
            f_Q                            = Overlap_points(i,14);
            Q_cart                         = Q_tail_cart*sin((1-Q_P)*DQ)/sin(DQ)+Q_head_cart*sin(f_Q*DQ)/sin(DQ);
            Q                              = cart2latlon(Q_cart);
            Q_l2                           = Haversine(Q_tail,Q,1);
            
            overlapping_region_length(i,3) = Q_l1+Q_l2;
            overlapping_region_length(i,4) = overlapping_region_length(i,3)/overall_length_Q*100;
        % Current overlapping region lies within at least three consecutive 
        % unitary squares
        else
            % Compute first contribution
            Q_tail                         = Q_red(in_us_Q,:);
            Q_head                         = Q_red(in_us_Q+1,:);
            Q_tail_cart                    = latlon2cart(Q_tail,1);
            Q_head_cart                    = latlon2cart(Q_head,1);
            DQ                             = Q_length(in_us_Q);
            f_Q                            = Overlap_points(i,6);
            Q_cart                         = Q_tail_cart*sin((1-f_Q)*DQ)/sin(DQ)+Q_head_cart*sin(f_Q*DQ)/sin(DQ);
            Q                              = cart2latlon(Q_cart);
            Q_l1                           = Haversine(Q,Q_head,1);
            % Compute second contribution
            Q_l2                           = sum(Q_length(in_us_Q+1:fin_us_Q-1));
            % Compute third contribution
            Q_tail                         = Q_red(fin_us_Q,:);
            Q_head                         = Q_red(fin_us_Q+1,:);
            Q_tail_cart                    = latlon2cart(Q_tail,1);
            Q_head_cart                    = latlon2cart(Q_head,1);
            DQ                             = Q_length(fin_us_Q);
            f_Q                            = Overlap_points(i,14);
            Q_cart                         = Q_tail_cart*sin((1-f_Q)*DQ)/sin(DQ)+Q_head_cart*sin(f_Q*DQ)/sin(DQ);
            Q                              = cart2latlon(Q_cart);
            Q_l3                           = Haversine(Q_tail,Q,1);
            
            overlapping_region_length(i,3) = Q_l1+Q_l2+Q_l3;
            overlapping_region_length(i,4) = overlapping_region_length(i,3)/overall_length_Q*100;
        end
    end
else
end

return
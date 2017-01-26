% Written by Alessandro Bombelli, May 29th 2016
% Comparing performance of undersampling routine for frecher Distance 
% Computation 

clc
clear all
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Re                           = 1;
eps                          = 0.001;
tolerance                    = 1;


P = [0.1 -10;0.1 -8;0.15 -6;0.25 -5;0.3 -4;0.25 -3;-0.05 -2;-0.1 -1;-0.15 0;-0.05 2;-0.03 4;0.1 5;0.2 6;0.15 7;0.1 9; 0.05 10];
Q = horzcat(zeros(11,1),(-10:2:10)'); 
P_red = P;

A = load(strcat(pwd,'/TRX_PROCESSED/2014070100_2014070124_1of3.mat'));

names = fieldnames(A.ACID_info);

for i=1:length(names)
    flight      = A.ACID_info.(names{i});
    origin      = flight.origin;
    destination = flight.destination;
    if strcmp(origin,'KJFK')==1
        if strcmp(destination,'KLAX')==1
            disp(['Flight ',num2str(i),' is good'])
        else
        end
    else
    end
end

 P = A.ACID_info.(names{146}).lat_lon;
 Q = [40.6398 -73.7789;40.2 -74.5;40 -76.5;38.33 -81.8;36.1 -86.65;35 -89.95;...
      34.66 -92.16;35.33 -97.6;35.25 -101.6;34.6 -112.5;34.05 -115.7;...
      33.95 -117.5;33.9425 -118.4081];
headingP = traj_heading(P);

max_delta_heading   = 4;
delta_heading1      = abs(headingP(2:end)-headingP(1:end-1));
dummy1              = delta_heading1<=max_delta_heading;
idx1                = find(dummy1==1)+1;
P_red               = P;
P_red(idx1,:)       = [];


[cd_all,FD,sol,no_sol,Lmin,Lmax,Bmin,Bmax,A,sequence,Node_ID] = FD_computation_v03(P_red,Q,1);


figure()
hold on
h1 = plot(P_red(:,2),P_red(:,1),'Color','b','Linewidth',2,'Marker','*','Markersize',8);
h2 = plot(Q(:,2),Q(:,1),'Color','k','Linewidth',2,'Marker','*','Markersize',8);
xlabel('Longitude [deg]','Fontname','Avantgarde','Fontsize',14)
ylabel('Latitude [deg]','Fontname','Avantgarde','Fontsize',14)
l1 = legend([h1,h2],'P polygonal curve','Q polygonal curve');
set(l1,'Fontname','Avantgarde','Fontsize',14)
grid on


tic

alpha = 1/4;
[Lmin_eps,Lmax_eps,Bmin_eps,Bmax_eps,A_list_eps,Node_ID_eps] = free_space_shape_v02(P_red,Q,Re,alpha*FD);
plot_FS(Bmin_eps,Bmax_eps,Lmin_eps,Lmax_eps)

[P_length,Q_length,dist_vector,NODES] = overlap_index_v01(P_red,Q,sequence,Node_ID,Lmin,Lmax,Bmin,Bmax);
idx                                   = find(dist_vector<=alpha*FD);
seq2                                  = NODES(idx(2:end),1)-NODES(idx(1:end-1),1);
check_vector                          = 1;
cont                                  = 1;
BLOCKS                                = cell(0,0);
cont_block                            = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while check_vector ~= 0
    
    if cont > numel(idx)
        check_vector = 0;
        break
    else
    end
    
    if cont == numel(idx)
        this_block           = idx(cont);
        BLOCKS{cont_block,1} = this_block;
        cont                 = cont+1;
        cont_block           = cont_block+1;
        check_vector         = 0;
    else
        % Sequence is formed by a single index
        if idx(cont+1)-idx(cont)~=1
            this_block           = idx(cont);
            BLOCKS{cont_block,1} = this_block;
            cont                 = cont+1;
            cont_block           = cont_block+1;
            % Sequence is formed by more than one index
        else
            
            this_block = idx(cont);
            while idx(cont+1)-idx(cont)==1
                this_block = horzcat(this_block,idx(cont+1));
                cont       = cont+1;
                if cont == numel(idx)
                    check_vector = 0;
                    break
                else
                end
            end
            BLOCKS{cont_block,1} = this_block;
            cont                 = cont+1;
            cont_block           = cont_block+1;
        end
        
    end
end

time1 = toc;

tic

if ~isempty(BLOCKS)
    mu_R      = 3;
    mu_C      = 3;
    p         = length(P_red(:,1))-1;
    q         = length(Q(:,1))-1;
    Sequences = [];
    cont      = 1;
    
    for i=1:length(BLOCKS)
        this_sequence      = BLOCKS{i,1};
        % Initial startpoint and endpoint used as inputs for the DFS
        
        if numel(this_sequence)>1
            initial_startpoint      = NODES(this_sequence(1),2);
            initial_startpoint_type = NODES(this_sequence(1),3);
            initial_startpoint_R    = NODES(this_sequence(1),4);
            initial_startpoint_C    = NODES(this_sequence(1),5);
            initial_endpoint        = NODES(this_sequence(2),2);
            initial_endpoint_type   = NODES(this_sequence(2),3);
            initial_endpoint_R      = NODES(this_sequence(2),4);
            initial_endpoint_C      = NODES(this_sequence(2),5);
         else
            initial_startpoint      = NODES(this_sequence(1),2);
            initial_startpoint_type = NODES(this_sequence(1),3);
            initial_startpoint_R    = NODES(this_sequence(1),4);
            initial_startpoint_C    = NODES(this_sequence(1),5);
            initial_endpoint        = initial_startpoint;
            initial_endpoint_type   = initial_startpoint_type;
            initial_endpoint_R      = initial_startpoint_R;
            initial_endpoint_C      = initial_startpoint_C;
        end
        
        Rin      = max([1,initial_startpoint_R-mu_R]);
        Cin      = max([1,initial_startpoint_C-mu_C]);
        Rfin     = min([initial_endpoint_R+mu_R,q+1]);
        Cfin     = min([initial_endpoint_C+mu_C,p+1]);
        Rvec_in  = (Rin:initial_startpoint_R)';
        Cvec_in  = (Cin:initial_startpoint_C)';
        Rvec_fin = (initial_endpoint_R:Rfin)';
        Cvec_fin = (initial_endpoint_C:Cfin)';
        
        % Creating set of starting points
        if initial_startpoint_type==1
            extra_starting = [];
            for j=1:numel(Rvec_in)
                if Rvec_in(j)~=initial_startpoint_R
                    dummy1 = (Rvec_in(j)-1)*(3*p+2)+(2*Cin-1:2*initial_startpoint_C-1)';
                    dummy2 = (Rvec_in(j)-1)*(3*p+2)+2*p+1+Cvec_in;
                else
                    dummy1 = (Rvec_in(j)-1)*(3*p+2)+(2*Cin-1:2*initial_startpoint_C-1)';
                    dummy2 = [];
                end
                extra_starting = [extra_starting;dummy1;dummy2];
            end
            
        elseif initial_startpoint_type==2
            extra_starting = [];
            for j=1:numel(Rvec_in)
                if Rvec_in(j)~=initial_startpoint_R
                    dummy1 = (Rvec_in(j)-1)*(3*p+2)+(2*Cin-1:2*initial_startpoint_C)';
                    dummy2 = (Rvec_in(j)-1)*(3*p+2)+2*p+1+Cvec_in;
                else
                    dummy1 = (Rvec_in(j)-1)*(3*p+2)+(2*Cin-1:2*initial_startpoint_C)';
                    dummy2 = [];
                end
                extra_starting = [extra_starting;dummy1;dummy2];
            end
        else
            extra_starting = [];
            for j=1:numel(Rvec_in)
                dummy1 = (Rvec_in(j)-1)*(3*p+2)+(2*Cin-1:2*initial_startpoint_C-1)';
                dummy2 = (Rvec_in(j)-1)*(3*p+2)+2*p+1+Cvec_in;
                extra_starting = [extra_starting;dummy1;dummy2];
            end
        end
        
        % Creating set of end points
        if initial_endpoint_type==1
            extra_ending = [];
            for j=1:numel(Rvec_fin)
                if Rvec_fin(j)~=Rfin
                    dummy1 = (Rvec_fin(j)-1)*(3*p+2)+(2*initial_startpoint_C-1:2*Cfin-1)';
                    dummy2 = (Rvec_fin(j)-1)*(3*p+2)+2*p+1+(initial_startpoint_C:Cfin)';
                else
                    dummy1 = (Rvec_fin(j)-1)*(3*p+2)+(2*initial_startpoint_C-1:2*Cfin-1)';
                    dummy2 = [];
                end
                extra_ending = [extra_ending;dummy1;dummy2];
            end
        elseif initial_endpoint_type==2
            extra_ending = [];
            for j=1:numel(Rvec_fin)
                if Rvec_fin(j)~=Rfin
                    dummy1 = (Rvec_fin(j)-1)*(3*p+2)+(2*initial_startpoint_C:2*Cfin+1)';
                    dummy2 = (Rvec_fin(j)-1)*(3*p+2)+2*p+2+(initial_startpoint_C:Cfin)';
                else
                    dummy1 = (Rvec_fin(j)-1)*(3*p+2)+(2*initial_startpoint_C:2*Cfin+1)';
                    dummy2 = [];
                end
                extra_ending = [extra_ending;dummy1;dummy2];
            end
        else
            extra_ending = [];
            for j=1:numel(Rvec_fin)
                dummy1 = (Rvec_fin(j)-1)*(3*p+2)+2*p+1+Cvec_fin;
                dummy2 = Rvec_fin(j)*(3*p+2)+(2*initial_startpoint_C-1:2*Cfin-1)';
                extra_ending = [extra_ending;dummy1;dummy2];
            end
        end
        
        
        for ii=1:numel(extra_starting)
            sp = extra_starting(ii);
            for jj=1:numel(extra_ending)
                ep = extra_ending(numel(extra_ending)-jj+1);
                [is_solution,~]    = DFS(sp,ep,A_list_eps);
                
                if is_solution
                    Sequences(cont,:) = horzcat(Node_ID_eps(sp,:),Node_ID_eps(ep,:));
                    cont              = cont+1;
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
else
end

time2 = toc;

tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Checking if some sequences can be merged %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Trying to merge sequences makes sense only if there are at least 2
% sequences in the original array
if length(Sequences(:,1))>=2
    
    Merged_Sequences = [];
    
    for i=1:length(Sequences(:,1))-1
        
        % The first 2 sequences are compared. Is there's a path from the
        % starting point of the first to the endpoint of the second, merge
        % them, otherwise store them as separate sequences
        if i==1
            last_sp = Sequences(1,4);
            last_ep = Sequences(1,8);
            new_sp  = Sequences(2,4);
            new_ep  = Sequences(2,8);
            
            [is_solution,~]    = DFS(last_sp,new_ep,A_list_eps);
            % In this case, there's a path. Sequences are merged
            if is_solution
                Merged_Sequences = [Merged_Sequences;Sequences(1,1:4) Sequences(2,5:8)];
                % Otherwise, store them separately
            else
                Merged_Sequences = [Merged_Sequences;Sequences(1,:);Sequences(2,:)];
            end
            
            % the next new sequence which will be checked is the third one
            cont = 3;
            
            % After i=1, we check the last sequence stored in the new array (which
            % could be a combination of sequences previously computed) with the
            % next old sequence that has not been checked yet
        else
            last_sp = Merged_Sequences(end,4);
            last_ep = Merged_Sequences(end,8);
            new_sp  = Sequences(cont,4);
            new_ep  = Sequences(cont,8);
            
            [is_solution,~]    = DFS(last_sp,new_ep,A_list_eps);
            % In this case, there's a path. Sequences are merged
            if is_solution
                Merged_Sequences = [Merged_Sequences(1:end-1,:);Merged_Sequences(end,1:4) Sequences(cont,5:8)];
                % Otherwise, store them separately
            else
                Merged_Sequences = [Merged_Sequences;Sequences(cont,:)];
            end
            
            % Increase counter
            cont = cont+1;
            
        end
        
        
        
    end
    
else
    Merged_Sequences = Sequences;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define for each block of Merged_Sequences the initial/final row  %%%
%%% and column of the subset of unitary squares that contain that    %%%
%%% specific sequence                                                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Merged_Sequences_rc = zeros(length(Merged_Sequences(:,1)),4);
% First column : row of the starting unitary square
% Second column: row of the ending unitary square
% Third column : column of the starting unitary square
% Fourth column: column of the ending unitary square

for i=1:length(Merged_Sequences(:,1))
    % Get information of initial/final nodes
    sp_type = Merged_Sequences(i,1);
    sp_R    = Merged_Sequences(i,2);
    sp_C    = Merged_Sequences(i,3);
    ep_type = Merged_Sequences(i,5);
    ep_R    = Merged_Sequences(i,6);
    ep_C    = Merged_Sequences(i,7);
    
    if sp_type == 1
        if sp_R == 1 && sp_C == 1
            Merged_Sequences_rc(i,1) = sp_R;
            Merged_Sequences_rc(i,3) = sp_C;
        elseif sp_R == 1
            Merged_Sequences_rc(i,1) = sp_R;
            if sp_C == p+1
                Merged_Sequences_rc(i,3) = p;
            else
                Merged_Sequences_rc(i,3) = sp_C;
            end
        elseif sp_C == 1
            Merged_Sequences_rc(i,3) = sp_C;
            if sp_R == q+1
                Merged_Sequences_rc(i,1) = q;
            else
                Merged_Sequences_rc(i,1) = sp_R;
            end
        else
            Merged_Sequences_rc(i,1) = sp_R;
            Merged_Sequences_rc(i,3) = sp_C;
        end
    elseif sp_type == 2
        if sp_R == q+1
            Merged_Sequences_rc(i,1) = q;
            Merged_Sequences_rc(i,3) = sp_C;
        else
            Merged_Sequences_rc(i,1) = sp_R-1;
            Merged_Sequences_rc(i,3) = sp_C;
        end
    else
        if sp_C == p+1
            Merged_Sequences_rc(i,1) = sp_R;
            Merged_Sequences_rc(i,3) = p;
        else
            Merged_Sequences_rc(i,1) = sp_R;
            Merged_Sequences_rc(i,3) = sp_C-1;
        end
    end
    
    if ep_type == 1
        if ep_R == q+1 && ep_C == p+1
            Merged_Sequences_rc(i,2) = q;
            Merged_Sequences_rc(i,4) = p;
        elseif ep_R == 1
            Merged_Sequences_rc(i,2) = ep_R;
            if ep_C == p+1
                Merged_Sequences_rc(i,4) = p;
            else
                Merged_Sequences_rc(i,4) = ep_C;
            end
        elseif ep_C == 1
            Merged_Sequences_rc(i,4) = sp_C;
            if ep_R == q+1
                Merged_Sequences_rc(i,2) = q;
            else
                Merged_Sequences_rc(i,2) = ep_R;
            end
        else
            Merged_Sequences_rc(i,2) = ep_R;
            Merged_Sequences_rc(i,4) = ep_C;
        end
    elseif ep_type == 2
        if ep_R == q+1
            Merged_Sequences_rc(i,2) = q;
            Merged_Sequences_rc(i,4) = sp_C;
        else
            Merged_Sequences_rc(i,2) = ep_R;
            Merged_Sequences_rc(i,4) = ep_C;
        end
    else
        if ep_C == p+1
            Merged_Sequences_rc(i,2) = ep_R;
            Merged_Sequences_rc(i,4) = p;
        else
            Merged_Sequences_rc(i,2) = ep_R;
            Merged_Sequences_rc(i,4) = ep_C;
        end
    end
    
end

time3 = toc;

tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Check if there are intersections between great circles that       %%%
%%% result in feasible regions within some unitary squares, without   %%%
%%% touching the sides of the unitary square                          %%%
%%% Ref:                                                              %%%
%%% http://enrico.spinielli.net/understanding-great-circle-arcs_57/   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int_matrix = zeros(q,p);

for i=1:q
    a0     = latlon2cart(Q(i,:),1);
    a1     = latlon2cart(Q(i+1,:),1);
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

time4 = toc;

tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% For each sequence, translate the starting and endpoint of the %%%
%%% sequence into the corresponding points on curves P and Q, and %%%
%%% store the values in a new matrix                              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Overlap_points  = zeros(length(Merged_Sequences(:,1)),8);
Overlap_points2 = zeros(length(Merged_Sequences(:,1)),16);
N_f             = 20;

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
        disp(['Starting node is corner. Row is ',num2str(r_sp),', column is ',num2str(c_sp)])
        if r_sp == 1 && c_sp == 1
            sp_P          = P_red(c_sp,:);
            sp_Q          = Q(r_sp,:);
            best_f_Q      = 0;
            best_f_P      = 0;
            unit_square_Q = r_sp;
            unit_square_P = c_sp;
        else
            f_vec_P     = linspace(0,1,N_f);
            f_vec_Q     = linspace(0,1,N_f);
            P_tail      = P_red(c_sp-1,:);
            P_head      = P_red(c_sp,:);
            Q_tail      = Q(r_sp-1,:);
            Q_head      = Q(r_sp,:);
            P_tail_cart = latlon2cart(P_tail,1);
            P_head_cart = latlon2cart(P_head,1);
            Q_tail_cart = latlon2cart(Q_tail,1);
            Q_head_cart = latlon2cart(Q_head,1);
            DP          = P_length(c_sp-1);
            DQ          = Q_length(r_sp-1);
            dist_mat    = zeros(numel(f_vec_Q),numel(f_vec_P));
            dist_tail       = zeros(size(dist_mat));
            for ii=1:numel(f_vec_Q)
                Q_cart   = Q_tail_cart*sin((1-f_vec_Q(ii))*DQ)/sin(DQ)+Q_head_cart*sin(f_vec_Q(ii)*DQ)/sin(DQ);
                Q_ll     = cart2latlon(Q_cart);
                for jj=1:numel(f_vec_P)
                    P_cart   = P_tail_cart*sin((1-f_vec_P(jj))*DP)/sin(DP)+P_head_cart*sin(f_vec_P(jj)*DP)/sin(DP);
                    P_ll     = cart2latlon(P_cart);
                    dist_mat(ii,jj) = Haversine(P_ll,Q_ll,1);
                    dist_tail(ii,jj)    = Haversine(P_ll,P_head,1)+Haversine(Q_ll,Q_head,1);
                end
            end
            [r,c]    = find(dist_mat<=alpha*FD);
            dist_tail_vec = zeros(numel(r));
            for ii=1:numel(r)
                dist_tail_vec(ii) = dist_tail(r(ii),c(ii));
            end
            [~,idx]   = max(dist_tail_vec);
            best_f_Q      = f_vec_Q(r(idx));
            best_f_P      = f_vec_P(c(idx));
            unit_square_Q = r_sp-1;
            unit_square_P = c_sp-1;
            P_cart_sp = P_tail_cart*sin((1-best_f_P)*DP)/sin(DP)+P_head_cart*sin(best_f_P*DP)/sin(DP);
            Q_cart_sp = Q_tail_cart*sin((1-best_f_Q)*DQ)/sin(DQ)+Q_head_cart*sin(best_f_Q*DQ)/sin(DQ);
            sp_P      = cart2latlon(P_cart_sp);
            sp_Q      = cart2latlon(Q_cart_sp);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%
        %%% Horizontal node %%%
        %%%%%%%%%%%%%%%%%%%%%%%
    elseif Type_sp == 2
        disp(['Starting node is horizontal. Row is ',num2str(r_sp),', column is ',num2str(c_sp)])
        if r_sp == 1
            P_tail      = P_red(c_sp,:);
            P_head      = P_red(c_sp+1,:);
            P_tail_cart = latlon2cart(P_tail,1);
            P_head_cart = latlon2cart(P_head,1);
            best_f_Q      = 0;
            best_f_P      = Bmin(r_sp,c_sp);
            unit_square_Q = r_sp;
            unit_square_P = c_sp;
            DP          = P_length(c_sp);
            P_cart_sp   = P_tail_cart*sin((1-best_f_P)*DP)/sin(DP)+P_head_cart*sin(best_f_P*DP)/sin(DP);
            sp_P        = cart2latlon(P_cart_sp);
            sp_Q        = Q(r_sp,:);
        else
            f_vec_P     = linspace(0,Bmax(r_sp,c_sp),N_f);
            f_vec_Q     = linspace(0,1,N_f);
            P_tail      = P_red(c_sp,:);
            P_head      = P_red(c_sp+1,:);
            Q_tail      = Q(r_sp-1,:);
            Q_head      = Q(r_sp,:);
            P_tail_cart = latlon2cart(P_tail,1);
            P_head_cart = latlon2cart(P_head,1);
            Q_tail_cart = latlon2cart(Q_tail,1);
            Q_head_cart = latlon2cart(Q_head,1);
            DP          = P_length(c_sp);
            DQ          = Q_length(r_sp-1);
            dist_mat    = zeros(numel(f_vec_Q),numel(f_vec_P));
            dist_tail       = zeros(size(dist_mat));
            for ii=1:numel(f_vec_Q)
                Q_cart   = Q_tail_cart*sin((1-f_vec_Q(ii))*DQ)/sin(DQ)+Q_head_cart*sin(f_vec_Q(ii)*DQ)/sin(DQ);
                for jj=1:numel(f_vec_P)
                    P_cart   = P_tail_cart*sin((1-f_vec_P(jj))*DP)/sin(DP)+P_head_cart*sin(f_vec_P(jj)*DP)/sin(DP);
                    P_ll     = cart2latlon(P_cart);
                    Q_ll     = cart2latlon(Q_cart);
                    dist_mat(ii,jj) = Haversine(P_ll,Q_ll,1);
                    dist_tail(ii,jj)    = Haversine(P_ll,P_head,1)+Haversine(Q_ll,Q_head,1);
                end
            end
            [r,c]    = find(dist_mat<=alpha*FD);
            dist_tail_vec = zeros(numel(r));
            for ii=1:numel(r)
                dist_tail_vec(ii) = dist_tail(r(ii),c(ii));
            end
            [~,idx] = max(dist_tail_vec);
            best_f_Q      = f_vec_Q(r(idx));
            best_f_P      = f_vec_P(c(idx));
            unit_square_Q = r_sp-1;
            unit_square_P = c_sp;
            P_cart_sp = P_tail_cart*sin((1-best_f_P)*DP)/sin(DP)+P_head_cart*sin(best_f_P*DP)/sin(DP);
            Q_cart_sp = Q_tail_cart*sin((1-best_f_Q)*DQ)/sin(DQ)+Q_head_cart*sin(best_f_Q*DQ)/sin(DQ);
            sp_P      = cart2latlon(P_cart_sp);
            sp_Q      = cart2latlon(Q_cart_sp);
        end
        
        %%%%%%%%%%%%%%%%%%%%%
        %%% Vertical node %%%
        %%%%%%%%%%%%%%%%%%%%%
    else
        disp(['Starting node is vertical. Row is ',num2str(r_sp),', column is ',num2str(c_sp)])
        if c_sp == 1
            Q_tail      = Q(r_sp,:);
            Q_head      = Q(r_sp+1,:);
            Q_tail_cart = latlon2cart(Q_tail,1);
            Q_head_cart = latlon2cart(Q_head,1);
            best_f_Q      = Lmin(r_sp,c_sp);
            best_f_P      = 0;
            unit_square_Q = r_sp;
            unit_square_P = c_sp;
            DQ          = Q_length(r_sp);
            Q_cart_sp   = Q_tail_cart*sin((1-best_f_Q)*DQ)/sin(DQ)+Q_head_cart*sin(best_f_Q*DQ)/sin(DQ);
            sp_P        = P_red(c_sp,:);
            sp_Q        = cart2latlon(Q_cart_sp);
        else
            f_vec_P     = linspace(0,1,N_f);
            f_vec_Q     = linspace(0,Lmax(r_sp,c_sp),N_f);
            P_tail      = P_red(c_sp-1,:);
            P_head      = P_red(c_sp,:);
            Q_tail      = Q(r_sp,:);
            Q_head      = Q(r_sp+1,:);
            P_tail_cart = latlon2cart(P_tail,1);
            P_head_cart = latlon2cart(P_head,1);
            Q_tail_cart = latlon2cart(Q_tail,1);
            Q_head_cart = latlon2cart(Q_head,1);
            DP          = P_length(c_sp-1);
            DQ          = Q_length(r_sp);
            dist_mat    = zeros(numel(f_vec_Q),numel(f_vec_P));
            dist_tail       = zeros(size(dist_mat));
            for ii=1:numel(f_vec_Q)
                Q_cart   = Q_tail_cart*sin((1-f_vec_Q(ii))*DQ)/sin(DQ)+Q_head_cart*sin(f_vec_Q(ii)*DQ)/sin(DQ);
                Q_ll     = cart2latlon(Q_cart);
                for jj=1:numel(f_vec_P)
                    P_cart   = P_tail_cart*sin((1-f_vec_P(jj))*DP)/sin(DP)+P_head_cart*sin(f_vec_P(jj)*DP)/sin(DP);
                    P_ll     = cart2latlon(P_cart);
                    dist_mat(ii,jj) = Haversine(P_ll,Q_ll,1);
                    dist_tail(ii,jj)    = Haversine(P_ll,P_head,1)+Haversine(Q_ll,Q_head,1);
                end
            end
            [r,c]    = find(dist_mat<=alpha*FD);
            dist_tail_vec = zeros(numel(r),1);
            for ii=1:numel(r)
                dist_tail_vec(ii) = dist_tail(r(ii),c(ii));
            end
            [~,idx] = max(dist_tail_vec);
            best_f_Q      = f_vec_Q(r(idx));
            best_f_P      = f_vec_P(c(idx));
            unit_square_Q = r_sp;
            unit_square_P = c_sp-1;
            P_cart_sp = P_tail_cart*sin((1-best_f_P)*DP)/sin(DP)+P_head_cart*sin(best_f_P*DP)/sin(DP);
            Q_cart_sp = Q_tail_cart*sin((1-best_f_Q)*DQ)/sin(DQ)+Q_head_cart*sin(best_f_Q*DQ)/sin(DQ);
            sp_P      = cart2latlon(P_cart_sp);
            sp_Q      = cart2latlon(Q_cart_sp);
        end
    end
    
    Overlap_points2(i,1:8) = horzcat(unit_square_P,best_f_P,sp_P,unit_square_Q,best_f_Q,sp_Q);
    
    %%%%%%%%%%%%%%%%%
    %%%  Endpoint %%%
    %%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%
    %%% Corner node %%%
    %%%%%%%%%%%%%%%%%%%
    if Type_ep == 1
        disp(['Ending node is corner. Row is ',num2str(r_ep),', column is ',num2str(c_ep)])
        if r_ep == q+1 && c_ep == p+1
            ep_P          = P_red(c_ep,:);
            ep_Q          = Q(r_ep,:);
            best_f_Q      = 1;
            best_f_P      = 1;
            unit_square_Q = q;
            unit_square_P = p;
        else
            f_vec_P     = linspace(0,1,N_f);
            f_vec_Q     = linspace(0,1,N_f);
            if r_ep == q+1
                disp('Case 1')
                P_tail      = P_red(c_ep,:)
                P_head      = P_red(c_ep+1,:)
                Q_tail      = Q(q,:)
                Q_head      = Q(q+1,:)
                DP          = P_length(c_ep);
                DQ          = Q_length(q);
                unit_square_Q = q
                unit_square_P = c_ep
                P_tail_cart = latlon2cart(P_tail,1);
            P_head_cart = latlon2cart(P_head,1);
            Q_tail_cart = latlon2cart(Q_tail,1);
            Q_head_cart = latlon2cart(Q_head,1);
            dist_mat    = zeros(numel(f_vec_Q),numel(f_vec_P));
            dist_tail   = zeros(size(dist_mat));
            for ii=1:numel(f_vec_Q)
                Q_cart   = Q_tail_cart*sin((1-f_vec_Q(ii))*DQ)/sin(DQ)+Q_head_cart*sin(f_vec_Q(ii)*DQ)/sin(DQ);
                for jj=1:numel(f_vec_P)
                    P_cart   = P_tail_cart*sin((1-f_vec_P(jj))*DP)/sin(DP)+P_head_cart*sin(f_vec_P(jj)*DP)/sin(DP);
                    P_ll     = cart2latlon(P_cart);
                    Q_ll     = cart2latlon(Q_cart);
                    dist_mat(ii,jj) = Haversine(P_ll,Q_ll,1);
                    dist_tail(ii,jj)    = Haversine(P_tail,P_ll,1)+Haversine(Q_head,Q_ll,1);
                end
            end
            [r,c]    = find(dist_mat<=alpha*FD);
            dist_tail_vec = zeros(numel(r),1);
            for ii=1:numel(r)
                dist_tail_vec(ii) = dist_tail(r(ii),c(ii));
            end
            [~,idx] = max(dist_tail_vec);
            best_f_Q      = f_vec_Q(r(idx));
            best_f_P      = f_vec_P(c(idx));
            P_cart_ep = P_tail_cart*sin((1-best_f_P)*DP)/sin(DP)+P_head_cart*sin(best_f_P*DP)/sin(DP);
            Q_cart_ep = Q_tail_cart*sin((1-best_f_Q)*DQ)/sin(DQ)+Q_head_cart*sin(best_f_Q*DQ)/sin(DQ);
            ep_P      = cart2latlon(P_cart_ep);
            ep_Q      = cart2latlon(Q_cart_ep);
            elseif c_ep == p+1
                disp('Case 2')
                P_tail      = P_red(c_ep-1,:)
                P_head      = P_red(c_ep,:)
                Q_tail      = Q(r_ep,:)
                Q_head      = Q(r_ep+1,:)
                DP          = P_length(p);
                DQ          = Q_length(r_ep);
                unit_square_Q = r_ep
                unit_square_P = p
                 P_tail_cart = latlon2cart(P_tail,1);
            P_head_cart = latlon2cart(P_head,1);
            Q_tail_cart = latlon2cart(Q_tail,1);
            Q_head_cart = latlon2cart(Q_head,1);
            dist_mat    = zeros(numel(f_vec_Q),numel(f_vec_P));
            dist_tail   = zeros(size(dist_mat),1);
            for ii=1:numel(f_vec_Q)
                Q_cart   = Q_tail_cart*sin((1-f_vec_Q(ii))*DQ)/sin(DQ)+Q_head_cart*sin(f_vec_Q(ii)*DQ)/sin(DQ);
                for jj=1:numel(f_vec_P)
                    P_cart   = P_tail_cart*sin((1-f_vec_P(jj))*DP)/sin(DP)+P_head_cart*sin(f_vec_P(jj)*DP)/sin(DP);
                    P_ll     = cart2latlon(P_cart);
                    Q_ll     = cart2latlon(Q_cart);
                    dist_mat(ii,jj) = Haversine(P_ll,Q_ll,1);
                    dist_tail(ii,jj)    = Haversine(P_tail,P_ll,1)+Haversine(Q_tail,Q_ll,1);
                end
            end
            [r,c]    = find(dist_mat<=alpha*FD);
            dist_tail_vec = zeros(numel(r),1);
            for ii=1:numel(r)
                dist_tail_vec(ii) = dist_tail(r(ii),c(ii));
            end
            [~,idx] = max(dist_tail_vec);
            best_f_Q      = f_vec_Q(r(idx));
            best_f_P      = f_vec_P(c(idx));
            P_cart_ep = P_tail_cart*sin((1-best_f_P)*DP)/sin(DP)+P_head_cart*sin(best_f_P*DP)/sin(DP);
            Q_cart_ep = Q_tail_cart*sin((1-best_f_Q)*DQ)/sin(DQ)+Q_head_cart*sin(best_f_Q*DQ)/sin(DQ);
            ep_P      = cart2latlon(P_cart_ep);
            ep_Q      = cart2latlon(Q_cart_ep);
            else
                disp('Case 3')
                P_tail      = P_red(c_ep,:)
                P_head      = P_red(c_ep+1,:)
                Q_tail      = Q(r_ep,:)
                Q_head      = Q(r_ep+1,:)
                DP          = P_length(c_ep);
                DQ          = Q_length(r_ep);
                unit_square_Q = r_ep
                unit_square_P = c_ep
                 P_tail_cart = latlon2cart(P_tail,1);
            P_head_cart = latlon2cart(P_head,1);
            Q_tail_cart = latlon2cart(Q_tail,1);
            Q_head_cart = latlon2cart(Q_head,1);
            dist_mat    = zeros(numel(f_vec_Q),numel(f_vec_P));
            dist_tail   = zeros(size(dist_mat));
            for ii=1:numel(f_vec_Q)
                Q_cart   = Q_tail_cart*sin((1-f_vec_Q(ii))*DQ)/sin(DQ)+Q_head_cart*sin(f_vec_Q(ii)*DQ)/sin(DQ);
                for jj=1:numel(f_vec_P)
                    P_cart   = P_tail_cart*sin((1-f_vec_P(jj))*DP)/sin(DP)+P_head_cart*sin(f_vec_P(jj)*DP)/sin(DP);
                    P_ll     = cart2latlon(P_cart);
                    Q_ll     = cart2latlon(Q_cart);
                    dist_mat(ii,jj) = Haversine(P_ll,Q_ll,1);
                    dist_tail(ii,jj)    = Haversine(P_tail,P_ll,1)+Haversine(Q_tail,Q_ll,1);
                end
            end
            [r,c]    = find(dist_mat<=alpha*FD);
            dist_tail_vec = zeros(numel(r),1);
            for ii=1:numel(r)
                dist_tail_vec(ii) = dist_tail(r(ii),c(ii));
            end
            [~,idx] = max(dist_tail_vec);
            best_f_Q      = f_vec_Q(r(idx));
            best_f_P      = f_vec_P(c(idx));
            P_cart_ep = P_tail_cart*sin((1-best_f_P)*DP)/sin(DP)+P_head_cart*sin(best_f_P*DP)/sin(DP);
            Q_cart_ep = Q_tail_cart*sin((1-best_f_Q)*DQ)/sin(DQ)+Q_head_cart*sin(best_f_Q*DQ)/sin(DQ);
            ep_P      = cart2latlon(P_cart_ep);
            ep_Q      = cart2latlon(Q_cart_ep);
            end
%             P_tail_cart = latlon2cart(P_tail,1);
%             P_head_cart = latlon2cart(P_head,1);
%             Q_tail_cart = latlon2cart(Q_tail,1);
%             Q_head_cart = latlon2cart(Q_head,1);
%             dist_mat    = zeros(numel(f_vec_Q),numel(f_vec_P));
%             dist_tail   = zeros(size(dist_mat));
%             for ii=1:numel(f_vec_Q)
%                 Q_cart   = Q_tail_cart*sin((1-f_vec_Q(ii))*DQ)/sin(DQ)+Q_head_cart*sin(f_vec_Q(ii)*DQ)/sin(DQ);
%                 for jj=1:numel(f_vec_P)
%                     P_cart   = P_tail_cart*sin((1-f_vec_P(jj))*DP)/sin(DP)+P_head_cart*sin(f_vec_P(jj)*DP)/sin(DP);
%                     P_ll     = cart2latlon(P_cart);
%                     Q_ll     = cart2latlon(Q_cart);
%                     dist_mat(ii,jj) = Haversine(P_ll,Q_ll,1);
%                     dist_tail(ii,jj)    = Haversine(P_tail,P_ll,1)+Haversine(Q_tail,Q_ll,1);
%                 end
%             end
%             [r,c]    = find(dist_mat<=alpha*FD);
%             dist_tail_vec = zeros(numel(r));
%             for ii=1:numel(r)
%                 dist_tail_vec(ii) = dist_tail(r(ii),c(ii));
%             end
%             [~,idx] = max(dist_tail_vec);
%             best_f_Q      = f_vec_Q(r(idx));
%             best_f_P      = f_vec_P(c(idx));
%             P_cart_ep = P_tail_cart*sin((1-best_f_P)*DP)/sin(DP)+P_head_cart*sin(best_f_P*DP)/sin(DP);
%             Q_cart_ep = Q_tail_cart*sin((1-best_f_Q)*DQ)/sin(DQ)+Q_head_cart*sin(best_f_Q*DQ)/sin(DQ);
%             ep_P      = cart2latlon(P_cart_ep);
%             ep_Q      = cart2latlon(Q_cart_ep);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%
        %%% Horizontal node %%%
        %%%%%%%%%%%%%%%%%%%%%%%
    elseif Type_ep == 2
        disp(['Ending node is horizontal. Row is ',num2str(r_ep),', column is ',num2str(c_ep)])
        if r_ep == q+1
            P_tail      = P_red(c_ep,:);
            P_head      = P_red(c_ep+1,:);
            P_tail_cart = latlon2cart(P_tail,1);
            P_head_cart = latlon2cart(P_head,1);
            best_f_Q      = 1;
            best_f_P      = Bmin(r_ep,c_ep);
            unit_square_Q = q;
            unit_square_P = c_ep;
            DP          = P_length(c_ep);
            P_cart_ep   = P_tail_cart*sin((1-best_f_P)*DP)/sin(DP)+P_head_cart*sin(best_f_P*DP)/sin(DP);
            ep_P        = cart2latlon(P_cart_ep);
            ep_Q        = Q(r_ep,:);
        else
            f_vec_P     = linspace(Bmin(r_ep,c_ep),1,N_f);
            f_vec_Q     = linspace(0,1,N_f);
            P_tail      = P_red(c_ep,:);
            P_head      = P_red(c_ep+1,:);
            Q_tail      = Q(r_ep,:);
            Q_head      = Q(r_ep+1,:);
            P_tail_cart = latlon2cart(P_tail,1);
            P_head_cart = latlon2cart(P_head,1);
            Q_tail_cart = latlon2cart(Q_tail,1);
            Q_head_cart = latlon2cart(Q_head,1);
            DP          = P_length(c_ep);
            DQ          = Q_length(r_ep);
            dist_mat    = zeros(numel(f_vec_Q),numel(f_vec_P));
            dist_tail       = zeros(size(dist_mat));
            for ii=1:numel(f_vec_Q)
                Q_cart   = Q_tail_cart*sin((1-f_vec_Q(ii))*DQ)/sin(DQ)+Q_head_cart*sin(f_vec_Q(ii)*DQ)/sin(DQ);
                for jj=1:numel(f_vec_P)
                    P_cart   = P_tail_cart*sin((1-f_vec_P(jj))*DP)/sin(DP)+P_head_cart*sin(f_vec_P(jj)*DP)/sin(DP);
                    P_ll     = cart2latlon(P_cart);
                    Q_ll     = cart2latlon(Q_cart);
                    dist_mat(ii,jj) = Haversine(P_ll,Q_ll,1);
                    dist_tail(ii,jj)    = Haversine(P_tail,P_ll,1)+Haversine(Q_tail,Q_ll,1);
                end
            end
            [r,c]    = find(dist_mat<=alpha*FD);
            dist_tail_vec = zeros(numel(r),1);
            for ii=1:numel(r)
                dist_tail_vec(ii) = dist_tail(r(ii),c(ii));
            end
            [~,idx] = max(dist_tail_vec);
            best_f_Q      = f_vec_Q(r(idx));
            best_f_P      = f_vec_P(c(idx));
            unit_square_Q = r_ep;
            unit_square_P = c_ep;
            P_cart_ep = P_tail_cart*sin((1-best_f_P)*DP)/sin(DP)+P_head_cart*sin(best_f_P*DP)/sin(DP);
            Q_cart_ep = Q_tail_cart*sin((1-best_f_Q)*DQ)/sin(DQ)+Q_head_cart*sin(best_f_Q*DQ)/sin(DQ);
            ep_P      = cart2latlon(P_cart_ep);
            ep_Q      = cart2latlon(Q_cart_ep);
        end
        
        %%%%%%%%%%%%%%%%%%%%%
        %%% Vertical node %%%
        %%%%%%%%%%%%%%%%%%%%%
    else
        disp(['Ending node is vertical. Row is ',num2str(r_ep),', column is ',num2str(c_ep)])
        if c_ep == p+1
            Q_tail      = Q(r_ep,:);
            Q_head      = Q(r_ep+1,:);
            Q_tail_cart = latlon2cart(Q_tail,1);
            Q_head_cart = latlon2cart(Q_head,1);
            best_f_Q      = Lmin(r_ep,c_ep);
            best_f_P      = 1;
            unit_square_Q = r_ep;
            unit_square_P = p;
            DQ          = Q_length(r_ep);
            Q_cart_ep   = Q_tail_cart*sin((1-best_f_Q)*DQ)/sin(DQ)+Q_head_cart*sin(best_f_Q*DQ)/sin(DQ);
            ep_P        = P_red(c_ep,:);
            ep_Q        = cart2latlon(Q_cart_ep);
        else
            f_vec_P     = linspace(0,1,N_f);
            f_vec_Q     = linspace(Lmin(r_ep,c_ep),1,N_f);
            P_tail      = P_red(c_ep,:);
            P_head      = P_red(c_ep+1,:);
            Q_tail      = Q(r_ep,:);
            Q_head      = Q(r_ep+1,:);
            DP          = P_length(c_ep);
            DQ          = Q_length(r_ep);
            
            P_tail_cart = latlon2cart(P_tail,1);
            P_head_cart = latlon2cart(P_head,1);
            Q_tail_cart = latlon2cart(Q_tail,1);
            Q_head_cart = latlon2cart(Q_head,1);
            dist_mat    = zeros(numel(f_vec_Q),numel(f_vec_P));
            dist_tail       = zeros(size(dist_mat));
            for ii=1:numel(f_vec_Q)
                Q_cart   = Q_tail_cart*sin((1-f_vec_Q(ii))*DQ)/sin(DQ)+Q_head_cart*sin(f_vec_Q(ii)*DQ)/sin(DQ);
                Q_ll     = cart2latlon(Q_cart);
                for jj=1:numel(f_vec_P)
                    P_cart   = P_tail_cart*sin((1-f_vec_P(jj))*DP)/sin(DP)+P_head_cart*sin(f_vec_P(jj)*DP)/sin(DP);
                    P_ll     = cart2latlon(P_cart);
                    dist_mat(ii,jj) = Haversine(P_ll,Q_ll,1);
                    dist_tail(ii,jj)    = Haversine(P_tail,P_ll,1)+Haversine(Q_tail,Q_ll,1);
                end
            end
            [r,c]    = find(dist_mat<=alpha*FD);
            dist_tail_vec = zeros(numel(r),1);
            for ii=1:numel(r)
                dist_tail_vec(ii) = dist_tail(r(ii),c(ii));
            end
            [~,idx] = max(dist_tail_vec);
            best_f_Q      = f_vec_Q(r(idx));
            best_f_P      = f_vec_P(c(idx));
            unit_square_Q = r_ep;
            unit_square_P = c_ep;
            P_cart_ep = P_tail_cart*sin((1-best_f_P)*DP)/sin(DP)+P_head_cart*sin(best_f_P*DP)/sin(DP);
            Q_cart_ep = Q_tail_cart*sin((1-best_f_Q)*DQ)/sin(DQ)+Q_head_cart*sin(best_f_Q*DQ)/sin(DQ);
            ep_P      = cart2latlon(P_cart_ep);
            ep_Q      = cart2latlon(Q_cart_ep);
        end
        
    end
    
    Overlap_points(i,:)     = horzcat(sp_P,sp_Q,ep_P,ep_Q);
    Overlap_points2(i,9:16) = horzcat(unit_square_P,best_f_P,ep_P,unit_square_Q,best_f_Q,ep_Q);
    
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    
end

time5 = toc;

tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Check is some intersection points were not recognized by the %%%
%%% algorithm, and detect the associated overlapping region if   %%%
%%% at least one point was not recognized                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[r_im,c_im] = find(int_matrix);


if ~isempty(r_im)
    
    check_int   = zeros(numel(r_im),1);
    
    if isempty(Merged_Sequences)
        check_int = zeros(numel(r_im),1);
    else
        
        
        for i=1:length(r_im)
            r_int        = r_im(i);
            c_int        = c_im(i);
            cont         = 1;
            is_contained = 0;
            while is_contained~=1
                % Intersection is already accounted for in one of the
                % overlapping sequences that we have identified
                if r_int>=Merged_Sequences_rc(cont,1) && ...
                        r_int<=Merged_Sequences_rc(cont,2) && ...
                        c_int>=Merged_Sequences_rc(cont,3) && ...
                        c_int<=Merged_Sequences_rc(cont,4)
                    disp(['Intersection in unitary square ',num2str(r_int),...
                        '-',num2str(c_int),' already accounted for'])
                    is_contained = 1;
                    break;
                else
                    cont = cont+1;
                    if cont == length(Merged_Sequences_rc(:,1))+1
                        disp(['Intersection in unitary square ',num2str(r_int),...
                            '-',num2str(c_int),' needs to be computed'])
                        check_int(i) = 1;
                        is_contained = 1;
                    else
                    end
                end
            end
        end
        
    end
end

if ~isempty(find(check_int,1))
    idx                 = find(check_int,1);
    new_Overlap_points2 = zeros(numel(idx),16);
    for i=1:numel(idx)
        f_vec_P     = linspace(0,1,N_f);
        f_vec_Q     = linspace(0,1,N_f);
        P_tail      = P_red(c_im(idx(i)),:);
        P_head      = P_red(c_im(idx(i))+1,:);
        Q_tail      = Q(r_im(idx(i)),:);
        Q_head      = Q(r_im(idx(i))+1,:);
        P_tail_cart = latlon2cart(P_tail,1);
        P_head_cart = latlon2cart(P_head,1);
        Q_tail_cart = latlon2cart(Q_tail,1);
        Q_head_cart = latlon2cart(Q_head,1);
        DP          = P_length(c_im(idx(i)));
        DQ          = Q_length(r_im(idx(i)));
        dist_mat    = zeros(numel(f_vec_Q),numel(f_vec_P));
        dist_tail   = zeros(size(dist_mat));
        dist_head   = zeros(size(dist_mat));
        for ii=1:numel(f_vec_Q)
            Q_cart   = Q_tail_cart*sin((1-f_vec_Q(ii))*DQ)/sin(DQ)+Q_head_cart*sin(f_vec_Q(ii)*DQ)/sin(DQ);
            Q_ll     = cart2latlon(Q_cart);
            for jj=1:numel(f_vec_P)
                P_cart   = P_tail_cart*sin((1-f_vec_P(jj))*DP)/sin(DP)+P_head_cart*sin(f_vec_P(jj)*DP)/sin(DP);
                P_ll     = cart2latlon(P_cart);
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
        [~,idx1]    = min(dist_tail_vec);
        [~,idx2]    = min(dist_head_vec);
        best_f_Q_1  = f_vec_Q(r(idx1));
        best_f_P_1  = f_vec_P(c(idx1));
        P_cart_ep_1 = P_tail_cart*sin((1-best_f_P_1)*DP)/sin(DP)+P_head_cart*sin(best_f_P_1*DP)/sin(DP);
        Q_cart_ep_1 = Q_tail_cart*sin((1-best_f_Q_1)*DQ)/sin(DQ)+Q_head_cart*sin(best_f_Q_1*DQ)/sin(DQ);
        sp_P        = cart2latlon(P_cart_ep_1);
        sp_Q        = cart2latlon(Q_cart_ep_1);
        best_f_Q_2  = f_vec_Q(r(idx2));
        best_f_P_2  = f_vec_P(c(idx2));
        P_cart_ep_2 = P_tail_cart*sin((1-best_f_P_2)*DP)/sin(DP)+P_head_cart*sin(best_f_P_2*DP)/sin(DP);
        Q_cart_ep_2 = Q_tail_cart*sin((1-best_f_Q_2)*DQ)/sin(DQ)+Q_head_cart*sin(best_f_Q_2*DQ)/sin(DQ);
        ep_P        = cart2latlon(P_cart_ep_2);
        ep_Q        = cart2latlon(Q_cart_ep_2);
        
        new_Overlap_points2(i,:) = horzcat(c_im(idx(i)),best_f_P_1,sp_P,r_im(idx(i)),best_f_Q_1,sp_Q,c_im(idx(i)),best_f_P_2,ep_P,r_im(idx(i)),best_f_Q_2,ep_Q);
        
    end
else
end

time6 = toc;

tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 if ~isempty(Overlap_points2) && exist('new_Overlap_points2','var')
     OVERLAP = vertcat(Overlap_points2,new_Overlap_points2);
     over_id = 1;
 elseif ~isempty(Overlap_points2) && ~exist('new_Overlap_points2','var')
     OVERLAP = Overlap_points2;
     over_id = 1;
 elseif isempty(Overlap_points2) && exist('new_Overlap_points2','var')
     OVERLAP = new_Overlap_points2;
     over_id = 1;
 else
     disp('No common regions compatible with the threshold were found')
 end
 
 if over_id
     this_length = zeros(length(OVERLAP(:,1)),4);
     all_P       = sum(P_length);
     all_Q       = sum(Q_length);
     for i=1:length(OVERLAP(:,1))
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %%% Overlapping region on P %%%
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
         % Overlapping region lies within a single great circle of P
         if OVERLAP(i,9)-OVERLAP(i,1) == 0
             this_length(i,1) = Haversine([OVERLAP(i,3),OVERLAP(i,4)],[OVERLAP(i,11),OVERLAP(i,12)],1);
             this_length(i,2) = this_length(i,1)/all_P; 
         % Overlapping region lies within 2 contiguous great circles of P
         elseif OVERLAP(i,9)-OVERLAP(i,1) == 1
             this_length(i,1) = Haversine([OVERLAP(i,3),OVERLAP(i,4)],P_red(OVERLAP(i,1)+1,:),1)+...
                           Haversine(P_red(OVERLAP(i,1)+1,:),[OVERLAP(i,11),OVERLAP(i,12)],1);
             this_length(i,2) = this_length(i,1)/all_P;          
         % Overlapping regions extends over 3 or more great circles of P    
         else
             this_length(i,1) = Haversine([OVERLAP(i,3),OVERLAP(i,4)],P_red(OVERLAP(i,1)+1,:),1)+...
                           sum(P_length(OVERLAP(i,1)+1:OVERLAP(i,9)-1))+...
                           Haversine(P_red(OVERLAP(i,9),:),[OVERLAP(i,11),OVERLAP(i,12)],1);
             this_length(i,2) = this_length(i,1)/all_P;          
         end
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %%% Overlapping region on Q %%%
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
         % Overlapping region lies within a single great circle of Q
         if OVERLAP(i,13)-OVERLAP(i,5) == 0
             this_length(i,3) = Haversine([OVERLAP(i,7),OVERLAP(i,8)],[OVERLAP(i,15),OVERLAP(i,16)],1);
             this_length(i,4) = this_length(i,3)/all_Q;
         % Overlapping region lies within 2 contiguous great circles of Q
         elseif OVERLAP(i,13)-OVERLAP(i,5) == 1
             this_length(i,3) = Haversine([OVERLAP(i,7),OVERLAP(i,8)],Q(OVERLAP(i,5)+1,:),1)+...
                           Haversine(Q(OVERLAP(i,5)+1,:),[OVERLAP(i,15),OVERLAP(i,16)],1);
             this_length(i,4) = this_length(i,3)/all_Q;          
         % Overlapping regions extends over 3 or more great circles of Q    
         else
             this_length(i,3) = Haversine([OVERLAP(i,7),OVERLAP(i,8)],Q(OVERLAP(i,5)+1,:),1)+...
                           sum(Q_length(OVERLAP(i,5)+1:OVERLAP(i,13)-1))+...
                           Haversine(Q(OVERLAP(i,13),:),[OVERLAP(i,15),OVERLAP(i,16)],1);
             this_length(i,4) = this_length(i,3)/all_Q;          
         end
     end
 else
 end
 
 time7 = toc;
 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if over_id
    
    N = 20;
    figure()
    hold on
    for i=1:numel(P_red(:,1))-1
        [lat,lon] = gcwaypts(P_red(i,1),P_red(i,2),P_red(i+1,1),P_red(i+1,2),N);
        h1 = plot(lon,lat,'Color','b','Linewidth',2,'Marker','none','Markersize',8);
    end
    for i=1:numel(Q(:,1))-1
        [lat,lon] = gcwaypts(Q(i,1),Q(i,2),Q(i+1,1),Q(i+1,2),N);
        h2 = plot(lon,lat,'Color','k','Linewidth',2,'Marker','none','Markersize',8);
    end
    plot(P_red(:,2),P_red(:,1),'Color','b','Linewidth',2,'Linestyle','none','Marker','*','Markersize',6)
    plot(Q(:,2),Q(:,1),'Color','k','Linewidth',2,'Linestyle','none','Marker','*','Markersize',6)
    plot(OVERLAP(:,4),OVERLAP(:,3),'Color','g','Linewidth',2.5,'Linestyle','none','Marker','x','Markersize',10)
    plot(OVERLAP(:,8),OVERLAP(:,7),'Color','g','Linewidth',2.5,'Linestyle','none','Marker','x','Markersize',10)
    plot(OVERLAP(:,12),OVERLAP(:,11),'Color','r','Linewidth',2.5,'Linestyle','none','Marker','x','Markersize',10)
    plot(OVERLAP(:,16),OVERLAP(:,15),'Color','r','Linewidth',2.5,'Linestyle','none','Marker','x','Markersize',10)
    xlabel('Longitude [deg]','Fontname','Avantgarde','Fontsize',14)
    ylabel('Latitude [deg]','Fontname','Avantgarde','Fontsize',14)
    l1 = legend([h1,h2],'P polygonal curve','Q polygonal curve');
    set(l1,'Fontname','Avantgarde','Fontsize',14)
    grid on
    
else
end



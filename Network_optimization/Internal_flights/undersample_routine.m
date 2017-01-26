% Written by Alessandro Bombelli, May 29th 2016
% Comparing performance of undersampling routine for frecher Distance 
% Computation 

clc
clear all
close all

% Flights = load(strcat(pwd,'/NETWORK/INT_traj_filt_info/65/65_34_int.txt'));
% 
% flight1  = Flights(Flights(:,1)==1,:);
% flight2  = Flights(Flights(:,1)==100,:);
% heading1 = traj_heading(flight1(:,4:5));
% heading2 = traj_heading(flight2(:,4:5));
% 
% max_delta_heading   = 7;
% delta_heading1      = abs(heading1(2:end)-heading1(1:end-1));
% dummy1              = delta_heading1<=max_delta_heading;
% idx1                = find(dummy1==1)+1;
% flight1_red         = flight1;
% flight1_red(idx1,:) = [];
% delta_heading2      = abs(heading2(2:end)-heading2(1:end-1));
% dummy2              = delta_heading2<=max_delta_heading;
% idx2                = find(dummy2==1)+1;
% flight2_red         = flight2;
% flight2_red(idx2,:) = [];
% 
% Re                           = 1;
% eps                          = 0.001;
% tolerance                    = 1;
% 
% disp('Begin case 1')
% tic
% [FD,no_sol,sol] = FD_computation_v02(flight1(:,4:5),flight2(:,4:5),Re,eps,tolerance);
% time = toc;
% 
% disp('Begin case 2')
% tic
% [no_solnew,solnew,cd_allnew,FDnew] = FD_computation_v03(flight1(:,4:5),flight2(:,4:5),Re);
% timenew = toc;
% 
% disp('Begin case 3')
% tic
% [FD2,no_sol2,sol2] = FD_computation_v02(flight1_red(:,4:5),flight2_red(:,4:5),Re,eps,0.0001);
% time2 = toc;
% 
% disp('Begin case 4')
% tic
% [cd_all,FD,sol,no_sol] = FD_computation_v03(flight1_red(:,4:5),flight2_red(:,4:5),Re);
% timenew2 = toc;
% 
% disp('Begin case 5')
% tic
% [cm, cSq, cD] = func_discrete_Frechet_lat_lon_NM_v02(flight1(:,4:5),flight2(:,4:5));
% time3 = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Re                           = 1;
eps                          = 0.001;
tolerance                    = 1;

% A = load(strcat(pwd,'/TRX_PROCESSED/2014070100_2014070124_1of3.mat'));
% 
% names = fieldnames(A.ACID_info);
% 
% for i=1:length(names)
%     flight      = A.ACID_info.(names{i});
%     origin      = flight.origin;
%     destination = flight.destination;
%     if strcmp(origin,'KJFK')==1
%         if strcmp(destination,'KLAX')==1
%             disp(['Flight ',num2str(i),' is good'])
%         else
%         end
%     else
%     end
% end
% 
%  P = A.ACID_info.(names{146}).lat_lon;
%  Q = [40.6398 -73.7789;40.2 -74.5;40 -76.5;38.33 -81.8;36.1 -86.65;35 -89.95;...
%       34.66 -92.16;35.33 -97.6;35.25 -101.6;34.6 -112.5;34.05 -115.7;...
%       33.95 -117.5;33.9425 -118.4081];
  
  

% P = A.ACID_info.DAL2263_UID01.lat_lon;
% Q = [40.788 -111.98;41.5 -109; 41.8 -106; 41.9 -103.3;42.5 -99.9;43.6 -96.7;44.45 -95.09;44.64 -94.37;44.88 -93.22];
%

P = [0.1 -10;0.1 -8;0.15 -6;0.25 -5;0.3 -4;0.25 -3;-0.05 -2;-0.1 -1;-0.15 0;-0.05 2;-0.03 4;0.1 5;0.2 6;0.15 7;0.1 9; 0.05 10];
Q = horzcat(zeros(11,1),(-10:2:10)'); 

%headingP = traj_heading(P);

%max_delta_heading   = 4;
%delta_heading1      = abs(headingP(2:end)-headingP(1:end-1));
%dummy1              = delta_heading1<=max_delta_heading;
%idx1                = find(dummy1==1)+1;
%P_red               = P;
%P_red(idx1,:)       = [];

P_red = P;

% figure()
% hold on
% plot(1:numel(headingP),headingP,'Color','b','Linewidth',2,'Marker','*','Markersize',8)
% grid on
% 
% figure()
% hold on
% plot(P(:,2),P(:,1),'Color','b','Linewidth',2,'Marker','*','Markersize',8)
% plot(Q(:,2),Q(:,1),'Color','r','Linewidth',2,'Marker','*','Markersize',8)
% grid on
% 
% disp('Begin case 1')
% tic
% [cd_all,FD,sol,no_sol] = FD_computation_v03(P,Q,1);
% timenew2 = toc;
% 
%disp('Begin case 2')
%tic
%[FD2,no_sol2,sol2] = FD_computation_v02(P_red,Q,Re,0.001,0.0001);
%time2 = toc;
% 
% disp('Begin case 3')
% tic
% [cm, cSq, cD] = func_discrete_Frechet_lat_lon_NM_v02(P,Q);
% time3 = toc;

disp('Begin case 4')
tic
[cd_all,FD,sol,no_sol,Lmin,Lmax,Bmin,Bmax,A,sequence,Node_ID] = FD_computation_v03(P_red,Q,1);
time4 = toc;

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

[Lmin_eps,Lmax_eps,Bmin_eps,Bmax_eps,A_list_eps,Node_ID_eps] = free_space_shape_v02(P_red,Q,Re,FD/3);
plot_FS(Bmin_eps,Bmax_eps,Lmin_eps,Lmax_eps)

[P_length,Q_length,dist_vector,NODES] = overlap_index_v01(P_red,Q,sequence,Node_ID,Lmin,Lmax,Bmin,Bmax);
idx                                   = find(dist_vector<=FD/3);
seq2                                  = NODES(idx(2:end),1)-NODES(idx(1:end-1),1);
check_vector                          = 1;
cont                                  = 1;
BLOCKS                                = cell(0,0);
cont_block                            = 1;

while check_vector~=0
    
    % the whole vector has been scanned
    if cont >= numel(seq2)
        check_vector = 0;
    else
    end
    
    % just increase the counter
    if seq2(cont)~=1
        cont = cont+1;
        % we have a valid sequence
    else
        block = [];
        
        if cont == numel(seq2)
            block = cont;
            cont  = cont+1;
            
        elseif seq2(cont+1)-seq2(cont)~=0
            block = cont;
            cont = cont+1;
        else
            block     = cont;
            this_cont = cont;
            counter   = 1;
            next      = 1;
            while next==1
                block   = [block this_cont+counter];
                counter = counter+1;
                if this_cont+counter>numel(seq2)
                    next         = 0;
                    check_vector = 0;
                elseif seq2(this_cont+counter)~=1
                    next = 0;
                else
                end
            end
            cont = this_cont+counter;
        end
        BLOCKS{cont_block,1} = block;
        cont_block           = cont_block+1;
    end
end

if ~isempty(BLOCKS)
    mu_R      = 3;
    mu_C      = 5;
    p         = length(P_red(:,1))-1;
    q         = length(Q(:,1))-1;
    Sequences = [];
    cont      = 1;
    
    for i=1:length(BLOCKS)
        this_sequence      = BLOCKS{i,1};
        % Initial startpoint and endpoint used as inputs for the DFS
        initial_startpoint      = NODES(idx(this_sequence(1)),2);
        initial_startpoint_type = NODES(idx(this_sequence(1)),3);
        initial_startpoint_R    = NODES(idx(this_sequence(1)),4);
        initial_startpoint_C    = NODES(idx(this_sequence(1)),5);
        initial_endpoint        = NODES(idx(this_sequence(end)+1),2);
        initial_endpoint_type   = NODES(idx(this_sequence(end)+1),3);
        initial_endpoint_R      = NODES(idx(this_sequence(end)+1),4);
        initial_endpoint_C      = NODES(idx(this_sequence(end)+1),5);
        
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
                    disp(['Path between node ',num2str(sp),' and node ',num2str(ep)])
                    Node_ID_eps(sp,:)
                    Node_ID_eps(ep,:)
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
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    end
else
end

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% For each sequence, translate the starting and endpoint of the %%%
%%% sequence into the corresponding points on curves P and Q, and %%%
%%% store the values in a new matrix                              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Overlap_points = zeros(length(Merged_Sequences(:,1)),12);

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
    
    % Corner node
    if Type_sp == 1
        sp_P = P_red(c_sp,:);
        sp_Q = Q(r_sp,:);
        % Horizontal node
    elseif Type_sp == 2
        % Select the minimum value for f, since this is the starting point
        f           = Bmin(r_sp,c_sp);
        P_tail      = P_red(c_sp,:);
        P_head      = P_red(c_sp+1,:);
        P_tail_cart = latlon2cart(P_tail,1);
        P_head_cart = latlon2cart(P_head,1);
        D           = P_length(c_sp);
        sp_P_cart   = P_tail_cart*sin((1-f)*D)/sin(D)+P_head_cart*sin(f*D)/sin(D);
        sp_P        = cart2latlon(sp_P_cart);
        sp_Q        = Q(r_sp,:);
        % Vertical node
    else
        % Select the minimum value for f, since this is the starting point
        f           = Lmin(r_sp,c_sp);
        sp_P        = P_red(c_sp,:);
        Q_tail      = Q(r_sp,:);
        Q_head      = Q(r_sp+1,:);
        Q_tail_cart = latlon2cart(Q_tail,1);
        Q_head_cart = latlon2cart(Q_head,1);
        D           = Q_length(r_sp);
        sp_Q_cart   = Q_tail_cart*sin((1-f)*D)/sin(D)+Q_head_cart*sin(f*D)/sin(D);
        sp_Q        = cart2latlon(sp_Q_cart);
    end
    
    %%%%%%%%%%%%%%%%%
    %%%  Endpoint %%%
    %%%%%%%%%%%%%%%%%
    
    % Corner node
    if Type_ep == 1
        ep_P = P_red(c_ep,:);
        ep_Q = Q(r_ep,:);
        % Horizontal node
    elseif Type_ep == 2
        % Select the maximum value for f, since this is the endpoint
        f           = Bmax(r_ep,c_ep);
        P_tail      = P_red(c_ep,:);
        P_head      = P_red(c_ep+1,:);
        P_tail_cart = latlon2cart(P_tail,1);
        P_head_cart = latlon2cart(P_head,1);
        D           = P_length(c_ep);
        sp_P_cart   = P_tail_cart*sin((1-f)*D)/sin(D)+P_head_cart*sin(f*D)/sin(D);
        ep_P        = cart2latlon(sp_P_cart);
        ep_Q        = Q(r_ep,:);
        % Vertical node
    else
        % Select the minimum value for f, since this is the starting point
        f           = Lmax(r_ep,c_ep);
        ep_P        = P_red(c_ep,:);
        Q_tail      = Q(r_ep,:);
        Q_head      = Q(r_ep+1,:);
        Q_tail_cart = latlon2cart(Q_tail,1);
        Q_head_cart = latlon2cart(Q_head,1);
        D           = Q_length(r_ep);
        ep_Q_cart   = Q_tail_cart*sin((1-f)*D)/sin(D)+Q_head_cart*sin(f*D)/sin(D);
        ep_Q        = cart2latlon(ep_Q_cart);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%% Computing overall length of the sequence for both trajectories %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% starting point: corner node, endpoint: corner node %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if Type_sp==1 && Type_ep==1
        % Length of sequence for P
        lseq_P = zeros(c_ep-c_sp,1);
        for ii=1:c_ep-c_sp
            lseq_P(ii) = Haversine(P_red(c_sp+ii-1,:),P_red(c_sp+ii,:),1);
        end
        lseq_P = sum(lseq_P);
        % Length of sequence for Q
        lseq_Q = zeros(r_ep-r_sp,1);
        for ii=1:r_ep-r_sp
            lseq_Q(ii) = Haversine(Q(r_sp+ii-1,:),Q(r_sp+ii,:),1);
        end
        lseq_Q = sum(lseq_Q);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% starting point: corner node, endpoint: horizontal node %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif Type_sp==1 && Type_ep==2
        % Length of sequence for P
        if c_sp == c_ep
            lseq_P = Haversine(sp_P,ep_P,1);
        else
            lseq1 = zeros(c_ep-c_sp,1);
            for ii=1:c_ep-c_sp
                lseq1(ii) = Haversine(P_red(c_sp+ii-1,:),P_red(c_sp+ii,:),1);
            end
            lseq1  = sum(lseq1);
            lseq2  = Haversine(P_red(c_ep,:),ep_P,1);
            lseq_P = lseq1+lseq2;
        end
        % Length of sequence for Q
        if r_sp == r_ep
            lseq_Q = 0;
        else
            lseq1 = zeros(r_ep-r_sp,1);
            for ii=1:r_ep-r_sp
                lseq1(ii) = Haversine(Q(r_sp+ii-1,:),Q(r_sp+ii,:),1);
            end
            lseq1  = sum(lseq1);
            lseq_Q = lseq1;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% starting point: corner node, endpoint: vertical node %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif Type_sp==1 && Type_ep==3
        % Length of sequence for P
        if c_sp == c_ep
            lseq_P = 0;
        else
            lseq1 = zeros(c_ep-c_sp,1);
            for ii=1:c_ep-c_sp
                lseq1(ii) = Haversine(P_red(c_sp+ii-1,:),P_red(c_sp+ii,:),1);
            end
            lseq1  = sum(lseq1);
            lseq_P = lseq1;
        end
        % Length of sequence for Q
        if r_sp == r_ep
            lseq_Q = Haversine(sp_Q,ep_Q,1);
        else
            lseq1 = zeros(r_ep-r_sp,1);
            for ii=1:r_ep-r_sp
                lseq1(ii) = Haversine(Q(r_sp+ii-1,:),Q(r_sp+ii,:),1);
            end
            lseq1  = sum(lseq1);
            lseq2  = Haversine(Q(r_ep,:),ep_Q,1);
            lseq_Q = lseq1+lseq2;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% starting point: horizontal node, endpoint: corner node %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif Type_sp==2 && Type_ep==1
        % Length of sequence for P
        if c_ep-c_sp == 1
            lseq_P = Haversine(sp_P,P_red(c_ep,:),1);
        else
            lseq1 = zeros(c_ep-c_sp-1,1);
            for ii=1:c_ep-c_sp-1
                lseq1(ii) = Haversine(P_red(c_sp+ii,:),P_red(c_sp+ii+1,:),1);
            end
            lseq1  = sum(lseq1);
            lseq2  = Haversine(sp_P,P_red(c_sp+1,:),1);
            lseq_P = lseq1+lseq2;
        end
        % Length of sequence for Q
        if r_sp == r_ep
            lseq_Q = 0;
        else
            lseq1 = zeros(r_ep-r_sp,1);
            for ii=1:r_ep-r_sp
                lseq1(ii) = Haversine(Q(r_sp+ii-1,:),Q(r_sp+ii,:),1);
            end
            lseq1  = sum(lseq1);
            lseq_Q = lseq1;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        %%% starting point: horizontal node, endpoint: horizontal node %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif Type_sp==2 && Type_ep==2
        % Length of sequence for P
        if c_ep-c_sp==1
            lseq1 = Haversine(sp_P,P_red(c_sp+1,:),1);
            lseq2 = Haversine(P_red(c_sp+1,:),ep_P,1);
        else
            lseq1 = zeros(c_ep-c_sp-1,1);
            for ii=1:c_ep-c_sp-1
                lseq1(ii) = Haversine(P_red(c_sp+ii,:),P_red(c_sp+ii+1,:),1);
            end
            lseq1  = sum(lseq1);
            lseq2  = Haversine(sp_P,P_red(c_sp+1,:),1);
            lseq3  = Haversine(P_red(c_ep,:),ep_Q,1);
            lseq_P = lseq1+lseq2+lseq3;
        end
        % Length of sequence for Q
        if r_ep == r_sp
            lseq_Q = 0;
        else
            lseq1 = zeros(r_ep-r_sp,1);
            for ii=1:r_ep-r_sp
                lseq1(ii) = Haversine(Q(r_sp+ii-1,:),Q(r_sp+ii,:),1);
            end
            lseq1  = sum(lseq1);
            lseq_Q = lseq1;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% starting point: horizontal node, endpoint: vertical node %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif Type_sp==2 && Type_ep==3
        % Length of sequence for P
        if c_ep-c_sp==1
            lseq_P = Haversine(sp_P,ep_P,1);
        else
            lseq1 = zeros(c_ep-c_sp-1,1);
            for ii=1:c_ep-c_sp-1
                lseq1(ii) = Haversine(P_red(c_sp+ii,:),P_red(c_sp+ii+1,:),1);
            end
            lseq1  = sum(lseq1);
            lseq2  = Haversine(sp_P,P_red(c_sp+1,:),1);
            lseq_P = lseq1+lseq2;
        end
        % Length of sequence for Q
        if r_ep == r_sp
            lseq_Q = Haversine(Q(r_sp,:),ep_Q,1);
        else
            lseq1 = zeros(r_ep-r_sp,1);
            for ii=1:r_ep-r_sp
                lseq1(ii) = Haversine(Q(r_sp+ii-1,:),Q(r_sp+ii,:),1);
            end
            lseq1  = sum(lseq1);
            lseq2  = Haversine(Q(r_ep,:),ep_Q,1);
            lseq_Q = lseq1+lseq2;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% starting point: vertical node, endpoint: corner node %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif Type_sp==3 && Type_ep==1
        % Length of sequence for P
        if c_ep == c_sp
            lseq_P = 0;
        else
            lseq1 = zeros(c_ep-c_sp,1);
            for ii=1:c_ep-c_sp
                lseq1(ii) = Haversine(P_red(c_sp+ii-1,:),P_red(c_sp+ii,:),1);
            end
            lseq1  = sum(lseq1);
            lseq_P = lseq1;
        end
        % Length of sequence for Q
        if r_ep - r_sp == 1
            lseq_Q = Haversine(sp_Q,Q(r_ep,:),1);
        else
            lseq1 = zeros(r_ep-r_sp-1,1);
            for ii=1:r_ep-r_sp-1
                lseq1(ii) = Haversine(Q(r_sp+ii,:),Q(r_sp+ii+1,:),1);
            end
            lseq1  = sum(lseq1);
            lseq2  = Haversine(sp_Q,Q(r_sp+1,:),1);
            lseq_Q = lseq1+lseq2;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% starting point: vertical node, endpoint: horizontal node %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif Type_sp==3 && Type_ep==2
        % Length of sequence for P
        if c_sp == c_ep
            lseq_P = Haversine(sp_P,ep_P,1);
        else
            lseq1 = zeros(c_ep-c_sp,1);
            for ii=1:c_ep-c_sp
                lseq1(ii) = Haversine(P_red(c_sp+ii-1,:),P_red(c_sp+ii,:),1);
            end
            lseq2  = Haversine(P_red(c_ep,:),ep_P,1);
            lseq_P = lseq1+lseq2;
        end
        % Length of sequence for Q
        if r_ep - r_sp == 1
            lseq_Q = Haversine(sp_Q,ep_Q,1);
        else
            lseq1 = zeros(r_ep-r_sp-1,1);
            for ii=1:r_ep-r_sp-1
                lseq1(ii) = Haversine(Q(r_sp+ii,:),Q(r_sp+ii+1,:),1);
            end
            lseq2  = Haversine(sp_Q,Q(r_sp+1),1);
            lseq_Q = lseq1+lseq2; 
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% starting point: vertical node, endpoint: vertical node %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        % Length of sequence for P
        if c_ep == c_sp
            lseq_P = 0;
        else
            lseq1 = zeros(c_ep-c_sp,1);
            for ii=1:c_ep-c_sp
                lseq1(ii) = Haversine(P_red(c_sp+ii-1,:),P_red(c_sp+ii,:),1);
            end
            lseq1  = sum(lseq1);
            lseq_P = lseq1;
        end
        % Length of sequence for Q
        if r_ep-r_sp==1
            lseq1 = Haversine(sp_Q,Q(r_sp+1,:),1);
            lseq2 = Haversine(Q(r_sp+1,:),ep_Q,1);
        else
            lseq1 = zeros(r_ep-r_sp-1,1);
            for ii=1:r_ep-r_sp-1
                lseq1(ii) = Haversine(Q(r_sp+ii,:),Q(r_sp+ii+1,:),1);
            end
            lseq1  = sum(lseq1);
            lseq2  = Haversine(sp_Q,Q(r_sp+1,:),1);
            lseq3  = Haversine(Q(r_ep,:),ep_Q,1);
            lseq_Q = lseq1+lseq2+lseq3;
        end
    end
    
    lseq_P_perc = lseq_P/sum(P_length);
    lseq_Q_perc = lseq_Q/sum(Q_length);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     r_vec = r_sp:r_ep;
%     c_vec = c_sp:c_ep;
%         
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%% Length of sequence for P %%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if numel(c_vec)==1
%         lseq_P = Haversine(sp_P,ep_P,1);
%     elseif numel(c_vec)==2
%         lseq1 = Haversine(sp_P,P(c_ep,:),1);
%         lseq2 = Haversine(P(c_ep,:),ep_P,1);
%         lseq_P  = lseq1+lseq2;
%     else
%         lseq1 = Haversine(sp_P,P(c_sp+1,:),1);
%         lseq2 = Haversine(P(c_ep,:),ep_P,1);
%         lseq3 = zeros(numel(c_vec)-2,1);
%         for ii=1:numel(c_vec)-2
%             lseq3(ii) = Haversine(P(c_sp+ii,:),P(c_sp+ii+1,:),1);
%         end
%         lseq_P = lseq1+lseq2+sum(lseq3);
%     end
%     
%     lseq_P_perc = lseq_P/sum(P_length);
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%% Length of sequence for Q %%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if numel(r_vec)==1
%         lseq_Q = Haversine(sp_Q,ep_Q,1);
%     elseif numel(r_vec)==2
%         lseq1 = Haversine(sp_Q,Q(r_ep,:),1);
%         lseq2 = Haversine(Q(r_ep,:),ep_Q,1);
%         lseq_Q  = lseq1+lseq2;
%     else
%         lseq1 = Haversine(sp_Q,Q(r_sp+1,:),1);
%         lseq2 = Haversine(Q(r_ep,:),ep_Q,1);
%         lseq3 = zeros(numel(r_vec)-2,1);
%         for ii=1:numel(r_vec)-2
%             lseq3(ii) = Haversine(Q(r_sp+ii,:),Q(r_sp+ii+1,:),1);
%         end
%         lseq_Q = lseq1+lseq2+sum(lseq3);
%     end
%     
%     lseq_Q_perc = lseq_Q/sum(Q_length);
% 
     Overlap_points(i,:) = horzcat(sp_P,sp_Q,ep_P,ep_Q,lseq_P,lseq_P_perc,lseq_Q,lseq_Q_perc);
    
end

time = toc;

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
plot(P_red(:,2),P_red(:,1),'Color','b','Linewidth',2,'Marker','*','Markersize',6)
plot(Q(:,2),Q(:,1),'Color','k','Linewidth',2,'Marker','*','Markersize',6)
plot(Overlap_points([1 3 4],2),Overlap_points([1 3 4],1),'Color','g','Linewidth',2.5,'Linestyle','none','Marker','x','Markersize',10)
plot(-2.473,0.092,'Color','g','Linewidth',2.5,'Linestyle','none','Marker','x','Markersize',10)
plot(Overlap_points(:,4),Overlap_points(:,3),'Color','g','Linewidth',2.5,'Linestyle','none','Marker','x','Markersize',10)
plot(Overlap_points(:,6),Overlap_points(:,5),'Color','r','Linewidth',2.5,'Linestyle','none','Marker','x','Markersize',10)
plot(Overlap_points(:,8),Overlap_points(:,7),'Color','r','Linewidth',2.5,'Linestyle','none','Marker','x','Markersize',10)
xlabel('Longitude [deg]','Fontname','Avantgarde','Fontsize',14)
ylabel('Latitude [deg]','Fontname','Avantgarde','Fontsize',14)
l1 = legend([h1,h2],'P polygonal curve','Q polygonal curve');
set(l1,'Fontname','Avantgarde','Fontsize',14)
grid on




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R        = 3;
% C        = 6;
% mu_R     = 1;
% mu_C     = 4;
% p        = 8;
% q        = 4;
% Rin      = max([1,R-mu_R]);
% Cin      = max([1,C-mu_C]);
% Rfin     = min([R+mu_R,q+1]);
% Cfin     = min([C+mu_C,p+1]);
% Rvec_in  = (Rin:R)';
% Cvec_in  = (Cin:C)';
% Rvec_fin = (R:Rfin)';
% Cvec_fin = (C:Cfin)';
% 
% %%%%%%%%%%%%%%%%%%%%%
% %%% Vertical node %%%
% %%%%%%%%%%%%%%%%%%%%%
% extra_starting = [];
% for i=1:numel(Rvec_in)
%     dummy1 = (Rvec_in(i)-1)*(3*p+2)+(2*Cin-1:2*C-1)';
%     dummy2 = (Rvec_in(i)-1)*(3*p+2)+2*p+1+Cvec_in;
%     extra_starting = [extra_starting;dummy1;dummy2];
% end
% extra_ending = [];
% for i=1:numel(Rvec_fin)
%     dummy1 = (Rvec_fin(i)-1)*(3*p+2)+2*p+1+Cvec_fin;
%     dummy2 = Rvec_fin(i)*(3*p+2)+(2*C-1:2*Cfin-1)';
%     extra_ending = [extra_ending;dummy1;dummy2];
% end
% 
% R        = 4;
% C        = 5;
% mu_R     = 2;
% mu_C     = 2;
% p        = 8;
% q        = 4;
% Rin      = max([1,R-mu_R]);
% Cin      = max([1,C-mu_C]);
% Rfin     = min([R+mu_R,q+1]);
% Cfin     = min([C+mu_C,p+1]);
% Rvec_in  = (Rin:R)';
% Cvec_in  = (Cin:C)';
% Rvec_fin = (R:Rfin)';
% Cvec_fin = (C:Cfin)';
% 
% %%%%%%%%%%%%%%%%%%%%%%%
% %%% Horizontal node %%%
% %%%%%%%%%%%%%%%%%%%%%%%
% extra_starting = [];
% for i=1:numel(Rvec_in)
%     if Rvec_in(i)~=R
%         dummy1 = (Rvec_in(i)-1)*(3*p+2)+(2*Cin-1:2*C)';
%         dummy2 = (Rvec_in(i)-1)*(3*p+2)+2*p+1+Cvec_in;
%     else
%         dummy1 = (Rvec_in(i)-1)*(3*p+2)+(2*Cin-1:2*C)';
%         dummy2 = [];
%     end
%     extra_starting = [extra_starting;dummy1;dummy2];
% end
% extra_ending = [];
% for i=1:numel(Rvec_fin)
%     if Rvec_fin(i)~=Rfin
%         dummy1 = (Rvec_fin(i)-1)*(3*p+2)+(2*C:2*Cfin+1)';
%         dummy2 = (Rvec_fin(i)-1)*(3*p+2)+2*p+2+(C:Cfin)';
%     else
%         dummy1 = (Rvec_fin(i)-1)*(3*p+2)+(2*C:2*Cfin+1)';
%         dummy2 = [];
%     end
%     extra_ending = [extra_ending;dummy1;dummy2];
% end
% 
% R        = 3;
% C        = 5;
% mu_R     = 1;
% mu_C     = 3;
% p        = 8;
% q        = 4;
% Rin      = max([1,R-mu_R]);
% Cin      = max([1,C-mu_C]);
% Rfin     = min([R+mu_R,q+1]);
% Cfin     = min([C+mu_C,p+1]);
% Rvec_in  = (Rin:R)';
% Cvec_in  = (Cin:C)';
% Rvec_fin = (R:Rfin)';
% Cvec_fin = (C:Cfin)';
% 
% %%%%%%%%%%%%%%%%%%%
% %%% Corner node %%%
% %%%%%%%%%%%%%%%%%%%
% extra_starting = [];
% for i=1:numel(Rvec_in)
%     if Rvec_in(i)~=R
%         dummy1 = (Rvec_in(i)-1)*(3*p+2)+(2*Cin-1:2*C-1)';
%         dummy2 = (Rvec_in(i)-1)*(3*p+2)+2*p+1+Cvec_in;
%     else
%         dummy1 = (Rvec_in(i)-1)*(3*p+2)+(2*Cin-1:2*C-1)';
%         dummy2 = [];
%     end
%     extra_starting = [extra_starting;dummy1;dummy2];
% end
% extra_ending = [];
% for i=1:numel(Rvec_fin)
%     if Rvec_fin(i)~=Rfin
%         dummy1 = (Rvec_fin(i)-1)*(3*p+2)+(2*C-1:2*Cfin-1)';
%         dummy2 = (Rvec_fin(i)-1)*(3*p+2)+2*p+1+(C:Cfin)';
%     else
%         dummy1 = (Rvec_fin(i)-1)*(3*p+2)+(2*C-1:2*Cfin-1)';
%         dummy2 = [];
%     end
%     extra_ending = [extra_ending;dummy1;dummy2];
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% tic;
% [is_solution,A]                       = DFS(1,508,A_list_eps);
% time=toc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%startpoint = 1;
%endpoint   = 9;
%G          = {[2 4 5];[3 5 6];6;[5 7 8];[6 8 9];9;8;9;NaN};
%[is_solution,A] = DFS(startpoint,endpoint,G);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of points interior to the great circle segment between the tail
% and the head nodes
% N_points = 10;
% 
% oversampled_flight_plan = [];
% 
% for i=1:length(Q(:,1))-1
%     tail       = Q(i,:);
%     head       = Q(i+1,:);
%     new_points = gcwaypts(tail(1),tail(2),head(1),head(2),N_points+2);
%      
%      if i~=length(Q(:,1))-1
%          oversampled_flight_plan = vertcat(oversampled_flight_plan,new_points(1:end-1,:));
%      else
%          oversampled_flight_plan = vertcat(oversampled_flight_plan,new_points);
%      end
%     
% end
% 
% disp('Begin case 3')
% tic
% [cm, cSq, cD] = func_discrete_Frechet_lat_lon_NM_v02(P,oversampled_flight_plan);
% time3 = toc;
% 
% figure()
% hold on
% plot(P(:,2),P(:,1),'Color','b','Linewidth',2,'Marker','*','Markersize',8)
% plot(oversampled_flight_plan(:,2),oversampled_flight_plan(:,1),'Color','r','Linewidth',2,'Marker','*','Markersize',8)
% grid on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure()
% hold on
% plot(flight1(:,5),flight1(:,4),'Color','b','Linestyle','-','Linewidth',2)
% plot(flight2(:,5),flight2(:,4),'Color','r','Linestyle','-','Linewidth',2)
% grid on
% axis equal
% 
% flight_plan = [47.4 -122.3;46 -122.2;40.9 -120.9;39.8 -120.85;...
%                37 -120.3; 35 -119.4; 34 -118.4];
% 
% tic
% [FD3,no_sol3,sol3] = FD_computation_v02(flight1(:,4:5),flight_plan,Re,eps,tolerance);
% time4 = toc; 
% 
% tic
% [cm2, cSq2, cD2] = func_discrete_Frechet_lat_lon_NM_v02(flight1(:,4:5),flight_plan);
% time5 = toc;
% 
% % Number of points interior to the great circle segment between the tail
% % and the head nodes
% N_points = 15;
% 
% oversampled_flight_plan = [];
% 
% for i=1:length(flight_plan(:,1))-1
%     tail       = flight_plan(i,:);
%     head       = flight_plan(i+1,:);
%     new_points = gcwaypts(tail(1),tail(2),head(1),head(2),N_points+2);
%      
%      if i~=length(flight_plan(:,1))-1
%          oversampled_flight_plan = vertcat(oversampled_flight_plan,new_points(1:end-1,:));
%      else
%          oversampled_flight_plan = vertcat(oversampled_flight_plan,new_points);
%      end
%     
% end
% 
% tic
% [cm3, cSq3, cD3] = func_discrete_Frechet_lat_lon_NM_v02(flight1(:,4:5),oversampled_flight_plan);
% time6 = toc;

% 
% max_delta_heading   = 2;
% delta_heading1      = abs(heading1(2:end)-heading1(1:end-1));
% dummy1              = delta_heading1<=max_delta_heading;
% idx1                = find(dummy1==1)+1;
% flight1_red         = flight1;
% flight1_red(idx1,:) = [];

% startpoint  = 1;
% endpoint    = 9;
% G           = cell(9,1);
% G{1,1}      = [2 4 5];
% G{2,1}      = [3 5 6];
% G{3,1}      = 6;
% G{4,1}      = [5 7 8];
% G{5,1}      = [6 8 9];
% G{6,1}      = 9;
% G{7,1}      = 8;
% G{8,1}      = 9;
% G{9,1}      = [];
% is_solution = DFS(startpoint,endpoint,G);

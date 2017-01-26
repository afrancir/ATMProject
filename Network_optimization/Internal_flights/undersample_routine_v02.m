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

alpha = 1/3;
[Lmin_eps,Lmax_eps,Bmin_eps,Bmax_eps,A_list_eps,Node_ID_eps] = free_space_shape_v02(P_red,Q,Re,alpha*FD);
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

Overlap_points = zeros(length(Merged_Sequences(:,1)),8);
N_f            = 100;

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
        
        
        % NEW PART GOES INSIDE HERE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if Lmax(r_sp,c_sp)==1 && Bmax(r_sp+1,c_sp-1)==1
            disp(['Starting node is vertical. Row is ',num2str(r_sp),', column is ',num2str(c_sp)])
            %f_vec_P     = linspace(Bmin(r_sp+1,c_sp-1),Bmax(r_sp+1,c_sp-1),N_f);
            %f_vec_Q     = linspace(Lmin(r_sp,c_sp),Lmax(r_sp,c_sp),N_f);
            f_vec_P     = linspace(0.05,1,N_f);
            f_vec_Q     = linspace(0.05,1,N_f);
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
            dist2       = zeros(size(dist_mat));
            for ii=1:numel(f_vec_Q)
                Q_cart   = Q_tail_cart*sin((1-f_vec_Q(ii))*DQ)/sin(DQ)+Q_head_cart*sin(f_vec_Q(ii)*DQ)/sin(DQ);
                for jj=1:numel(f_vec_P)
                    P_cart   = P_tail_cart*sin((1-f_vec_P(jj))*DP)/sin(DP)+P_head_cart*sin(f_vec_P(jj)*DP)/sin(DP);
                    P_ll     = cart2latlon(P_cart);
                    Q_ll     = cart2latlon(Q_cart);
                    dist_mat(ii,jj) = Haversine(P_ll,Q_ll,1);
                    dist2(ii,jj)    = Haversine(P_ll,P_head,1)+Haversine(Q_ll,Q_head,1);
                end
            end   
            [r,c]    = find(dist_mat<=alpha*FD);
            dist_vec = zeros(numel(r));
            for ii=1:numel(r)
                dist_vec(ii) = dist2(r(ii),c(ii));
            end
            idx       = find(max(dist_vec));
            best_f_Q  = f_vec_Q(r(idx));
            best_f_P  = f_vec_P(c(idx));
            P_cart_sp = P_tail_cart*sin((1-best_f_P)*DP)/sin(DP)+P_head_cart*sin(best_f_P*DP)/sin(DP);
            Q_cart_sp = Q_tail_cart*sin((1-best_f_Q)*DQ)/sin(DQ)+Q_head_cart*sin(best_f_Q*DQ)/sin(DQ);
            sp_P      = cart2latlon(P_cart_sp)
            sp_Q      = cart2latlon(Q_cart_sp)
        else
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%         f           = Lmin(r_sp,c_sp);
%         sp_P        = P_red(c_sp,:);
%         Q_tail      = Q(r_sp,:);
%         Q_head      = Q(r_sp+1,:);
%         Q_tail_cart = latlon2cart(Q_tail,1);
%         Q_head_cart = latlon2cart(Q_head,1);
%         D           = Q_length(r_sp);
%         sp_Q_cart   = Q_tail_cart*sin((1-f)*D)/sin(D)+Q_head_cart*sin(f*D)/sin(D);
%         sp_Q        = cart2latlon(sp_Q_cart);
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % NEW PART GOES INSIDE HERE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if Lmin(r_ep,c_ep)==0 && Bmin(r_ep,c_ep)==0
            disp(['Ending node is vertical. Row is ',num2str(r_ep),', column is ',num2str(c_ep)])
            f_vec_P     = linspace(0,1,N_f);
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
            dist2       = zeros(numel(f_vec_Q),numel(f_vec_P));
            for ii=1:numel(f_vec_Q)
                Q_cart   = Q_tail_cart*sin((1-f_vec_Q(ii))*DQ)/sin(DQ)+Q_head_cart*sin(f_vec_Q(ii)*DQ)/sin(DQ);
                Q_ll     = cart2latlon(Q_cart);
                for jj=1:numel(f_vec_P)
                    P_cart          = P_tail_cart*sin((1-f_vec_P(jj))*DP)/sin(DP)+P_head_cart*sin(f_vec_P(jj)*DP)/sin(DP);
                    P_ll            = cart2latlon(P_cart);
                    dist_mat(ii,jj) = Haversine(P_ll,Q_ll,1);
                    dist2(ii,jj)    = Haversine(P_tail,P_ll,1)+Haversine(Q_tail,Q_ll,1);
                end
            end
            [r,c]    = find(dist_mat<=alpha*FD);
            dist_vec = zeros(numel(r),1);
            for ii=1:numel(r)
                dist_vec(ii) = dist2(r(ii),c(ii));
            end
            idx       = find(max(dist_vec));
            best_f_Q  = f_vec_Q(r(idx));
            best_f_P  = f_vec_P(c(idx));
            P_cart_sp = P_tail_cart*sin((1-best_f_P)*DP)/sin(DP)+P_head_cart*sin(best_f_P*DP)/sin(DP);
            Q_cart_sp = Q_tail_cart*sin((1-best_f_Q)*DQ)/sin(DQ)+Q_head_cart*sin(best_f_Q*DQ)/sin(DQ);
            ep_P      = cart2latlon(P_cart_sp)
            ep_Q      = cart2latlon(Q_cart_sp)
        else
            f           = Lmax(r_sp,c_sp);
            ep_P        = P_red(c_sp,:);
            Q_tail      = Q(r_sp,:);
            Q_head      = Q(r_sp+1,:);
            Q_tail_cart = latlon2cart(Q_tail,1);
            Q_head_cart = latlon2cart(Q_head,1);
            D           = Q_length(r_sp);
            ep_Q_cart   = Q_tail_cart*sin((1-f)*D)/sin(D)+Q_head_cart*sin(f*D)/sin(D);
            ep_Q        = cart2latlon(ep_Q_cart);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Select the minimum value for f, since this is the starting point
        %f           = Lmax(r_ep,c_ep);
        %ep_P        = P_red(c_ep,:);
        %Q_tail      = Q(r_ep,:);
        %Q_head      = Q(r_ep+1,:);
        %Q_tail_cart = latlon2cart(Q_tail,1);
        %Q_head_cart = latlon2cart(Q_head,1);
        %D           = Q_length(r_ep);
        %ep_Q_cart   = Q_tail_cart*sin((1-f)*D)/sin(D)+Q_head_cart*sin(f*D)/sin(D);
        %ep_Q        = cart2latlon(ep_Q_cart);
    end
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%     %%% Computing overall length of the sequence for both trajectories %%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%% starting point: corner node, endpoint: corner node %%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if Type_sp==1 && Type_ep==1
%         % Length of sequence for P
%         lseq_P = zeros(c_ep-c_sp,1);
%         for ii=1:c_ep-c_sp
%             lseq_P(ii) = Haversine(P_red(c_sp+ii-1,:),P_red(c_sp+ii,:),1);
%         end
%         lseq_P = sum(lseq_P);
%         % Length of sequence for Q
%         lseq_Q = zeros(r_ep-r_sp,1);
%         for ii=1:r_ep-r_sp
%             lseq_Q(ii) = Haversine(Q(r_sp+ii-1,:),Q(r_sp+ii,:),1);
%         end
%         lseq_Q = sum(lseq_Q);
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%% starting point: corner node, endpoint: horizontal node %%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     elseif Type_sp==1 && Type_ep==2
%         % Length of sequence for P
%         if c_sp == c_ep
%             lseq_P = Haversine(sp_P,ep_P,1);
%         else
%             lseq1 = zeros(c_ep-c_sp,1);
%             for ii=1:c_ep-c_sp
%                 lseq1(ii) = Haversine(P_red(c_sp+ii-1,:),P_red(c_sp+ii,:),1);
%             end
%             lseq1  = sum(lseq1);
%             lseq2  = Haversine(P_red(c_ep,:),ep_P,1);
%             lseq_P = lseq1+lseq2;
%         end
%         % Length of sequence for Q
%         if r_sp == r_ep
%             lseq_Q = 0;
%         else
%             lseq1 = zeros(r_ep-r_sp,1);
%             for ii=1:r_ep-r_sp
%                 lseq1(ii) = Haversine(Q(r_sp+ii-1,:),Q(r_sp+ii,:),1);
%             end
%             lseq1  = sum(lseq1);
%             lseq_Q = lseq1;
%         end
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%% starting point: corner node, endpoint: vertical node %%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     elseif Type_sp==1 && Type_ep==3
%         % Length of sequence for P
%         if c_sp == c_ep
%             lseq_P = 0;
%         else
%             lseq1 = zeros(c_ep-c_sp,1);
%             for ii=1:c_ep-c_sp
%                 lseq1(ii) = Haversine(P_red(c_sp+ii-1,:),P_red(c_sp+ii,:),1);
%             end
%             lseq1  = sum(lseq1);
%             lseq_P = lseq1;
%         end
%         % Length of sequence for Q
%         if r_sp == r_ep
%             lseq_Q = Haversine(sp_Q,ep_Q,1);
%         else
%             lseq1 = zeros(r_ep-r_sp,1);
%             for ii=1:r_ep-r_sp
%                 lseq1(ii) = Haversine(Q(r_sp+ii-1,:),Q(r_sp+ii,:),1);
%             end
%             lseq1  = sum(lseq1);
%             lseq2  = Haversine(Q(r_ep,:),ep_Q,1);
%             lseq_Q = lseq1+lseq2;
%         end
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%% starting point: horizontal node, endpoint: corner node %%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     elseif Type_sp==2 && Type_ep==1
%         % Length of sequence for P
%         if c_ep-c_sp == 1
%             lseq_P = Haversine(sp_P,P_red(c_ep,:),1);
%         else
%             lseq1 = zeros(c_ep-c_sp-1,1);
%             for ii=1:c_ep-c_sp-1
%                 lseq1(ii) = Haversine(P_red(c_sp+ii,:),P_red(c_sp+ii+1,:),1);
%             end
%             lseq1  = sum(lseq1);
%             lseq2  = Haversine(sp_P,P_red(c_sp+1,:),1);
%             lseq_P = lseq1+lseq2;
%         end
%         % Length of sequence for Q
%         if r_sp == r_ep
%             lseq_Q = 0;
%         else
%             lseq1 = zeros(r_ep-r_sp,1);
%             for ii=1:r_ep-r_sp
%                 lseq1(ii) = Haversine(Q(r_sp+ii-1,:),Q(r_sp+ii,:),1);
%             end
%             lseq1  = sum(lseq1);
%             lseq_Q = lseq1;
%         end
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%         %%% starting point: horizontal node, endpoint: horizontal node %%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     elseif Type_sp==2 && Type_ep==2
%         % Length of sequence for P
%         if c_ep-c_sp==1
%             lseq1 = Haversine(sp_P,P_red(c_sp+1,:),1);
%             lseq2 = Haversine(P_red(c_sp+1,:),ep_P,1);
%         else
%             lseq1 = zeros(c_ep-c_sp-1,1);
%             for ii=1:c_ep-c_sp-1
%                 lseq1(ii) = Haversine(P_red(c_sp+ii,:),P_red(c_sp+ii+1,:),1);
%             end
%             lseq1  = sum(lseq1);
%             lseq2  = Haversine(sp_P,P_red(c_sp+1,:),1);
%             lseq3  = Haversine(P_red(c_ep,:),ep_Q,1);
%             lseq_P = lseq1+lseq2+lseq3;
%         end
%         % Length of sequence for Q
%         if r_ep == r_sp
%             lseq_Q = 0;
%         else
%             lseq1 = zeros(r_ep-r_sp,1);
%             for ii=1:r_ep-r_sp
%                 lseq1(ii) = Haversine(Q(r_sp+ii-1,:),Q(r_sp+ii,:),1);
%             end
%             lseq1  = sum(lseq1);
%             lseq_Q = lseq1;
%         end
%         
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%% starting point: horizontal node, endpoint: vertical node %%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     elseif Type_sp==2 && Type_ep==3
%         % Length of sequence for P
%         if c_ep-c_sp==1
%             lseq_P = Haversine(sp_P,ep_P,1);
%         else
%             lseq1 = zeros(c_ep-c_sp-1,1);
%             for ii=1:c_ep-c_sp-1
%                 lseq1(ii) = Haversine(P_red(c_sp+ii,:),P_red(c_sp+ii+1,:),1);
%             end
%             lseq1  = sum(lseq1);
%             lseq2  = Haversine(sp_P,P_red(c_sp+1,:),1);
%             lseq_P = lseq1+lseq2;
%         end
%         % Length of sequence for Q
%         if r_ep == r_sp
%             lseq_Q = Haversine(Q(r_sp,:),ep_Q,1);
%         else
%             lseq1 = zeros(r_ep-r_sp,1);
%             for ii=1:r_ep-r_sp
%                 lseq1(ii) = Haversine(Q(r_sp+ii-1,:),Q(r_sp+ii,:),1);
%             end
%             lseq1  = sum(lseq1);
%             lseq2  = Haversine(Q(r_ep,:),ep_Q,1);
%             lseq_Q = lseq1+lseq2;
%         end
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%% starting point: vertical node, endpoint: corner node %%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     elseif Type_sp==3 && Type_ep==1
%         % Length of sequence for P
%         if c_ep == c_sp
%             lseq_P = 0;
%         else
%             lseq1 = zeros(c_ep-c_sp,1);
%             for ii=1:c_ep-c_sp
%                 lseq1(ii) = Haversine(P_red(c_sp+ii-1,:),P_red(c_sp+ii,:),1);
%             end
%             lseq1  = sum(lseq1);
%             lseq_P = lseq1;
%         end
%         % Length of sequence for Q
%         if r_ep - r_sp == 1
%             lseq_Q = Haversine(sp_Q,Q(r_ep,:),1);
%         else
%             lseq1 = zeros(r_ep-r_sp-1,1);
%             for ii=1:r_ep-r_sp-1
%                 lseq1(ii) = Haversine(Q(r_sp+ii,:),Q(r_sp+ii+1,:),1);
%             end
%             lseq1  = sum(lseq1);
%             lseq2  = Haversine(sp_Q,Q(r_sp+1,:),1);
%             lseq_Q = lseq1+lseq2;
%         end
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%% starting point: vertical node, endpoint: horizontal node %%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     elseif Type_sp==3 && Type_ep==2
%         % Length of sequence for P
%         if c_sp == c_ep
%             lseq_P = Haversine(sp_P,ep_P,1);
%         else
%             lseq1 = zeros(c_ep-c_sp,1);
%             for ii=1:c_ep-c_sp
%                 lseq1(ii) = Haversine(P_red(c_sp+ii-1,:),P_red(c_sp+ii,:),1);
%             end
%             lseq2  = Haversine(P_red(c_ep,:),ep_P,1);
%             lseq_P = lseq1+lseq2;
%         end
%         % Length of sequence for Q
%         if r_ep - r_sp == 1
%             lseq_Q = Haversine(sp_Q,ep_Q,1);
%         else
%             lseq1 = zeros(r_ep-r_sp-1,1);
%             for ii=1:r_ep-r_sp-1
%                 lseq1(ii) = Haversine(Q(r_sp+ii,:),Q(r_sp+ii+1,:),1);
%             end
%             lseq2  = Haversine(sp_Q,Q(r_sp+1),1);
%             lseq_Q = lseq1+lseq2; 
%         end
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%% starting point: vertical node, endpoint: vertical node %%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     else
%         % Length of sequence for P
%         if c_ep == c_sp
%             lseq_P = 0;
%         else
%             lseq1 = zeros(c_ep-c_sp,1);
%             for ii=1:c_ep-c_sp
%                 lseq1(ii) = Haversine(P_red(c_sp+ii-1,:),P_red(c_sp+ii,:),1);
%             end
%             lseq1  = sum(lseq1);
%             lseq_P = lseq1;
%         end
%         % Length of sequence for Q
%         if r_ep-r_sp==1
%             lseq1 = Haversine(sp_Q,Q(r_sp+1,:),1);
%             lseq2 = Haversine(Q(r_sp+1,:),ep_Q,1);
%         else
%             lseq1 = zeros(r_ep-r_sp-1,1);
%             for ii=1:r_ep-r_sp-1
%                 lseq1(ii) = Haversine(Q(r_sp+ii,:),Q(r_sp+ii+1,:),1);
%             end
%             lseq1  = sum(lseq1);
%             lseq2  = Haversine(sp_Q,Q(r_sp+1,:),1);
%             lseq3  = Haversine(Q(r_ep,:),ep_Q,1);
%             lseq_Q = lseq1+lseq2+lseq3;
%         end
%     end
%     
%     lseq_P_perc = lseq_P/sum(P_length);
%     lseq_Q_perc = lseq_Q/sum(Q_length);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

     Overlap_points(i,:) = horzcat(sp_P,sp_Q,ep_P,ep_Q);
    
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
plot(Overlap_points(:,2),Overlap_points(:,1),'Color','g','Linewidth',2.5,'Linestyle','none','Marker','x','Markersize',10)
plot(Overlap_points(:,4),Overlap_points(:,3),'Color','g','Linewidth',2.5,'Linestyle','none','Marker','x','Markersize',10)
plot(Overlap_points(:,6),Overlap_points(:,5),'Color','r','Linewidth',2.5,'Linestyle','none','Marker','x','Markersize',10)
plot(Overlap_points(:,8),Overlap_points(:,7),'Color','r','Linewidth',2.5,'Linestyle','none','Marker','x','Markersize',10)
xlabel('Longitude [deg]','Fontname','Avantgarde','Fontsize',14)
ylabel('Latitude [deg]','Fontname','Avantgarde','Fontsize',14)
l1 = legend([h1,h2],'P polygonal curve','Q polygonal curve');
set(l1,'Fontname','Avantgarde','Fontsize',14)
grid on






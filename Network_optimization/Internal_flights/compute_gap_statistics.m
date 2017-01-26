clc
clear all
close all

origin_ID      = 33;
destination_ID = 10;

% Loading Frechet Distance Matrix
FM = load(strcat(pwd,'/NETWORK/INT_FDMatrix/',num2str(origin_ID),...
    '/',num2str(origin_ID),'_',num2str(destination_ID),'_red.txt'));

%FM = [0 10 0.5 1 0.7;0 0 8 9 9;0 0 0 1 0.8;0 0 0 0 0.6;0 0 0 0 0];
%
%FM = [0 1 0.95 1 0.95;0 0 1 1 1;0 0 0 1 0.95;0 0 0 0 0.95;0 0 0 0 0];
%
%FM = [0 1 1 1 1;0 0 sqrt(2) 2 sqrt(2);0 0 0 sqrt(2) 2;0 0 0 0 sqrt(2);0 0 0 0 0];

N_clust = 1:6;
W       = zeros(numel(N_clust),1);

for i=1:numel(N_clust)
    
    % Number of clusters
    N_c     = N_clust(i);
    values  = zeros(N_c,2);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Compute real dispersion %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % Single cluster case
    if N_c == 1
        n_r  = numel(FM(:,1));
        D_r  = sum(FM(FM~=0));
        W(i) = D_r/(2*n_r);
        % More than 1 cluster
    else
        W_temp = zeros(N_c,1);
        N_traj = numel(FM(:,1));
        dist_vec_km = [];
        for ii=1:N_traj-1
            dist_vec_km= horzcat(dist_vec_km,FM(ii,ii+1:end));
        end
        dendrogram = linkage(dist_vec_km,'average');
        T          = cluster(dendrogram,'maxclust',N_c);
        for ii=1:N_c
            % find rows/columns associated with the current cluster
            r_c  = find(T==ii);
            n_r  = numel(r_c);
            % if cluster has only one element, set intra-cluster distance to
            % zero
            if numel(r_c) == 1
                D_r  = 0;
                W_temp(ii) = D_r/(2*n_r);
            else
                % isolate sub-matrix associated with current cluster
                FM_local = FM(r_c,r_c);
                % compute average intra-cluster distance
                D_r   = mean(FM_local(FM_local~=0));
                W_temp(ii) = D_r/(2*n_r);
            end
        end
        W(i) = sum(W_temp);
    end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute reference dispersion %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B     = 2000;
ref_W = zeros(numel(N_clust),B);


for i=1:numel(N_clust)
    
    % Number of clusters
    N_c     = N_clust(i);
    
    for j=1:B
        unif_dist = unifrnd(min(FM(FM~=0)),max(FM(FM~=0)),numel(FM(:,1))*(numel(FM(:,1))-1)/2,1);
        FM_B      = triu(squareform(unif_dist));
        if N_c == 1
            n_r        = numel(FM_B(:,1));
            D_r        = sum(FM_B(FM_B~=0));
            ref_W(i,j) = D_r/(2*n_r);
            % More than 1 cluster
        else
            W_temp = zeros(N_c,1);
            N_traj = numel(FM_B(:,1));
            dist_vec_km = [];
            for ii=1:N_traj-1
                dist_vec_km= horzcat(dist_vec_km,FM_B(ii,ii+1:end));
            end
            dendrogram = linkage(dist_vec_km,'average');
            T          = cluster(dendrogram,'maxclust',N_c);
            for ii=1:N_c
                % find rows/columns associated with the current cluster
                r_c  = find(T==ii);
                n_r  = numel(r_c);
                % if cluster has only one element, set intra-cluster distance to
                % zero
                if numel(r_c) == 1
                    D_r          = 0;
                    W_temp(ii,1) = D_r/(2*n_r);
                else
                    % isolate sub-matrix associated with current cluster
                    FM_local = FM_B(r_c,r_c);
                    % compute average intra-cluster distance
                    D_r          = mean(FM_local(FM_local~=0));
                    W_temp(ii,1) = D_r/(2*n_r);
                end
            end
            ref_W(i,j) = (sum(W_temp));
        end
    end
    
end

W     = log10(W);
ref_W = log10(ref_W);

GAP = zeros(numel(N_clust),1);
sk  = zeros(numel(N_clust),1);

for i=1:numel(N_clust)
    GAP(i) = 1/B*sum(ref_W(i,:))-W(i);
    dummy  = 1/B*sum(ref_W(i,:));
    sdk    = sqrt(1/B*sum((ref_W(i,:)-dummy).^2));
    sk(i)  = sdk*sqrt(1+1/B);
end

index = zeros(numel(N_clust)-1,1);

for i=1:numel(N_clust)-1
    index(i) = GAP(i)-GAP(i+1)+sk(i+1);
end


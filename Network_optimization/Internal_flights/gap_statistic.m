% Written by Alessandro Bombelli, June 09, 2016
% Function that implements tha gap statistic method to estimate the optimal
% number of clusters, given a (N x p) matrix of observations

function opt_num_clust = gap_statistic(X,B,N_clust)

[~,~,V] = svd(X);
Xp      = X*V;
Zp(:,1) = unifrnd(min(Xp(:,1)),max(Xp(:,1)),numel(Xp(:,1))*100,1);
Zp(:,2) = unifrnd(min(Xp(:,2)),max(Xp(:,2)),numel(Xp(:,2))*100,1);
Z       = Zp*V';

alpha = 30;
figure()
hold on
plot(Z(:,1),Z(:,2),'Color','r','Marker','o','Linestyle','none')
plot(X(:,1),X(:,2),'Color','b','Marker','x','Linestyle','none')
quiver(0,0,alpha*V(1,1),alpha*V(2,1))
quiver(0,0,alpha*V(1,2),alpha*V(2,2))
axis equal
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute Distance Matrix %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FM = zeros(numel(X(:,1)),numel(X(:,1)));

for i=1:numel(X(:,1))-1
    for j=i+1:numel(X(:,1))
        FM(i,j) = sqrt((X(i,1)-X(j,1))^2+((X(i,2)-X(j,2))^2));
    end
end

dist_vec = squareform(FM+FM');
W        = zeros(numel(N_clust),1);

for i=1:numel(N_clust)
    
    % Number of clusters
    N_c     = N_clust(i);
    
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
        W_temp     = zeros(N_c,1);
        dendrogram = linkage(dist_vec,'average');
        T          = cluster(dendrogram,'maxclust',N_c);
        
        for ii=1:N_c
            % find rows/columns associated with the current cluster
            r_c  = find(T==ii);
            n_r  = numel(r_c);
            % if cluster has only one element, set intra-cluster distance to
            % zero
            if numel(r_c) == 1
                D_r        = 0;
                W_temp(ii) = 0;
            else
                % isolate sub-matrix associated with current cluster
                FM_local = FM(r_c,r_c);
                % compute average intra-cluster distance
                D_r        = sum(FM_local(FM_local~=0));
                W_temp(ii) = D_r/(2*n_r);
            end
        end
        W(i) = sum(W_temp);
    end
end

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Computing reference dispersions')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute reference dispersion %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ref_W = zeros(numel(N_clust),B);

for i=1:numel(N_clust)
    
    % Number of clusters
    N_c     = N_clust(i);
    
    for j=1:B
        idx_perm = randperm(numel(Z(:,1)));
        X_B      = Z(idx_perm(1:numel(X(:,1))),:);
        FM_B     = zeros(numel(X_B(:,1)),numel(X_B(:,1)));
        
        for ii=1:numel(X_B(:,1))-1
            for jj=ii+1:numel(X_B(:,1))
                FM_B(ii,jj) = sqrt((X_B(ii,1)-X_B(jj,1))^2+((X_B(ii,2)-X_B(jj,2))^2));
            end
        end
        
        dist_vec = squareform(FM_B+FM_B');
        if N_c == 1
            n_r        = numel(FM_B(:,1));
            D_r        = sum(FM_B(FM_B~=0));
            ref_W(i,j) = D_r/(2*n_r);
            % More than 1 cluster
        else
            W_temp = zeros(N_c,1);
            dendrogram = linkage(dist_vec,'average');
            T          = cluster(dendrogram,'maxclust',N_c);
            
            
            for ii=1:N_c
                % find rows/columns associated with the current cluster
                r_c  = find(T==ii);
                n_r  = numel(r_c);
                % if cluster has only one element, set intra-cluster distance to
                % zero
                if numel(r_c) == 1
                    W_temp(ii,1) = 0;
                else
                    % isolate sub-matrix associated with current cluster
                    FM_local = FM_B(r_c,r_c);
                    % compute average intra-cluster distance
                    D_r          = sum(FM_local(FM_local~=0));
                    W_temp(ii,1) = D_r/(2*n_r);
                end
            end
            ref_W(i,j) = (sum(W_temp));
        end
    end
    
end

Wlog     = log10(W);
ref_Wlog = log10(ref_W);
Wlog_exp = zeros(size(Wlog));

GAP = zeros(numel(N_clust),1);
sk  = zeros(numel(N_clust),1);

for i=1:numel(N_clust)
    GAP(i)      = 1/B*sum(ref_Wlog(i,:))-Wlog(i);
    dummy       = 1/B*sum(ref_Wlog(i,:));
    Wlog_exp(i) = dummy;
    sdk         = sqrt(1/B*sum((ref_Wlog(i,:)-dummy).^2));
    sk(i)       = sdk*sqrt(1+1/B);
end

index = zeros(numel(N_clust)-1,1);

for i=1:numel(N_clust)-1
    index(i) = GAP(i)-GAP(i+1)+sk(i+1);
end

opt_num_clust = 1;

return
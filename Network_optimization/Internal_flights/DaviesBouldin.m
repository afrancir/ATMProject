function DB = DaviesBouldin(FM,N_clust,T,flag)

intra_cluster_distance = zeros(N_clust,1);
extra_cluster_distance = zeros(N_clust,N_clust);
R                      = zeros(N_clust,N_clust);
D                      = zeros(N_clust,1);
FM                     = FM+FM';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CASE 1: we compute the mean and use as a sort of centroid %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag
    disp('Considering mean')
    % Compute all intra-cluster distances
    for ii=1:N_clust
        % find rows/columns associated with the current cluster
        r_c  = find(T==ii);
        % if cluster has only one element, set intra-cluster distance to
        % zero
        if numel(r_c) == 1
            intra_cluster_distance(ii) = 0;
        else
            % isolate sub-matrix associated with current cluster
            FM_local = FM(r_c,r_c);
            % compute average intra-cluster distance
            if isempty(FM_local(FM_local~=0))
                intra_cluster_distance(ii) = 0;
            else
                mu       = mean(FM_local(FM_local~=0));
                [r,c]    = find(FM_local~=0);
                disp_vec = zeros(numel(r),1);
                for j=1:numel(r)
                    disp_vec(j) = abs(FM_local(r(j),c(j))-mu);
                end
                intra_cluster_distance(ii) = sum(disp_vec)/numel(r);
            end
        end
    end
    
    % Compute all extra-cluster distances
    for ii=1:N_clust
        for jj=1:N_clust
            if jj==ii
            else
                rows     = T==ii;
                cols     = T==jj;
                FM_local =  FM(rows,cols);
                if isempty(FM_local(FM_local~=0))
                    extra_cluster_distance(ii,jj) = 0;
                else
                    mu       = mean(FM_local(FM_local~=0));
                    [r,c]    = find(FM_local~=0);
                    disp_vec = zeros(numel(r),1);
                    for j=1:numel(r)
                        disp_vec(j) = abs(FM_local(r(j),c(j))-mu);
                    end
                    extra_cluster_distance(ii,jj) = sum(disp_vec)/numel(r);
                end
            end
        end
    end
    
    for ii=1:N_clust
        for jj=1:N_clust
            if ii==jj
            else
                R(ii,jj) = (intra_cluster_distance(ii)+intra_cluster_distance(jj))/...
                    (extra_cluster_distance(ii,jj));
            end
        end
    end
    
    for ii=1:N_clust
        D(ii) = max(R(ii,1:end ~= ii));
    end
    
    DB = sum(D)/N_clust;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CASE 2: we average the distances without computing the mean %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
else
    
    % Compute all intra-cluster distances
    for ii=1:N_clust
        % find rows/columns associated with the current cluster
        r_c  = find(T==ii);
        % if cluster has only one element, set intra-cluster distance to
        % zero
        if numel(r_c) == 1
            intra_cluster_distance(ii) = 0;
        else
            % isolate sub-matrix associated with current cluster
            FM_local = FM(r_c,r_c);
            % compute average intra-cluster distance
            if isempty(FM_local(FM_local~=0))
                intra_cluster_distance(ii) = 0;
            else
                intra_cluster_distance(ii) = mean(FM_local(FM_local~=0));
            end
        end
    end
    
    % Compute all extra-cluster distances
    for ii=1:N_clust
        for jj=1:N_clust
            if jj==ii
            else
                rows     = T==ii;
                cols     = T==jj;
                FM_local =  FM(rows,cols);
                if isempty(FM_local(FM_local~=0))
                    extra_cluster_distance(ii,jj) = 0;
                else
                    extra_cluster_distance(ii,jj) = mean(FM_local(FM_local~=0));
                end
            end
        end
    end
    
    for ii=1:N_clust
        for jj=1:N_clust
            if ii==jj
            else
                R(ii,jj) = (intra_cluster_distance(ii)+intra_cluster_distance(jj))/...
                    (extra_cluster_distance(ii,jj));
            end
        end
    end
    
    for ii=1:N_clust
        D(ii) = max(R(ii,1:end ~= ii));
    end
    
    DB = sum(D)/N_clust;
    
end

return
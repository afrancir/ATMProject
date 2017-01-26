% Written by Alessandro Bombelli, June 07, 2016
% Testing the gap statistics for a set of lat-lon pairs

clc
clear all
close all
tic

xc   = 60;
yc   = 20;

r_lb = 30;
r_ub = 30;

theta1_lb = 30*pi/180;
theta1_ub = 50*pi/180;

theta2_lb = 120*pi/180;
theta2_ub = 140*pi/180;

theta3_lb = 230*pi/180;
theta3_ub = 250*pi/180;

theta4_lb = 290*pi/180;
theta4_ub = 310*pi/180;

N = 100;

C1 = zeros(N,2);
C2 = zeros(N,2);
C3 = zeros(N,2);
C4 = zeros(N,2);

for i=1:N
    r1      = rand(1);
    r2      = rand(1);
    r3      = rand(1);
    r4      = rand(1);
    t1      = rand(1);
    t2      = rand(1);
    t3      = rand(1);
    t4      = rand(1);
    C1(i,1) = xc+((r_ub-r_lb)*r1+r_lb)*cos((theta1_ub-theta1_lb)*t1+theta1_lb);
    C1(i,2) = yc+((r_ub-r_lb)*r1+r_lb)*sin((theta1_ub-theta1_lb)*t1+theta1_lb);
    C2(i,1) = xc+((r_ub-r_lb)*r2+r_lb)*cos((theta2_ub-theta2_lb)*t2+theta2_lb);
    C2(i,2) = yc+((r_ub-r_lb)*r2+r_lb)*sin((theta2_ub-theta2_lb)*t2+theta2_lb);
    C3(i,1) = xc+((r_ub-r_lb)*r3+r_lb)*cos((theta3_ub-theta3_lb)*t3+theta3_lb);
    C3(i,2) = yc+((r_ub-r_lb)*r3+r_lb)*sin((theta3_ub-theta3_lb)*t3+theta3_lb);
    C4(i,1) = xc+((r_ub-r_lb)*r4+r_lb)*cos((theta4_ub-theta4_lb)*t4+theta4_lb);
    C4(i,2) = yc+((r_ub-r_lb)*r4+r_lb)*sin((theta4_ub-theta4_lb)*t4+theta4_lb);
end

figure()
hold on
plot(C1(:,1),C1(:,2),'Color','b','Marker','*','Linestyle','none','Linewidth',2)
plot(C2(:,1),C2(:,2),'Color','b','Marker','*','Linestyle','none','Linewidth',2)
plot(C3(:,1),C3(:,2),'Color','b','Marker','*','Linestyle','none','Linewidth',2)
plot(C4(:,1),C4(:,2),'Color','b','Marker','*','Linestyle','none','Linewidth',2)
grid on


%X       = [C1;C2;C3;C4];
X       = C1;
X       = [X(:,1)-mean(X(:,1)) X(:,2)-mean(X(:,2))];
[U,D,V] = svd(X);
Xp      = X*V;
Zp(:,1) = unifrnd(min(Xp(:,1)),max(Xp(:,1)),numel(Xp(:,1))*100,1);
Zp(:,2) = unifrnd(min(Xp(:,2)),max(Xp(:,2)),numel(Xp(:,2))*100,1);
Z       = Zp*V';

alpha = 30;
figure()
hold on
%plot(Z(:,1),Z(:,2),'Color','r','Marker','o','Linestyle','none')
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

N_clust = 1:6;
W       = zeros(numel(N_clust),1);


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
        N_traj     = numel(FM(:,1));
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
B     = 2000;
ref_W = zeros(numel(N_clust),B);

for i=1:numel(N_clust)
    
    % Number of clusters
    N_c     = N_clust(i);
    
    disp(['Considering the case with ',num2str(N_c),' cluster(s)'])
    
    for j=1:B
        % Generate Monte Carlo sample
        %Zp(:,1) = unifrnd(min(Xp(:,1)),max(Xp(:,1)),numel(Xp(:,1)),1);
        %Zp(:,2) = unifrnd(min(Xp(:,2)),max(Xp(:,2)),numel(Xp(:,2)),1);
        %X_B     = Zp*V';
        idx_perm = randperm(numel(Z(:,1)));
        X_B      = Z(idx_perm(1:numel(X(:,1))),:);

        
        FM_B     = zeros(numel(X_B(:,1)),numel(X_B(:,1)));
        
        for ii=1:numel(X_B(:,1))-1
            for jj=ii+1:numel(X_B(:,1))
                FM_B(ii,jj) = sqrt((X_B(ii,1)-X_B(jj,1))^2+((X_B(ii,2)-X_B(jj,2))^2));
            end
        end
        
        dist_vec_ref   = squareform(FM_B+FM_B');
        if N_c == 1
            n_r        = numel(FM_B(:,1));
            D_r        = sum(FM_B(FM_B~=0));
            ref_W(i,j) = D_r/(2*n_r);
            % More than 1 cluster
        else
            W_temp = zeros(N_c,1);
            N_traj = numel(FM_B(:,1));
            dendrogram_ref = linkage(dist_vec_ref,'average');
            T              = cluster(dendrogram_ref,'maxclust',N_c);
            
            
            for ii=1:N_c
                % find rows/columns associated with the current cluster
                r_c  = find(T==ii);
                n_r  = numel(r_c);
                % if cluster has only one element, set intra-cluster distance to
                % zero
                if numel(r_c) == 1
                    D_r          = 0;
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

%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
hold on
plot(N_clust,W,'Color','b','Marker','o','Linestyle','-','Linewidth',2)
grid on
xlabel('Number of clusters','Fontname','Avantgarde','Fontsize',14)
ylabel('Average intra-cluster distance','Fontname','Avantgarde','Fontsize',14)
%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
hold on
plot(N_clust,Wlog,'Color','b','Marker','o','Linestyle','-','Linewidth',2)
plot(N_clust,Wlog_exp,'Color','r','Marker','o','Linestyle','-','Linewidth',2)
grid on
xlabel('Number of clusters','Fontname','Avantgarde','Fontsize',14)
ylabel('Obs. and exp. log(Wk)','Fontname','Avantgarde','Fontsize',14)
%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
hold on
errorbar(N_clust,GAP,sk,'Color','b','Marker','o','Linestyle','-','Linewidth',2)
grid on
xlabel('Number of clusters','Fontname','Avantgarde','Fontsize',14)
ylabel('Gap','Fontname','Avantgarde','Fontsize',14)

[~,idx]     = max(index);
opt_n_clust = N_clust(idx);
T           = cluster(dendrogram,'maxclust',opt_n_clust);

figure()
hold on
for i=1:opt_n_clust
    idx_points = find(T==i);
    rgb        = [rand(1) rand(1) rand(1)];
    plot(X(idx_points,1),X(idx_points,2),'Color',rgb,'Marker','x','Markersize',8,'Linestyle','none')
end




time = toc;
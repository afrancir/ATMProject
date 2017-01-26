function plot_unitary_square(Bmin,Bmax,Lmin,Lmax,P,Q,P_length,Q_length,P_us,Q_us,alpha,FD,tail_or_head)

figure()
hold on

p = length(Lmin(1,:))-1;
q = length(Bmin(:,1))-1;

for i=1:q+1
    plot([0 p],[i-1 i-1],'Color','k','Linewidth',1)
end
for i=1:p+1
    plot([i-1 i-1],[0 q],'Color','k','Linewidth',1)
end


for i=1:numel(Bmin(:,1))
    for j=1:numel(Bmin(1,:))
        if Bmin(i,j)==-1
        else
            plot([Bmin(i,j),Bmax(i,j)]+j-1,[numel(Bmin(:,1))-i,numel(Bmin(:,1))-i],'Color','g','Linewidth',5)
        end
    end
end

for i=1:numel(Lmin(:,1))
    for j=1:numel(Lmin(1,:))
        if Lmin(i,j)==-1
        else
            plot([j-1,j-1],[numel(Lmin(:,1))-i+(1-Lmin(i,j)),numel(Lmin(:,1))-i+(1-Lmax(i,j))],'Color','g','Linewidth',5)
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q_tail           = Q(Q_us,:);
Q_head           = Q(Q_us+1,:);
P_tail           = P(P_us,:);
P_head           = P(P_us+1,:);
Q_tail_cart      = latlon2cart(Q_tail,1);
Q_head_cart      = latlon2cart(Q_head,1);
P_tail_cart      = latlon2cart(P_tail,1);
P_head_cart      = latlon2cart(P_head,1);
DQ               = Q_length(Q_us);
DP               = P_length(P_us);
N_f              = 40;
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

for ii=1:numel(f_vec_Q)
    for jj=1:numel(f_vec_P)
        if dist_mat(ii,jj) <= alpha*FD
            plot(P_us+f_vec_P(jj)-1,numel(Q_length)-Q_us-f_vec_Q(ii)+1,'Color','g','Marker','*','Markersize',2,'Linewidth',1.3)
        else
            plot(P_us+f_vec_P(jj)-1,numel(Q_length)-Q_us-f_vec_Q(ii)+1,'Color','r','Marker','*','Markersize',2,'Linewidth',1.3)
        end
    end
end

[r,c]    = find(dist_mat<=alpha*FD);

if ~isempty(r)
    dist_tail_vec = zeros(numel(r),1);
    dist_head_vec = zeros(numel(r),1);
    for ii=1:numel(r)
        dist_tail_vec(ii) = dist_tail(r(ii),c(ii));
        dist_head_vec(ii) = dist_head(r(ii),c(ii));
    end
    
    % we are interested in the max distance from the head point (thus the
    % point we are considering is the initial point of a sequence)
    if tail_or_head
        [~,idx]       = max(dist_head_vec);
    % we are interested in the max distance from the tail point (thus the
    % point we are considering is the endpoint of a sequence)    
    else
        [~,idx]       = max(dist_tail_vec);
    end
    best_f_Q      = f_vec_Q(r(idx));
    best_f_P      = f_vec_P(c(idx));
    [~,idx_f_Q]   = find(f_vec_Q==best_f_Q);
    [~,idx_f_P]   = find(f_vec_P==best_f_P);
    plot(P_us+f_vec_P(idx_f_P)-1,numel(Q_length)-Q_us-f_vec_Q(idx_f_Q)+1,'Color','b','Marker','s','Markersize',4,'Linewidth',2)
else
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



set(gca,'XTick',0.5:p-0.5)
set(gca,'YTick',0.5:q-0.5)
XTickVec = cell(q,1);
for i=1:p
    XTickVec{i} = num2str(i);
end
YTickVec = cell(q+1,1);
for i=1:q
    YTickVec{i} = num2str(q+1-i);
end
set(gca,'XTickLabel',XTickVec)
set(gca,'YTickLabel',YTickVec)
axis([0 p 0 q])
set(gca,'TickLength',[0 0])
set(gca,'xaxisLocation','top')
xlabel('Flight track P','Fontname','Avantgarde','Fontsize',14)
ylabel('Flight track Q','Fontname','Avantgarde','Fontsize',14)

return
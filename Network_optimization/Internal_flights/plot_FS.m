function plot_FS(Bmin,Bmax,Lmin,Lmax)

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

% text(3,3,num2str(99),'Fontname','Avantgarde','Fontsize',14)
% text(3,3,num2str(100),'Fontname','Avantgarde','Fontsize',14)
% text(3,3,num2str(128),'Fontname','Avantgarde','Fontsize',14)
% text(3,3,num2str(129),'Fontname','Avantgarde','Fontsize',14)

%X = [0:13 13.5 14 15];
%Y = [10 9 8 7.5 7 6.5 6 5.5 5 4 3 2.5 2 1.5 1 0.5 0];
%plot(X,Y,'Color','r','Linewidth',2.5,'Marker','x','Markersize',10)

%X = [0 0.5 0.6 0.7 1 2 3 4 5 6 7 7.1 7.6 7.62 7.63 7.7 7.72 7.75 7.8 7.82 7.88 7.9 7.95 8 8.5 9 10 11 11.7 12 13];
%Y = [29 28 27 26 25 24 23 22 21 20 19.5 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1 0];
%plot(X,Y,'Color','r','Linewidth',2.5)

%X = [0 1 2 3 4 5 6 7 7.1 7.15 7.2 7.3 7.35 7.4 7.5 7.55 7.75 7.88 7.9 7.95 8 9 10 11 12 13];
%Y = [22 21.5 21 20 19 18 17.5 17.4 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1 0];
%plot(X,Y,'Color','r','Linewidth',2.5)

% set(gca,'XTick',0:p)
% set(gca,'YTick',0:q)
% YTickVec = cell(q+1,1);
% for i=1:q+1
%     YTickVec{i} = num2str(q+1-i);
% end
% set(gca,'YTickLabel',YTickVec)
% axis([0 p 0 q])

set(gca,'XTick',0.5:p-0.5)
set(gca,'YTick',0.5:q-0.5)
XTickVec = cell(q,1);
% for i=1:q
%     XTickVec{i} = num2str(i);
% end
% YTickVec = cell(q+1,1);
% for i=1:q+1
%     YTickVec{i} = num2str(q+1-i);
% end
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
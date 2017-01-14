%%%
XMIN = min(min(weather_longitude));
YMIN = min(min(weather_latitude));
XMAX = max(max(weather_longitude));
YMAX = max(max(weather_latitude));

Polygons  = allLatlonCH{end-5};
nPolygons = length(Polygons);

figure(1899)
clf
axis([XMIN XMAX YMIN YMAX])
grid on 
hold on

for pp = 1:nPolygons
    Polygons{pp} = [Polygons{pp};Polygons{pp}(1,:)]; %for a closed polygon plot
%     plot(Polygons{pp}(:,2),Polygons{pp}(:,1),'LineStyle','-.')
    if size(Polygons{pp},1) < 4 %for 1 cell avoidance. remember head and tail are the same so 2 are only 1 point.
        % 3 are 2 points but they display poorly if plotted
        scatter(Polygons{pp}(:,2),Polygons{pp}(:,1),40,[0 0 0],'filled')
        scatter(Polygons{pp}(:,2),Polygons{pp}(:,1),10,[1 1 0],'filled')
        hold on
        plot(Polygons{pp}(:,2),Polygons{pp}(:,1),'Color',[0 0 0],'linewidth',3)
        plot(Polygons{pp}(:,2),Polygons{pp}(:,1),'k--', 'Color',[1 1 0],'linewidth',2)
        hold on
    else 
        plot(Polygons{pp}(:,2),Polygons{pp}(:,1),'Color',[0 0 0],'linewidth',3)
        plot(Polygons{pp}(:,2),Polygons{pp}(:,1),'k--', 'Color',[1 1 0],'linewidth',2)
        hold on
    end 
end 

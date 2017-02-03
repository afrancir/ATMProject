XMIN = LON_lim(1);
YMIN = LAT_lim(1);
XMAX = LON_lim(2);
YMAX = LAT_lim(2);

Polygons  = allLatlonCH{1};
nPolygons = length(Polygons);

figure(1899)
clf
axis([XMIN XMAX YMIN YMAX])
grid on 
hold on

for pp = 1:nPolygons
    %%%plot(lon,lat)
    plot(Polygons{pp}(:,2),Polygons{pp}(:,1),'Color',[0 0 0],'linewidth',3)
    plot(Polygons{pp}(:,2),Polygons{pp}(:,1),'k--', 'Color',[1 1 0],'linewidth',2)
    hold on
end 

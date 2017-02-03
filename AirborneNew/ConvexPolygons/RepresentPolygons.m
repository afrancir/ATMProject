function pngname = RepresentPolygons(weather_longitude,weather_latitude,jj,Polygons)
  figure(jj) 
    clf
    contour(weather_longitude,weather_latitude,Polygons,'Linewidth',2)
    hold on
    colormap jet
    c=colorbar;
    ylabel(c,'Normalized Weather Product','Fontsize',14,'Fontname','Avantgarde')
    axis equal
    grid on        
    title(['T =  ',num2str(jj)])
%     title('Initial Stage')
    xlabel('Longitude [deg]','Fontsize',14,'Fontname','Avantgarde')
    ylabel('Latitude [deg]','Fontsize',14,'Fontname','Avantgarde')
    pngname = ['C:\Users\Arnau\Desktop\Matlab figures\Timestep',num2str(jj),'.png'];
%     pngname = ['C:\Users\Arnau\Desktop\Matlab figures\111Initial_Stage.png'];
%    saveas(gcf,pngname);
end 

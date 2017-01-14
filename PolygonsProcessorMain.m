%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%$$               Arnau Francí-Rodon  11-8-2016                        $$%                                                                    $$%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% description: runs PolygonClassifier and PolygonsConvexHull until the 
% polygons are stabilized and the new iterations do not alter their shape. 
%
% INPUT: weatherPolygons is an nxm binary matrix with values of weather
% product forecast > threshold.
% OUTPUT: definitiveHullPolygons is the definitive nxm binary matrix with
% convex weather polygons that we will use for the planning algorithms


%%% p = genpath(pwd) at this to the beginning of the MAIn file so every
%%% subfolder is accessible and ConvexPolygons can be run from the parent
%%% folder (and keep the files organized that way)

weatherPolygons = weather_polygons(:,:,10); % or weather_polygons(:,:,x)
polygonsMatrix     = weatherPolygons;
pastPolygonsMatrix = weatherPolygons.*0; %initialize it 
nIterations        = 0;

while ~ isequal(polygonsMatrix,pastPolygonsMatrix)
   PolygonsClassified = PolygonClassifierEvolution(polygonsMatrix);
   PolygonsConvex     = PolygonsConvexHull(PolygonsClassified);
   
   pastPolygonsMatrix = polygonsMatrix;
   polygonsMatrix     = PolygonsConvex;
   
   nIterations        = nIterations + 1; %check how many does it take
end 
%%% output the version with the Polygons separated from the definitive
%%% iteration.
[~,polygonsMatrixSeparated] = PolygonsConvexHull(PolygonsClassified);

%% Returning Polygons to lat/lon and take Convex Hull
latlonCH = cell(size(polygonsMatrixSeparated,3),1);

for i = 1: size(polygonsMatrixSeparated,3)
    weatherLatitudeConvex  = weather_latitude(polygonsMatrixSeparated(:,:,i)==1); %x
    weatherLongitudeConvex = weather_longitude(polygonsMatrixSeparated(:,:,i)==1);%y
    
    % avoid error of calculating CH of a single point or a pair
    if numel(weatherLatitudeConvex) < 3
        latlonCH{i}        = [weatherLatitudeConvex,weatherLongitudeConvex];  
    else 
        K = convhull(weatherLatitudeConvex,weatherLongitudeConvex);
        weatherLatitudeCH  = weatherLatitudeConvex(K);
        weatherLongitudeCH = weatherLongitudeConvex(K);    
        latlonCH{i}        = [weatherLatitudeCH,weatherLongitudeCH];
    end 
end 
    
%% GRAPHICS Plot and displaying results 
display(nIterations)

figure(177)
clf
subplot(1,2,1) 
contour(weather_longitude,weather_latitude,weatherPolygons,'Linewidth',2)
hold on
colormap jet
c=colorbar;
ylabel(c,'Normalized Weather Product','Fontsize',14,'Fontname','Avantgarde')
axis equal
grid on        
title('Initial Polygons')
xlabel('Longitude [deg]','Fontsize',14,'Fontname','Avantgarde')
ylabel('Latitude [deg]','Fontsize',14,'Fontname','Avantgarde')
subplot(1,2,2) 
contour(weather_longitude,weather_latitude,polygonsMatrix,'Linewidth',2)
hold on
colormap jet
c=colorbar;
ylabel(c,'Normalized Weather Product','Fontsize',14,'Fontname','Avantgarde')
axis equal
grid on      
title('Final Polygons')
xlabel('Longitude [deg]','Fontsize',14,'Fontname','Avantgarde')
ylabel('Latitude [deg]','Fontsize',14,'Fontname','Avantgarde')

%%% In case we want to overlap the convex hull to check it is coincident
%%% with the contour of the original figure: 
% % % figure(gcf)
% % % hold on
% % % for qq = 1:numel(latlonCH)
% % %     plot(latlonCH{qq}(:,2),latlonCH{qq}(:,1),'o','Color',[1 0 1/2],'MarkerSize',10)
% % %     hold on 
% % % end 


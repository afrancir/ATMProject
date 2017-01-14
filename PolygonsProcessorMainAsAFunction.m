%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%$$               Arnau Francí-Rodon  8-29-2016                        $$%                                                                    $$%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function allLatlonCH = PolygonsProcessorMainAsAFunction(weatherPolygons, ...
%     weather_latitude, weather_longitude, graphicFlag) 
%
% Description: runs PolygonClassifier and PolygonsConvexHull until the 
% polygons are stabilized and the new iterations do not alter their shape. 
%
% INPUTS: 
%         - weatherPolygons is an nxm binary matrix with values of weather
%           product forecast > threshold.
%         - graphicFlag TRUE or FALSE for graphic displaying
%
% OUTPUT: 
%         - allLatlonCH. All M convex polygons defined by lat/lon pairs sets 
%           for each of the N timesteps the horizon has.


%%% p = genpath(pwd) at this to the beginning of the Main file so every
%%% subfolder is accessible and ConvexPolygons can be run from the parent
%%% folder (and keep the files organized that way)

function allLatlonCH = PolygonsProcessorMainAsAFunction(weatherPolygons, ...
    weather_latitude, weather_longitude, graphicFlag) 

nTimeSteps = size(weatherPolygons,3);
allLatlonCH  = cell(1,nTimeSteps);

for jj = 1:nTimeSteps
    
    polygonsMatrix     = weatherPolygons(:,:,jj);
    pastPolygonsMatrix = weatherPolygons.*0; %initialize it 
    nIterations        = 0;

    while ~ isequal(polygonsMatrix,pastPolygonsMatrix)
       PolygonsClassified = PolygonClassifierEvolution2(polygonsMatrix);
       PolygonsConvex     = PolygonsConvexHull2(PolygonsClassified);

       pastPolygonsMatrix = polygonsMatrix;
       polygonsMatrix     = PolygonsConvex;
       
       nIterations        = nIterations + 1; %check how many does it take
       RepresentPolygons(weather_longitude,weather_latitude,nIterations,PolygonsConvex);
    end 
    %%% output the version with the Polygons separated from the definitive
    %%% iteration.
    [~,polygonsMatrixSeparated] = PolygonsConvexHull2(PolygonsClassified);

    %% Returning Polygons to lat/lon and take only the needed points for representing Convex Hull
    latlonCH = cell(size(polygonsMatrixSeparated,3),1);

    for i = 1: size(polygonsMatrixSeparated,3)
        weatherLatitudeConvex  = weather_latitude(polygonsMatrixSeparated(:,:,i)==1); %x
        weatherLongitudeConvex = weather_longitude(polygonsMatrixSeparated(:,:,i)==1);%y

        % avoid error of calculating CH of a single point or a pair. Notice
        % that to calculate collinearity weatherLongitudeConvex could equally be used 
        if (numel(weatherLatitudeConvex) <= 3 || areCollinear(polygonsMatrixSeparated(:,:,i)))
            latlonCH{i}        = [weatherLatitudeConvex,weatherLongitudeConvex];  
        else 
            K = convhull(weatherLatitudeConvex,weatherLongitudeConvex);
            weatherLatitudeCH  = weatherLatitudeConvex(K);
            weatherLongitudeCH = weatherLongitudeConvex(K);    
            latlonCH{i}        = [weatherLatitudeCH,weatherLongitudeCH];
        end 
    end 
    
    allLatlonCH{1,jj} = latlonCH;
    %% GRAPHICS Plot and displaying results 
    display(nIterations)
%     graphicFlag = false;
    
    if graphicFlag
        figure(jj) %#ok<UNRCH>
        clf
        contour(weather_longitude,weather_latitude,weatherPolygons(:,:,jj),'Linewidth',2)
        hold on
        colormap jet
        c=colorbar;
        ylabel(c,'Normalized Weather Product','Fontsize',14,'Fontname','Avantgarde')
        axis equal
        grid on        
        title(['Iteration  ',num2str(jj)])
        xlabel('Longitude [deg]','Fontsize',14,'Fontname','Avantgarde')
        ylabel('Latitude [deg]','Fontsize',14,'Fontname','Avantgarde')
        pngname = ['C:\Users\Arnau\Desktop\Matlab figures\',num2str(jj),'.png'];
        saveas(gcf,pngname);
% %         figure(jj+100) %#ok<UNRCH>
% %         clf
% %         subplot(1,2,1) 
% %         contour(weather_longitude,weather_latitude,weatherPolygons(:,:,jj),'Linewidth',2)
% %         hold on
% %         colormap jet
% %         c=colorbar;
% %         ylabel(c,'Normalized Weather Product','Fontsize',14,'Fontname','Avantgarde')
% %         axis equal
% %         grid on        
% %         title('Initial Polygons')
% %         xlabel('Longitude [deg]','Fontsize',14,'Fontname','Avantgarde')
% %         ylabel('Latitude [deg]','Fontsize',14,'Fontname','Avantgarde')
% %         subplot(1,2,2) 
% %         contour(weather_longitude,weather_latitude,polygonsMatrix,'Linewidth',2)
% %         hold on
% %         colormap jet
% %         c=colorbar;
% %         ylabel(c,'Normalized Weather Product','Fontsize',14,'Fontname','Avantgarde')
% %         axis equal
% %         grid on      
% %         title('Final Polygons')
% %         xlabel('Longitude [deg]','Fontsize',14,'Fontname','Avantgarde')
% %         ylabel('Latitude [deg]','Fontsize',14,'Fontname','Avantgarde')
% %         figname = ['C:\Users\Arnau\Desktop\Matlab figures\',num2str(jj+100),'.fig'];
% %         savefig(figname);
% %         pngname = ['C:\Users\Arnau\Desktop\Matlab figures\',num2str(jj+100),'.png'];
% %         saveas(gcf,pngname);
        %%% In case we want to overlap the convex hull to check it is coincident
        %%% with the contour of the original figure: 
        % % % figure(gcf)
        % % % hold on
        % % % for qq = 1:numel(latlonCH)
        % % %     plot(latlonCH{qq}(:,2),latlonCH{qq}(:,1),'o','Color',[1 0 1/2],'MarkerSize',10)
        % % %     hold on 
        % % % end 
    end %endIf    
end %endFor
end %function
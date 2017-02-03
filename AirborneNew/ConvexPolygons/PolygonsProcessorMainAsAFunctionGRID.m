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

function allLatlonCH = PolygonsProcessorMainAsAFunctionGRID(weatherPolygons, ...
    GRID) 

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
    end 
    %%% output the version with the Polygons separated from the definitive
    %%% iteration.
    [~,polygonsMatrixSeparated] = PolygonsConvexHull2(PolygonsClassified);

    %% Returning Polygons to lat/lon and take Convex Hull
    latlonCH = cell(size(polygonsMatrixSeparated,3),1);

    for i = 1: size(polygonsMatrixSeparated,3)
        convex1Indices          =  find(polygonsMatrixSeparated(:,:,i)==1);
        convex1IndicesLength    =  numel(convex1Indices);
        weatherLatitudeConvex  =  zeros(4*convex1IndicesLength,1); %4 different lat/lon for every cell element
        weatherLongitudeConvex =  zeros(4*convex1IndicesLength,1);
        
        for w = 4:4:4*convex1IndicesLength
            weatherLatitudeConvex(w-3)    = GRID{convex1Indices(w/4)}(1,2);
            weatherLongitudeConvex(w-3)   = GRID{convex1Indices(w/4)}(1,1);
            weatherLatitudeConvex(w-2)  = GRID{convex1Indices(w/4)}(2,2);
            weatherLongitudeConvex(w-2) = GRID{convex1Indices(w/4)}(2,1);
            weatherLatitudeConvex(w-1)  = GRID{convex1Indices(w/4)}(3,2);
            weatherLongitudeConvex(w-1) = GRID{convex1Indices(w/4)}(3,1);
            weatherLatitudeConvex(w)  = GRID{convex1Indices(w/4)}(4,2);
            weatherLongitudeConvex(w) = GRID{convex1Indices(w/4)}(4,1);          
        end 
        K = convhull(weatherLatitudeConvex,weatherLongitudeConvex);
        weatherLatitudeCH  = weatherLatitudeConvex(K);
        weatherLongitudeCH = weatherLongitudeConvex(K);    
        latlonCH{i}        = [weatherLatitudeCH,weatherLongitudeCH];
    end 
    
    allLatlonCH{1,jj} = latlonCH;
end %endFor
end %function
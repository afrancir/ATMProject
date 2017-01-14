function [PolygonsConvex,PolygonsConvexApart] = PolygonsConvexHull2(polygonsClassified)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%$$               Arnau Francí-Rodon  8-9-2016 (m. 8-23)               $$%                                                                    $$%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function PolygonsConvex = PolygonsConvexHull2(polygonsClassified)
% description: It calculates the convex hull of every polygon labeled as
% independent separately, that is each polygon is isolated and then all
% them are overlapped back in their original position
%
% INPUT:  PolygonsClassified is an nxm integer matrix with numbers
% as polygons'Poly ID and 0 for non-polygonal pixels.
% OUTPUT: PolygonsConvex is an nxm binary matrix with the convex hull of
% all weather avoidance zones.
%
%
% Example:
%   PolygonsClassified = [1 0 2; 0 0 2; 2 2 2] 
%   --->  PolygonsConvex = [1 0 0; 0 0 0; 0 0 0] + [0 0 1; 0 1 1; 1 1 1] 
%                        = [1 0 1; 0 1 1; 1 1 1]

%%%This first part is considerably inefficient. In a future it should be
%%%improved changing the part of the classCounter in the
%%%PolygonClassifierEvolution function .

differentPolygonsClass = unique(polygonsClassified) ; 
differentPolygonsClass = differentPolygonsClass(2:end); %0 doesn't count
nPolygons = numel(differentPolygonsClass);
%%% Then we will take the convexhull of every one of the different polygons

PolygonsApart    = zeros([size(polygonsClassified),nPolygons]);
ConvexHullsTogether = zeros([size(polygonsClassified),1]); %for graphics purposes
ConvexHullsApart = zeros([size(polygonsClassified),nPolygons]); 

for k = 1:nPolygons
    PolygonsApart(:,:,k) = (polygonsClassified == differentPolygonsClass(k));
    ConvexHullsApart(:,:,k) = bwconvhull(PolygonsApart(:,:,k));
    ConvexHullsTogether = or(ConvexHullsTogether,ConvexHullsApart(:,:,k));
end 

PolygonsConvex = ConvexHullsTogether; %output
PolygonsConvexApart = ConvexHullsApart; %output2

end 
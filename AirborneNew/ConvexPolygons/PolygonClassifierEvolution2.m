function PolygonsClassified = PolygonClassifierEvolution2(polygonsMatrix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%$$               Arnau Francí-Rodon  8-25-2016                        $$                                                                    $$%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function PolygonsClassified = PolygonClassifier(polygonsMatrix)
% description: Closed Polygon detector, similar to flood fill algorithms.
%
% INPUT: polygonsMatrix is an nxm binary matrix 
% OUTPUT: PolygonsClassified is an nxm integer matrix with numbers
% as polygons'Poly ID an NaN for non-polygonal pixels.
%
%
% Example:
%   polygonsMatrix = [1 0 1; 0 0 1; 1 1 1] 
%   --->  PolygonsClassified = [1 0 2; 0 0 2; 2 2 2] 

[row,col]          = size(polygonsMatrix);
PolygonsClassified = zeros(row,col);
classCounter       = 0; 

for i = 1:row
    for j= 1:col
        if(polygonsMatrix(i,j)) %that is, the element is 1 and not 0         
            top = 0; left = 0;  %top, left initialization         
            % checks to avoid exceeding matrix boundaries
            topCorner      = (i==1);
            firstColumn    = (j==1);            
            if ~(topCorner) 
                top  = polygonsMatrix(i-1,j);
            end             
            if ~(firstColumn) 
                left = polygonsMatrix(i,j-1);
            end 
                                    
            % check if any neighbour is not null 
                % make a construct for the switch case. Notice that 0 is
                % not possible since it would not enter the condition of if
                %(anyCorner)
                switcher = top + 2*left;
                
                switch switcher    
                    case 0 %corners were not visited yet. New polygon
                         classCounter = classCounter + 1;
                         PolygonsClassified(i,j) = classCounter;
                    case 1 %only the top cell is 1, left is 0
                        topClass  = PolygonsClassified(i-1,j);
                        PolygonsClassified(i,j) = topClass;
                    case 2 %left cell is 1, top is 0
                        leftClass = PolygonsClassified(i,j-1);
                        PolygonsClassified(i,j) = leftClass;                        
                    case 3 % both are 1
                        leftClass = PolygonsClassified(i,j-1);
                        topClass  = PolygonsClassified(i-1,j);
                        
                        if(leftClass==topClass)
                            PolygonsClassified(i,j) = leftClass; %or topClass
                        else %two polygons are in fact connected so same
                            minClass = min([topClass,leftClass]);   
                            PolygonsClassified(i,j) = minClass;

                            lS    = 0; %leftShift
                            upS   = 0; %upShift
                            % correct the polygon Class mismatchment
                            while ~isequal(topClass,leftClass,minClass)
                                [~,index] = max([topClass,leftClass,minClass]);
                                switch index
                                    case 1 %move up
                                        upS = upS + 1;
                                    case 2 %move leftwards
                                        lS  = lS + 1;
                                    case 3 %we reached a border
                                        %little check. Both MUST be 0. 
                                        if any([topClass,leftClass])
                                            display('error in the code 2')
                                        end%if
                                        break %break the while loop  
                                end                                 
                                PolygonsClassified(i-upS,j-lS) = minClass;
                                topClass  = PolygonsClassified(i-(1+upS),j-lS); 
                                leftClass = PolygonsClassified(i-upS,j-(1+lS));
                            end 
                        end 
                    otherwise
                        display('error in the code')                                        
                end                                      
        end     	
    end 
end 
end % end function



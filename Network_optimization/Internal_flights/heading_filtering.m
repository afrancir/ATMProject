function good_idx = heading_filtering(P,angle)

% Input:  - P     : (Nx2) matrix with lat/lon pairs (deg) that define a
%                   flight track
%         - angle : threshold (deg) below which three consecutive lat/lon
%                   pairs are assumed to form a single great circle segments
% Output: - new_P : (nx2) with n<=N matrix with lat/lon pairs (deg) that
%                   define the reduced trajectory

% Switch from deg to rad
P = P*pi/180;

% Initialize the counter and the vector that stores the indices we want to
% keep
cont = 1;
good_idx = zeros(length(P(:,1)),1);

track_is_complete = 1;

while track_is_complete~=0
    
    good_idx(cont) = 1;
    
    % Current lat/lon pair and next lat/lon pair. We compute the heading
    % that defines this great cirle, and that will be used as a reference
    lat_a = P(cont,1);
    lon_a = P(cont,2);
    lat_b = P(cont+1,1);
    lon_b = P(cont+1,2);
    
    X           = cos(lat_b)*sin(lon_b-lon_a);
    Y           = cos(lat_a)*sin(lat_b)-sin(lat_a)*cos(lat_b)*cos(lon_b-lon_a);
    ref_heading = atan2(X,Y)*180/pi;
    
    
    
    % Normalizing heading to be in the [0,360 deg] interval
    if ref_heading >= 360
        ref_heading = ref_heading-360*floor(ref_heading/360);
    else
    end
    % Normalizing heading to be in the [0,360 deg] interval
    if ref_heading < 0
        ref_heading = ref_heading+360;
    else
    end
    
    is_below_threshold = 1;
    cont2              = 1;
    
    while is_below_threshold~=0
        lat_b       = P(cont+1+cont2,1);
        lon_b       = P(cont+1+cont2,2);
        X           = cos(lat_b)*sin(lon_b-lon_a);
        Y           = cos(lat_a)*sin(lat_b)-sin(lat_a)*cos(lat_b)*cos(lon_b-lon_a);
        new_heading = atan2(X,Y)*180/pi;
        
        % Normalizing heading to be in the [0,360 deg] interval
        if new_heading >= 360
            new_heading = new_heading-360*floor(new_heading/360);
        else
        end
        % Normalizing heading to be in the [0,360 deg] interval
        if new_heading < 0
            new_heading = new_heading+360;
        else
        end
        
        % The heading between these entries exceeds the threshold. The last
        % point checked is the new starting point
        if abs(new_heading-ref_heading) > angle
            is_below_threshold = 0;
            cont               = cont+cont2;
            cont2              = 1;
            if cont == length(P(:,1))-1
                good_idx(cont)    = 1;
                good_idx(cont+1)  = 1;
                track_is_complete = 0;
            else
            end
            % Otherwise, we search if the next point still satisfies the
            % requirement or not
        else
            if cont+1+cont2 == length(P(:,1))
                good_idx(cont)    = 1;
                good_idx(end)     = 1;
                track_is_complete = 0;
                break;
            else
                cont2 = cont2+1;
            end
        end
        
    end
    
end

return
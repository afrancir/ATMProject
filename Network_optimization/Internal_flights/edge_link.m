function [ROW]=edge_link(ID,p,q,Lmin,Lmax,Bmin,Bmax)

N   = (2*p+1)*(q+1)+(p+1)*q;
ROW = zeros(1,N);

BLOCK = ceil(ID/(3*p+2));
POS   = ID-(BLOCK-1)*(3*p+2);
R     = ceil(ID/(3*p+2));
C     = ceil(POS/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Corner nodes of last row (apart from destination node) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ID>q*(3*p+2) && mod(POS,2)~=0 && ID~=(2*p+1)*(q+1)+(p+1)*q
    % Next horizontal node to the right
    if Bmin(R,C)==0 && Bmax(R,C)>0
        ROW(1,ID+1) = 1;
    else
    end
    % Next corner node to the right
    if Bmin(R,C)==0 && Bmax(R,C)==1
        ROW(1,ID+2) = 1;
    else
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Corner nodes of last column (apart from destination node) %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif POS==2*p+1 && BLOCK~=q+1
    if Lmin(R,C)==0 && Lmax(R,C)>0
        ROW(1,ID+p+1) = 1;
    else
    end
    if Lmin(R,C)==0 && Lmax(R,C)==1
        ROW(1,ID+3*p+2) = 1;
    else
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    %%% Destination node %%%
    %%%%%%%%%%%%%%%%%%%%%%%%
elseif ID==(2*p+1)*(q+1)+(p+1)*q
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Any internal corner node %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    
    % Next horizontal node to the right
    if Bmin(R,C)==0 && Bmax(R,C)>0
        ROW(1,ID+1) = 1;
    else
    end
    % Next vertical node below
    if Lmin(R,C)==0 && Lmax(R,C)>0
        ROW(1,ID+2*p+1-(C-1)) = 1;
    else
    end
    % Next vertical node to the right
    if (Lmax(R,C+1)>0) && (Lmin(R,C)==0 || Bmin(R,C)==0)
        ROW(1,ID+2*p+1-(C-1)+1) = 1;
    else
    end
    % Next horizontal node below
    if (Bmin(R,C)==0 || Lmin(R,C)==0) && (Bmax(R+1,C)>0)
        ROW(1,ID+3*p+3) = 1;
    else
    end
    % Corner node diagonally opposite
    if (Bmax(R+1,C)==1 || Lmax(R,C+1)==1) && (Bmin(R,C)==0 || Lmin(R,C)==0)
        ROW(1,ID+3*p+4) = 1;
    else
    end
    % Next horizontal corner node
    if (Bmin(R,C)==0 && Bmax(R,C)==1)
        ROW(1,ID+2) = 1;
    else
    end
    % Next vertical corner node
    if (Lmin(R,C)==0 && Lmax(R,C)==1)
        ROW(1,ID+3*p+2) = 1;
    else
    end
end

return
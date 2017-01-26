function [ROW]=hor_node_link(ID,p,q,Lmin,Lmax,Bmin,Bmax)

N   = (2*p+1)*(q+1)+(p+1)*q;
ROW = zeros(1,N);

BLOCK = ceil(ID/(3*p+2));
POS   = ID-(BLOCK-1)*(3*p+2);
R     = ceil(ID/(3*p+2));
C     = ceil(POS/2);

if ID>q*(3*p+2) && mod(POS,2)==0
    %disp(['Node ',num2str(ID), ' is an internal node along a horizontal egde on the lower boundary of the grid'])
    if Bmax(R,C)==1
        %disp(['Connection between node ',num2str(ID),' and node ',num2str(ID+1)])
        ROW(1,ID+1) = 1;
    else
    end
else
    %disp(['Node ',num2str(ID), ' is an internal node along a horizontal egde'])
    %% Connection 1
    
    %if Bmax(R,C)==1
    if Bmin(R,C)>=0 && Bmax(R,C)==1    
        %disp(['Connection between node ',num2str(ID),' and node ',num2str(ID+1)])
        ROW(1,ID+1) = 1;
    else
    end
    
    %% Connection 2
    
    %if (Lmax(R,C+1)~=0 || Lmin(R,C+1)~=1) && (Lmin(R,C+1)>=0) && Bmin(R,C)>=0
    if Bmin(R,C)>=0 && Lmin(R,C+1)>=0 && Lmax(R,C+1)>=0 
        %disp(['Connection between node ',num2str(ID),' and node ',num2str(ID+2*p+2-POS/2)])
        ROW(1,ID+2*p+2-POS/2) = 1;
    else
    end
    
    %% Connection 3
    
    %if (Bmax(R+1,C)>=Bmax(R,C) && Bmax(R+1,C)>0 && Bmin(R+1,C)<1 && Bmin(R,C)>=0)
    if Bmin(R,C)>=0 && Bmax(R+1,C)>=Bmin(R,C) && Bmin(R+1,C)>=0
        %disp(['Connection between node ',num2str(ID),' and node ',num2str(ID+3*p+2)])
        ROW(1,ID+3*p+2) = 1;
    else
    end
    
    %% Connection 4
    
    %if (Bmax(R+1,C)==1 || Lmax(R,C+1)==1) && ((Bmin(R,C)>=0 && Bmax(R,C)>Bmin(R,C)) || (Bmin(R,C)==Bmax(R,C) && (Bmin(R,C)~=0 && Bmin(R,C)~=1) && Bmin(R,C)>=0))
    if Bmin(R,C)>=0 && (Bmax(R+1,C)==1 || Lmax(R,C+1)==1) 
        %disp(['Connection between node ',num2str(ID),' and node ',num2str(ID+3*p+3)])
        ROW(1,ID+3*p+3) = 1;
    else
    end
    
end

return
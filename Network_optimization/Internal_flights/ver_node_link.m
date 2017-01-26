function [ROW]=ver_node_link(ID,p,q,Lmin,Lmax,Bmin,Bmax)

N   = (2*p+1)*(q+1)+(p+1)*q;
ROW = zeros(1,N);

BLOCK = ceil(ID/(3*p+2));
POS   = ID-(BLOCK-1)*(3*p+2);
R     = ceil(ID/(3*p+2));
C     = ID-(BLOCK-1)*(3*p+2)-(2*p+1);

if POS==3*p+2
    %disp(['Node ',num2str(ID), ' is an internal node along a vertical edge on the right boundary of the grid'])
    if Lmax(R,C)==1 && Lmin(R,C)<1
        %disp(['Connection between node ',num2str(ID),' and node ',num2str(ID+2*p+1)])
        ROW(1,ID+2*p+1) = 1;
    else
    end
else
    %disp(['Node ',num2str(ID), ' is an internal node along a vertical edge'])
    %% Connection 1
    
    %if (Lmax(R,C+1)>=Lmax(R,C) && Lmax(R,C+1)>0 && Lmin(R,C+1)<1 && Lmin(R,C)>=0 && Lmin(R,C+1)~=-1)
    if  Lmin(R,C)>=0 && Lmax(R,C+1)>=Lmin(R,C) && Lmin(R,C+1)>=0
        %disp(['Connection between node ',num2str(ID),' and node ',num2str(ID+1)])
        ROW(1,ID+1) = 1;
    else
    end
    
    %% Connection 2
    
    if Lmin(R,C)>=0 && Lmax(R,C)==1 
        %disp(['Connection between node ',num2str(ID),' and node ',num2str(ID+POS-p-1)])
        ROW(1,ID+POS-p-1) = 1;
    else
    end
    
    %% Connection 3
    
    %if (Bmax(R+1,C)~=0 || Bmin(R+1,C)~=1) && (Bmin(R+1,C)>=0) && Lmin(R,C)>=0
    if  Lmin(R,C)>=0  && Bmin(R+1,C)>=0
        %disp(['Connection between node ',num2str(ID),' and node ',num2str(ID+POS-p)])
        ROW(1,ID+POS-p) = 1;
    else
    end
    
    
    
    %% Connection 4
    
    if Lmin(R,C)>=0 && (Bmax(R+1,C)==1 || Lmax(R,C+1)==1)
        %disp(['Connection between node ',num2str(ID),' and node ',num2str(ID+POS-p+1)])
        ROW(1,ID+POS-p+1) = 1;
    else
    end
    
end

return
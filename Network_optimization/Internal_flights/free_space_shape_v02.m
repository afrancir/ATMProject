function [Lmin,Lmax,Bmin,Bmax,A_list,Node_ID] = free_space_shape_v02(P,Q,Re,eps)

% Building the 2 trajectories
p = length(P(:,1))-1;
q = length(Q(:,1))-1;

Lmin   = zeros(q,p+1);
Lmax   = zeros(q,p+1);
Bmin   = zeros(q+1,p);
Bmax   = zeros(q+1,p);

        
for i=1:p+1
    for j=1:q
        PC        = P(i,:);
        P1        = Q(j,:);
        P2        = Q(j+1,:);
        T         = gc_int(PC,P1,P2,Re,eps);
        Lmin(j,i) = T(1);
        Lmax(j,i) = T(2);
    end
end

for i=1:q+1
    for j=1:p
        PC        = Q(i,:);
        P1        = P(j,:);
        P2        = P(j+1,:);
        T         = gc_int(PC,P1,P2,Re,eps);
        Bmin(i,j) = T(1);
        Bmax(i,j) = T(2);
    end
end

% Getting inference matrix A

N            = (2*p+1)*(q+1)+(p+1)*q;
A_list       = cell(N,1);
Node_ID      = zeros(N,4);
Node_ID(N,1) = 1;
Node_ID(N,2) = q+1;
Node_ID(N,3) = p+1;
Node_ID(N,4) = N;

        for i=1:N
            
            BLOCK = ceil(i/(3*p+2));
            POS   = i-(BLOCK-1)*(3*p+2);
            
            % First row of the block, and not last element on the right
            if mod(POS,(3*p+2))<=2*p+1 && mod(POS,(3*p+2))>=1
                % Edge node (not last node on the right though)
                if mod(POS,2)~=0 && POS~=2*p+1
                    ROW          = edge_link(i,p,q,Lmin,Lmax,Bmin,Bmax);
                    Node_ID(i,1) = 1;
                    Node_ID(i,2) = ceil(i/(3*p+2));
                    Node_ID(i,3) = ceil(POS/2);
                    Node_ID(i,4) = i;
                % Last node on the right    
                elseif mod(POS,2)~=0 && POS==2*p+1
                    ROW          = edge_link(i,p,q,Lmin,Lmax,Bmin,Bmax);
                    Node_ID(i,1) = 1;
                    Node_ID(i,2) = ceil(i/(3*p+2));
                    Node_ID(i,3) = ceil(POS/2);
                    Node_ID(i,4) = i;
                % Internal nodes on a horizontal edge    
                else
                    ROW          = hor_node_link(i,p,q,Lmin,Lmax,Bmin,Bmax);
                    Node_ID(i,1) = 2;
                    Node_ID(i,2) = ceil(i/(3*p+2));
                    Node_ID(i,3) = ceil(POS/2);
                    Node_ID(i,4) = i;
                end
                
            elseif mod(POS,(3*p+2))>2*p+1
                ROW          = ver_node_link(i,p,q,Lmin,Lmax,Bmin,Bmax);
                Node_ID(i,1) = 3;
                Node_ID(i,2) = ceil(i/(3*p+2));
                Node_ID(i,3) = i-(BLOCK-1)*(3*p+2)-(2*p+1);
                Node_ID(i,4) = i;
             else
                ROW        = ver_node_link(i,p,q,Lmin,Lmax,Bmin,Bmax);
                Node_ID(i,1) = 3;
                Node_ID(i,2) = ceil(i/(3*p+2));
                Node_ID(i,3) = i-(BLOCK-1)*(3*p+2)-(2*p+1);
                Node_ID(i,4) = i;
            end
            
            if ~isempty(find(ROW,1))
                A_list{i,1} = find(ROW);
            else
            end
            
        end

return
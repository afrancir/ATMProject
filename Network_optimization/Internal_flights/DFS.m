function [is_solution,A] = DFS(startpoint,endpoint,G)

% Initializing the visited vector with all zeros. It means that initially
% all nodes have not been visite yet
visited = zeros(length(G),1);
% Initialization of the stack vector. The first element is the start point
stack   = startpoint;
% Initialization of the A matrix
A       = ones(3,1);
% Initialization of the level index
level   = 1;
% Initialization of the is+solution index
is_solution = 0;

while ~isempty(stack)
    v          = stack(end);
    stack(end) = [];
    
    % Current node has not been visited yet
    if visited(v)==0
        % Set this node as visited
        visited(v) = 1; 
        
        if v ~= startpoint
            level = A(3,A(1,:)==v);
        else
        end
        
        if v == endpoint
            is_solution = 1;
            break;
        else
        end
        
        for i=1:numel(G{v,1})
            stack = vertcat(stack,G{v,1}(i));
        end
        
        dummy      = zeros(3,numel(G{v,1}));
        dummy(1,:) = G{v,1};
        dummy(2,:) = v*ones(1,numel(G{v,1}));
        dummy(3,:) = (level+1)*ones(1,numel(G{v,1}));
        A          = horzcat(A,dummy);
    else
    end
end
    
    
    

return
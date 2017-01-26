function [FD,no_sol,sol] = FD_computation(ROUTE1,ROUTE2,r,alpha0,tolerance)

is_solution = 1;
no_sol      = [];
sol         = [];
while is_solution % keep iterating until a solution is found
    [~,~,~,~,A,~] = free_space_shape(ROUTE1,ROUTE2,r,alpha0);
    [d,~,~,~]     = dfs(A,1,[],numel(A(:,1)));
    
    
    if d(numel(A(:,1)),1)==-1 % no solution found
        no_sol = vertcat(no_sol,alpha0);
    else
        sol = vertcat(sol,alpha0);
    end
    
    if isempty(no_sol) % all solutions so far are overestimation
        alpha0 = 0.5*alpha0;
    elseif isempty(sol) % all solutions so far are underestimation
        alpha0 = 2.0*alpha0;
    else % we have at least one element in both arrays
        error = (min(sol)-max(no_sol))/max(no_sol)*100;
        if error<=tolerance
            is_solution = 0;
            FD = min(sol);
            break
        else
            alpha0 = 0.5*(max(no_sol)+min(sol));
        end
    end
end

return
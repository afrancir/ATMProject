function [FD,no_sol,sol] = FD_computation_v02(ROUTE1,ROUTE2,r,alpha0,tolerance)

found_FD    = 1;
no_sol      = [];
sol         = [];
while found_FD % keep iterating until a solution is found
    [~,~,~,~,A_list] = free_space_shape_v02(ROUTE1,ROUTE2,r,alpha0);
    startpoint       = 1;
    endpoint         = length(A_list);
    is_solution      = DFS(startpoint,endpoint,A_list);
    
    
    if is_solution == 0 % no solution found
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
            found_FD = 0;
            FD       = min(sol);
            break
        else
            alpha0 = 0.5*(max(no_sol)+min(sol));
        end
    end
end

return
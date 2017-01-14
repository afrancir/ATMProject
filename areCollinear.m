function logical = areCollinear(Matrix)
% function logical = areCollinear(Matrix)
%
% this function can be optimized, namely instead of "one point jumps" it
% could go from two to two, since if 3 points are collinear you can check
% the next iteration only with the third one and 2 new and if they ar enot
% collinear then the array is not collinear.
[row,col] = find(Matrix);
N = numel(row);
logical = 1;
for n = 1:N
    p1 = [row(n),col(n)];
    p2 = [row(n+1),col(n+1)];
    p3 = [row(n+2),col(n+2)];
    mat = [p1(1)-p3(1) p1(2)-p3(2); ...
    p2(1)-p3(1) p2(2)-p3(2)];
    if (det(mat) == 0)
        continue
    else
        logical = 0;
        return 
    end 
end 

end %function


function [smol] = divvy(tol,interval)
[rows,cols] = size(tol);
newrows = ceil(rows/interval);
smol = zeros(newrows,cols);
for n = 1:newrows
    upper = ((n-1)*interval) + 1;
    lower = n*interval;
    for m = 1:cols
        smol(n,m) = mean(tol(upper:lower, m)); 
    end
end
end


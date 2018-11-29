function [c,r] = chebycenter(A,b)
%CHEBYCENTER Compute Chebyshev center of polytope Ax <= b.
%  The Chebyshev center of a polytope is the center of the largest
%  hypersphere enclosed by the polytope. 
%  Requires optimization toolbox.
[n,p] = size(A);
an = sqrt(sum(A.^2,2));
A1 = zeros(n,p+1);
A1(:,1:p) = A;
A1(:,p+1) = an;
f = zeros(p+1,1);
f(p+1) = -1;
options = optimset;
% options = optimset(options,'Display', 'off');
options = optimoptions('linprog','Display', 'off','Algorithm','dual-simplex','OptimalityTolerance',1e-9,'ConstraintTolerance',1e-9,'MaxIterations',1e7);
c = linprog(f,A1,b,[],[],[],[],[],options);

if length(c)<p+1
    error('linprog');
end
if any(A1*c>b+1e-8)
    error('linprog error');
end
r = c(p+1);
c = c(1:p);
if any(A*c>b+1e-8)
    error('A*c>b');
end
end

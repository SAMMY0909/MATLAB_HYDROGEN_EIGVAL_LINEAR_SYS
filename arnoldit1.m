%Function sanity check
%Program to test that the arnoldi implementation in my code is perfectly
%alright and has got no problems of its own
a=[1 2 -1 0; 2 4 -2 1;3  8  -3 2;4 16 -4 3];
q1=[1;1;1;1];
%
%direct calculation of eigenvalues of a
eigs(a)
m=4;
%implement arnoldi and then calculate eigenvalues of a
[Q,H] = arnoldit(a,q1,m);

%line 13 should give the same result as line 7
eigs(H(1:m,1:m))
function [Q,H] = arnolditer(hess,q1,m)

%   ARNOLDI    Arnoldi iteration
%   [Q,H] = ARNOLDI(A,q1,m) has M Arnoldi iterations 
%    with N-by-N matrix A and initial vector q1.  If m < n it produces an n-by-(m+1) matrix Q with
%   orthonormal columns and an (m+1)-by-m upper Hessenberg matrix H such
%   that A*Q(:,1:M) = Q(:,1:M)*H(1:M,1:M) + H(M+1,M)*Q(:,M+1)*E_M', where
%   E_M is the M'th column of the M-by-M identity matrix.
%  Begin Arnoldi
t = length(hess);
if nargin < 3, m = t; end
q1 = q1/norm(q1);
Q = zeros(t,m); 
Q(:,1) = q1;
H = zeros(m+1,t);

for k=1:m
    z = hess*Q(:,k);
    for i=1:k
        H(i,k) = Q(:,i)'*z;
        z = z - H(i,k)*Q(:,i);
    end
    if k < t
       H(k+1,k) = norm(z);
       if H(k+1,k) == 0, return, end
       Q(:,k+1) = z/H(k+1,k);
   end
end
end
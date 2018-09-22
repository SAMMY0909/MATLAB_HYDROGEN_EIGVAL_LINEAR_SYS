%The radial equation of hydrogen is discretised with 'n' number of
%gridpoints with grid spacing 'd' and it is converetd to a linear system of
%equations, and hence the resulting matrix is a tridiagonal(band-diagonal)
%matrix. The matrix is made to undergo Arnoldi iterations and the smallest
%eigenvalues are sorted and displayed both on screen and dumpedin a text
%file. The energy units are taken to be in electron volts (eV) and the
%units of the radial variable is put in picometers. 38178 and 1440 are the
%respective conversion factors wherever displayed.

n=1000;%maximum number of iterations
d=1; % grid spacing
l=0; % angular momentum, we can choose l=1,2 etc
hess=zeros(n,n);% pre allocate matrix with zeroes
h12=-38178.0/(d^2.0);% constant elements of band diagonal

%1,1 entry of the band diagonal discretisation matrix
hess(1,1)=+38178.0/(d^2.0)+38178.0*l*(l+1.0)/((1.0*d)^2.0)-(1438.0/d);
hess(1,2)=h12;

%loop to fill in the matrix cells with specified values according to the
%linear system
for i = 2:n
    hess(i,i)=+(38178.0*2.0)/(d^2.0)+38178.0*l*(l+1.0)/((i*d)^2.0)-(1438.0)/(i*d);
    hess(i,i-1)=h12;
    if i<n
        hess(i,i+1)=h12;
    end
end

%save the band diagonal matrix WARNING: USE YOUR FILE PATH
save '/home/soumyananda/Dropbox/565_COMPUTATION_PROJECT/MY_FINAL_565_UOFL_PROJECT/hessdat.txt' hess -ASCII

%for arnoldi iteration
 m=n; %no of iterations of arnoldi matrix
 q1=zeros(n,1); % initialise inital vector q1
 q1(1)=1.0;
 q1(2)=1.0;
 [Q,H]=arnolditer(hess,q1,m);% implement Arnoldi iteration
 
 %e1=eigs(hess); %directly use the eigs function of matlab on tridiagonal matrix
 %however, it costs more memory operations and is inefficient in terms of
 %number of operations
 
 
 %e1 %display eigenvalues inline in prompt
 
 %save Upper diagonal matrix obtained after arnoldi, WARNING: USE YOUR FILE PATH
 save '/home/soumyananda/Dropbox/565_COMPUTATION_PROJECT/MY_FINAL_565_UOFL_PROJECT/Hmat.txt' H  -ASCII
 
 %store the lowest magnitude eigenvalues in array e, +ve values DISCARDED,
 %-ve vales accepted. Code displays both. Now eigs function can be used
 % since the eigenvalues of the hessenberg matrix is eigenvalue of the tridiagonal matrix
 % by virtue of linear algebra. But eigenvalue calculation is more
 % efficient here and faster.
 e=eigs(H(1:m,1:m),10,'sm');
 e % display eigenvalues inline in prompt
save '/home/soumyananda/Dropbox/565_COMPUTATION_PROJECT/MY_FINAL_565_UOFL_PROJECT/hyd_eig.txt' e  -ASCII
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

clc;
clear;

% Size of the domain
L = 128;

% Max time
T = 400;

% Number of intervals in the x direction
n = 1024;

% Number of intervals in time (A small solution leads to instabilities
% related to the implicit - explicit scheme)
Nt = 2400;

% The mesh size h
h = 2 * L / n;
x = linspace(-L, L, n+1);

% Time
dt = T / (Nt);
t = linspace(0, T, Nt+1);

% One dimensional discrete first derivative (Central Diff)
e = ones(n,1);
F = spdiags([-e e],[-1 1],n,n) / (2 * h);
% Periodic boundary
F(1,n) = -1;
F(n,1) = 1;

% One dimensional discrete second derivative (\Delta)
D = - gallery("tridiag", n);
% Periodic Boundary
D(1,n) = 1; D(n,1) = 1;
D = D / (h^2);

% One dimensional discrete forth derivative (\Delta^2)
D2 = D'*D;

% Identity matrix
I = speye(n);

% Forming the (CN Implicit Linear part - Explicit Nonlinear part)
A = I+dt/2*D+dt/2*D2;
B = I-dt/2*D-dt/2*D2;

% Matrix to retain solution at each timestep
S = zeros(n,Nt+1);

% Initial condition
S(:,1) = cos(x(1:end-1)') + 0.15*cos(x(1:end-1)'/8).*(1+2*sin(x(1:end-1)'/8));

% Initialize the contribution of non-linear parts
Nn=(F*S(:,1)).*S(:,1);

% Iterative solution
for i=1:Nt
    No=Nn;
    Nn=(F*S(:,i)).*S(:,i);
    S(:,i+1)=A\(B*S(:,i)-1.5*dt*Nn+0.5*dt*No);
end

% Plot time evolution and account for periodic boundary
imagesc([S;S(1,:)]);colormap autumn;axis off;
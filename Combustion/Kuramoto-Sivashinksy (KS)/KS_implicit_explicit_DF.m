% Solution of Kuramoto - Sivashinsky equation
% In 1D periodic domain, using Finite Differences
% and Simple IMEX (Implicit - Explicit) integration scheme

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

% Forming the (Fully Implicit Linear part - Explicit Nonlinear part)
A = I+dt*D+dt*D2;

% Matrix to retain solution at each timestep
S = zeros(n,Nt+1);

% Initial condition
S(:,1) = cos(x(1:end-1)') + 0.15*cos(x(1:end-1)'/8).*(1+2*sin(x(1:end-1)'/8));

% Iterative solution
for i=1:Nt
    S(:,i+1)=A\(S(:,i)-dt*(F*S(:,i)).*S(:,i));
end

% Plot time evolution and account for periodic boundary
imagesc([S;S(1,:)]);colormap autumn;axis off;
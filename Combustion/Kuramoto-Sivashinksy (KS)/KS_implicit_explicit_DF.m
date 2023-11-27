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
h = L / n;
x = linspace(0, L, n+1);

% Time
dt = T / (Nt);
t = linspace(0, L, Nt+1);

% One dimensional discrete first derivative (Central Diff)
e = ones(n,1);
F = spdiags([-e e],[-1 1],n,n) / (2 * h);
% Periodic boundary
F(1,n) = -1;
F(n,1) = 1;

% One dimensional discrete second derivative (\Delta)
T = - gallery("tridiag", n);
% Periodic Boundary
T(1,n) = 1; T(n,1) = 1;
T = T / (h^2);

% One dimensional discrete forth derivative (\Delta^2)
T2 = T'*T;

% Identity matrix
I = speye(n);

% Forming the (Fully Implicit Linear part - Explicit Nonlinear part)
A = I+dt*T+dt*T2;

% Matrix to retain solution at each timestep
S = zeros(n,Nt+1);

% Initial condition
S(:,1) = cos(x(1:end-1)') + 0.1*cos(x(1:end-1)'/16).*(1+2*sin(x(1:end-1)'/16));

% Iterative solution
for i=1:Nt
    S(:,i+1)=A\(S(:,i)-dt*(F*S(:,i)).*S(:,i));
end

% Plot time evolution and account for periodic boundary
imagesc([S;S(1,:)]);colormap autumn;axis off;xlabel('time');ylabel('space');
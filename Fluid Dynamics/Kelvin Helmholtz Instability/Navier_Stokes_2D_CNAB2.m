% Incompressible Navier - Stokes 2D
% on doubly periodic domain
% Rayleigh Taylor instability
% using FFT and CNAB2

% Select plot type
% pltype = 0 -> Stream function
% pltype = 1 -> Tracer
pltype = 1;

% Size of domain
Lx = 2;
Ly = 2;

% Number of points for discretization
nx = 256;
ny = 256;
dx = Lx / nx;
dy = Ly / ny;

% Reynolds number
Re = 1e4;

% Maximum time
T = 4;

% Number of time steps
Nt = 3000;

% Frequencies
klinx = (2.0 * pi / Lx) .* [0:nx/2,-nx/2+1:-1];
kliny = (2.0 * pi / Ly) .* [0:ny/2,-ny/2+1:-1];
[kx, ky] = meshgrid(klinx, kliny);

% Space discretization
x = linspace(0,Lx,nx+1);
x(1)=[];
y = linspace(0,Ly,ny+1);
y(1)=[];
[x,y] = meshgrid(x,y);

% Time step
dt = T / Nt;

% 1st derivative
Dx = 1i * kx;
Dy = 1i * ky;
D = sqrt(Dx.^2+Dy.^2);

% Crank - Nicolson matrices
A = dt^-1 - 0.5/Re * D.^2;
B = dt^-1 + 0.5/Re * D.^2;
IA = A.^-1;

% Initialize matrices
% Velocities across x
Q = (0.5 + 0.5 * tanh(10-20*abs(1-2*y/Ly)));
um1 = fftn((Q).*(1+0.5*sin(Lx*pi*x)));
us = um1;
u = um1;
% Velocities across y
vm1 = fftn(zeros(size(D)));
vs = vm1;
v = vm1;

% Compute roots of unity
M = 32;
r = exp ( 1i * pi * ( (1:M) - 0.5 ) / M );

% Form matrix to compute the points required for the contour integrals
r = permute(reshape(repmat(r,ny*nx,1),nx,ny,[]),[2,1,3]);

% Fix D.^2 matrix through modified Cauchy integral equation
D2 = (mean(D.^2 + r,3)); 

% Time advancing loop
for i=1:Nt
    % Tracer equation update
    Q = Q - (tracer(Q, dt / dy * real(ifftn(v))) + tracer(Q', dt / dx * real(ifftn(u))')');
    
    % x - Direction
    Hxn = H(u,v,Dx,Dy);
    Hxn1 = H(um1,vm1,Dx,Dy);
    
    % y - Direction
    Hyn = H(v,u,Dy,Dx);
    Hyn1 = H(vm1,um1,Dy,Dx);
    
    % Compute intermediate solution
    us = u - us + IA .* (-1.5*Hxn+0.5*Hxn1+B.*u);
    vs = v - vs + IA .* (-1.5*Hyn+0.5*Hyn1+B.*v);
    
    % Retain solution of previous step
    um1 = u;
    vm1 = v;
    
    % Enforce conservation of mass
    phi = (Dx.*us + Dy.*vs) ./ D2;
    u = us - Dx .* phi;
    v = vs - Dy .* phi;
    
    % Plot result and drawnow
    if pltype == 0
        Q = -real(ifftn((Dy.*u-Dx.*v)./D2));
    end
    contour(x,y,Q);colormap hot;xlabel('x');ylabel('y');
    drawnow
end

% Convective terms
% Input:
%   u:   Velocities across x
%   v:   Velocities across y
%   Dx:  Spectral differenciation operator across x
%   Dy:  Spectral differenciation operator across y
% Output:
%   out: Value of the nonlinear function
function [out] = H(u,v,Dx,Dy)
    out = fftn(ifftn(u).*ifftn(Dx.*u) + ifftn(v).*ifftn(Dy.*u));
end

% Tracer equation
% Inputs:
%   Q:  Solution to the tracer equation
%   c:  Velocity vector
% Outputs:
%   dQ: Change in the solution with respect to direction of c
function [dQ] = tracer(Q,c)
    dQ = c .* (Q-circshift(Q,1)) - 0.5 * c .* (1 - c) .* ((Q-circshift(Q,1)) + (Q-circshift(Q,-1)));
end

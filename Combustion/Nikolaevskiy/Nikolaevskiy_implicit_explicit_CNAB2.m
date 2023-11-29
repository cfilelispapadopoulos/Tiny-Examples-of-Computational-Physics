% Solution of the Nikolaevskiy equation
% In 2D doubly periodic domain, using FFT
% and CNAB2 integration scheme

clc;
clear;

% Size of the domain
Lx = 80*pi;
Ly = 30*pi;

% Parameter r of the Nikolaevskiy equation
r = 0.3;

% Max time
T = 20;

% Number of intervals in the x direction
nx = 512;
ny = 512;

% Number of intervals in time (A small solution leads to instabilities
% related to the implicit - explicit scheme)
Nt = 1024;

% The mesh size h
klinx = (2.0 * pi / Lx) .* [0:nx/2,-nx/2+1:-1];
kliny = (2.0 * pi / Ly) .* [0:ny/2,-ny/2+1:-1];
[kx, ky] = meshgrid(klinx, kliny);

% Time
dt = T / (Nt);
t = linspace(0, T, Nt+1);

% 1st/2nd/4th Derivative Matrices
Dx = 1i * kx;
Dy = 1i * ky;
D = sqrt(Dx.^2 + Dy.^2);

% Linear Operator
LL = -((1 - r) * D.^2 + 2*D.^4 + D.^6);

% Matrix to retain solution - Initial condition
% Starting with rand gives rise to interesting
% Dynamics
Si = randn(ny,nx);

% Crank-Nicolson matrices and Corresponding Fourier
% Transforms
A = ones(ny,nx) + 0.5*dt*LL;
B = ones(ny,nx) - 0.5*dt*LL;
IA = A.^-1;
Si = fftn(Si);
Nn = -0.5*fftn(abs(ifftn(Dx.*Si).^2+(ifftn(Dy.*Si)).^2));

% Iterative solution
for i=1:Nt
    % Retain previous value
    No=Nn;
    % Compute Nonlinear part
    Nn = -0.5*fftn(abs(ifftn(Dx.*Si).^2+(ifftn(Dy.*Si)).^2));
    % Form right hand side
    rhs=(B.*Si+1.5*dt*Nn-0.5*dt*No);
    % Solve
    Si=IA.*rhs;
    
    % Plot time evolution and account for periodic boundary
    % Negate and rescale for aesthetics
    imagesc(imresize(-real(ifftn(Si)),6));colormap hot;pbaspect([2*Lx/max(Lx,Ly) 2*Ly/max(Lx,Ly) 1]);axis off;
    M = getframe;
end


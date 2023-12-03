% Solution of the Cahn-Hiliard
% In 2D doubly periodic domain, using FFT
% and ETDRK4 integration scheme

clc;
clear;

% Size of the domain
Lx = 8;
Ly = 8;

% Parameter r of the Cahn-Hiliard equation
epsilon2 = 0.01;
alpha = 1;

% Max time
T = 8;

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
LL = - epsilon2*D.^4+alpha*D.^2;

% Matrix to retain solution - Initial condition
% Starting with rand gives rise to interesting
% Dynamics
Si = -1+2*rand(ny,nx);

% ETD RK4 required quantities
E = exp ( dt * LL );
E2 = exp ( dt * LL / 2.0 );
M = 16;

% Compute roots of unity
r = exp ( 1i * pi * ( (1:M) - 0.5 ) / M );

% Fourier transform of initial condition
v = fftn(Si);

% Form matrix to compute the points required for the contour integrals
r = permute(reshape(repmat(r,ny*nx,1),nx,ny,[]),[2,1,3]);

% Linear operator multiplied by time step and integrated with
% respect to unit circle
LM = dt * LL + r;

% Estimate all coefficients for RK4 and the 
% matrix Q required for the RK4 terms
a = zeros(ny,nx);
fu = dt * real(mean((-4.0-LM+exp(LM).*(4.0-3.0*LM+LM.^2))./LM.^3,3));
fab = dt * real(mean((2.0+LM+exp(LM).*(-2.0+LM))./ LM.^3,3));
fc = dt * real(mean((-4.0-3.0*LM-LM.^2+exp(LM).*(4.0-LM))./LM.^3,3));
Q  = dt * real(mean((exp(LM/2.0)-1.0)./LM,3));

% Time stepping
for i = 1 : Nt
    % Compute the RK4 terms for different steps corresponding
    % to the nonlinear part of the Differential Equation
    Nv = (D.^2).*fftn((ifftn(v).^3-(1+alpha)*ifftn(v)));
    a = E2 .* v + Q .* Nv;
    Na = (D.^2).*fftn((ifftn(a).^3-(1+alpha)*ifftn(a)));
    b = E2 .* v + Q .* Na;
    Nb = (D.^2).*fftn((ifftn(b).^3-(1+alpha)*ifftn(b)));
    c = E2 .* a + Q .* ( 2.0 * Nb - Nv );
    Nc = (D.^2).*fftn((ifftn(c).^3-(1+alpha)*ifftn(c)));
    
    % Compute solution at next timestep
    v = E .* v + fu.*Nv + 2.0 * fab.*(Na+Nb) + fc.*Nc;
    
    % Plot time evolution and account for periodic boundary
    % Rescale for aesthetics
    imagesc(imresize(-real(ifftn(v)),4));colormap hot;
    pbaspect([2*Lx/max(Lx,Ly) 2*Ly/max(Lx,Ly) 1]);axis off;
    % Getframe is used to smooth out plotting
    MM = getframe;
 end


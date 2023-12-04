% Solution of the Cahn-Hiliard
% In 3D triply periodic domain, using FFT
% and ETDRK4 integration scheme

clc;
clear;

% Choose slice (in the z axis to print)
sl = 64;

% Size of the domain
Lx = 8;
Ly = 8;
Lz = 8;

% Parameter r of the Cahn-Hiliard equation
epsilon2 = 0.01;
alpha = 1;

% Max time
T = 8;

% Number of intervals in the x direction
nx = 128;
ny = 128;
nz = 128;

% Number of intervals in time (A small solution leads to instabilities
% related to the implicit - explicit scheme)
Nt = 1024;

% The mesh size h
klinx = (2.0 * pi / Lx) .* [0:nx/2,-nx/2+1:-1];
kliny = (2.0 * pi / Ly) .* [0:ny/2,-ny/2+1:-1];
klinz = (2.0 * pi / Lz) .* [0:nz/2,-nz/2+1:-1];
[kx, ky, kz] = meshgrid(klinx, kliny, klinz);

% Time
dt = T / (Nt);
t = linspace(0, T, Nt+1);

% 1st/2nd/4th Derivative Matrices
Dx = 1i * kx;
Dy = 1i * ky;
Dz = 1i * kz;
D = sqrt(Dx.^2 + Dy.^2 + Dz.^2);

% Linear Operator
LL = - epsilon2*D.^4+alpha*D.^2;

% Matrix to retain solution - Initial condition
% Starting with rand gives rise to interesting
% Dynamics
Si = -1+2*rand(ny,nx,nz);

% ETD RK4 required quantities
E = exp ( dt * LL );
E2 = exp ( dt * LL / 2.0 );
M = 32;

% Compute roots of unity
r = exp ( 1i * pi * ( (1:M) - 0.5 ) / M );

% Fourier transform of initial condition
v = fftn(Si);

% Form matrix to compute the points required for the contour integrals
r = permute(reshape(repmat(r,ny*nx*nz,1),nx,ny,nz,[]),[2,1,3,4]);

% Linear operator multiplied by time step and integrated with
% respect to unit circle
LM = dt * LL + r;

% Estimate all coefficients for RK4 and the 
% matrix Q required for the RK4 terms
a = zeros(ny,nx,nz);
fu = dt * real(mean((-4.0-LM+exp(LM).*(4.0-3.0*LM+LM.^2))./LM.^3,4));
fab = dt * real(mean((2.0+LM+exp(LM).*(-2.0+LM))./ LM.^3,4));
fc = dt * real(mean((-4.0-3.0*LM-LM.^2+exp(LM).*(4.0-LM))./LM.^3,4));
Q  = dt * real(mean((exp(LM/2.0)-1.0)./LM,4));

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
    S = ifftn(v);
    imagesc(imresize(real(S(:,:,sl)),3));colormap hot;
    pbaspect([2*Lx/max(Lx,Ly) 2*Ly/max(Lx,Ly) 1]);axis off;
    
    % Getframe is used to smooth out plotting
    MM = getframe;
 end


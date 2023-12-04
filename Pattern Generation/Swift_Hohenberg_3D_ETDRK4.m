% Solution of the Swift-Hohenberg
% In 3D triply periodic domain, using FFT
% and ETDRK4 integration scheme

clc;
clear;

% Choose slice along the z axis to plot
sl = 128;

% Size of the domain
Lx = 50;
Ly = 50;
Lz = 50;

% Parameter r of the Cahn-Hiliard equation
epsilon = 1;

% Max time
T = 100;

% Number of intervals in the x direction
nx = 256;
ny = 256;
nz = 256;

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
LL = -(D.^2+1).^2;

% Matrix to retain solution - Initial condition
% Starting with rand gives rise to interesting
% Dynamics
Si = rand(ny,nx,nz);

% ETD RK4 required quantities
E = exp ( dt * LL );
E2 = exp ( dt * LL / 2.0 );
M = 16;

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
a = zeros(ny,nx);
fu = dt * real(mean((-4.0-LM+exp(LM).*(4.0-3.0*LM+LM.^2))./LM.^3,4));
fab = dt * real(mean((2.0+LM+exp(LM).*(-2.0+LM))./ LM.^3,4));
fc = dt * real(mean((-4.0-3.0*LM-LM.^2+exp(LM).*(4.0-LM))./LM.^3,4));
Q  = dt * real(mean((exp(LM/2.0)-1.0)./LM,4));

% Time stepping
for i = 1 : Nt
    % Compute the RK4 terms for different steps corresponding
    % to the nonlinear part of the Differential Equation
    Nv = epsilon*v-fftn(ifftn(v).^5);
    a = E2 .* v + Q .* Nv;
    Na = epsilon*a-fftn(ifftn(a).^5);
    b = E2 .* v + Q .* Na;
    Nb = epsilon*b-fftn(ifftn(b).^5);
    c = E2 .* a + Q .* ( 2.0 * Nb - Nv );
    Nc = epsilon*c-fftn(ifftn(c).^5);
    
    % Compute solution at next timestep
    v = E .* v + fu.*Nv + 2.0 * fab.*(Na+Nb) + fc.*Nc;
    
    % Plot time evolution and account for periodic boundary
    % Rescale for aesthetics
    S = ifftn(v);
    imagesc(imresize(-real(S(:,:,sl)),4));colormap hot;
    pbaspect([2*Lx/max(Lx,Ly) 2*Ly/max(Lx,Ly) 1]);axis off;
    % Getframe is used to smooth out plotting
    MM = getframe;
 end


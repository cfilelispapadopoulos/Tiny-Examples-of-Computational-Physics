% Solution of the Navier - Stokes equation (Velocity - Vorticity)
% In 2D doubly periodic domain, using FFT
% and Crank - Nicolson / Adam - Bashforth
% integration scheme

clc;
clear;

% Size of the domain
Lx = 2*pi;
Ly = 2*pi;

% Reynolds number
Re = 1000;

% Max time
T = 1000;

% Number of intervals in the x direction
nx = 128;
ny = 128;

% Number of intervals in time (A small solution leads to instabilities
% related to the implicit - explicit scheme)
Nt = 10000;

% Create frequencies
klinx = (2.0 * pi / Lx) .* [0:nx/2,-nx/2+1:-1];
kliny = (2.0 * pi / Ly) .* [0:ny/2,-ny/2+1:-1];
[kx, ky] = meshgrid(klinx, kliny);

% Physical space coordinates
x = linspace(-Lx,Lx,nx+1);
% Remove duplicate
x(end)=[];
y = linspace(-Ly,Ly,ny+1);
% Remove duplicate
y(end)=[];
[xx, yy] = meshgrid(x,y); 

% Time
dt = T / (Nt);
t = linspace(0, T, Nt+1);

% 1st/2nd Derivative Matrices
Dx = 1i * kx;
Dy = 1i * ky;
D = sqrt(Dx.^2 + Dy.^2);

% Linear Operator
LL = D.^2;

% Matrix to retain solution - Initial condition
% Starting with rand gives rise to interesting
% Dynamics
w = -1+2*rand(ny,nx);

% Number of points for Cauchy integral formula
M = 32;

% Compute roots of unity
r = 0.5*exp ( 1i * pi * ( (1:M) - 0.5 ) / M );

% Form matrix to compute the points required for the contour integrals
r = permute(reshape(repmat(r,ny*nx,1),nx,ny,[]),[2,1,3]);

% Linear operator multiplied by time step and integrated with
% respect to unit circle
LM = mean(LL + r,3);
ILM = LM.^-1;

% Crank - Nicolson Matrices 
A = 1 / dt - 0.5 / Re * LL;
B = 1 / dt + 0.5 / Re * LL;
IA = A.^-1;

% Compute Stream Function
W = fft2(w);

% Compute Stream Function
Psi = - ILM .* W;

% Compute velocities from Psi
ux = real(ifft2(Dy .* Psi));
uy = real(ifft2(- Dx .* Psi));
wx = real(ifft2(Dx .* W));
wy = real(ifft2(Dy .* W));

% Compute nonlinearity
Nn = - (fft2(ux.*wx+uy.*wy));
Nn1 = Nn;

% Time stepping
for i = 1 : Nt    
    % Update vorticity
    W = IA .* (B .* W + 3/2 * Nn - 1/2 * Nn1);
    
    % Compute Stream Function
    Psi = - ILM .* W;

    % Compute velocities from Psi
    ux = real(ifft2(Dy .* Psi));
    uy = real(ifft2(- Dx .* Psi));
    wx = real(ifft2(Dx .* W));
    wy = real(ifft2(Dy .* W));
    
    Nn1 = Nn;
    
    % Compute nonlinearity
    Nn = - (fft2(ux.*wx+uy.*wy));
    
    % Plot time evolution and account for periodic boundary
    contour(xx,yy,real(ifftn(W)),50); colormap('hot');title(['Time: ' num2str(i*dt) ' sec'])
    drawnow;
 end


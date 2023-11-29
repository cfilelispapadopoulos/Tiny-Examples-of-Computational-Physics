% Solution of Kuramoto - Sivashinsky equation
% In 2D doubly periodic domain, using FFT
% and CNAB2 integration scheme

clc;
clear;

% Size of the domain
L = 80*pi;

% Max time
T = 200;

% Number of intervals in the x direction
n = 512;

% Number of intervals in time (A small solution leads to instabilities
% related to the implicit - explicit scheme)
Nt = 3*1024;

% The mesh size h
klin = (2.0 * pi / L) .* [0:n/2,-n/2+1:-1];
[kx, ky] = meshgrid(klin, klin);

% Time
dt = T / (Nt);
t = linspace(0, T, Nt+1);

% 1st/2nd/4th Derivative Matrices
Dx = 1i * kx;
Dy = 1i * ky;
D = sqrt(Dx.^2 + Dy.^2);
% Linear Operator
L =  D.^2 + D.^4;

% Matrix to retain solution - Initial condition
% Starting with rand gives rise to interesting
% Dynamics
Si = randn(n,n);

% Crank-Nicolson matrices and Corresponding Fourier
% Transforms
A = ones(n,n) + 0.5*dt*L;
B = ones(n,n) - 0.5*dt*L;
IA = A.^-1;
Si = fftn(Si);
Nn =  -fftn(0.5*abs(ifftn(Dx.*Si).^2+(ifftn(Dy.*Si)).^2));

% Iterative solution
for i=1:Nt
    % Retain previous value
    No=Nn;
    % Compute Nonlinear part
    Nn = -fftn(0.5*abs(ifftn(Dx.*Si).^2+(ifftn(Dy.*Si)).^2));
    % Form right hand side
    rhs=(B.*Si+1.5*dt*Nn-0.5*dt*No);
    % Solve
    Si=IA.*rhs;
    
    % Plot time evolution and account for periodic boundary
    % Negate and rescale for aesthetics
    imagesc(imresize(-real(ifftn(Si)),4));colormap hot;axis off;
    M = getframe;
end


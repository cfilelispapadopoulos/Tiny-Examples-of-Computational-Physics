% Simple visualization of waves

clear;
clc;

% Constants
% All units are SI
% Amplitude of wave
A = 0.2;
% Phase velocity
vp = 10;
% Phase of waves
phi = 0;
% Frequency
f = 0.5;
% Wavelength
lambda = vp / f;
% Wavenumber
k = 2*pi / lambda;
% Length of medium and sampling step
L = 50;
dl = 0.5;
% Maximum simulation time and time step
Tmax = 50;
dt = 0.01;
% Form all medium positions
x = 0:dl:L;

% Initialize time
t = 0;

% Simulation start
while (t<=Tmax)
% Compute position with respect to x and t
y=A*cos(k*(x-vp*t)+phi);

% Plot results with stem to be visibly more pleasing
stem(x,y);grid;xlabel('x');ylabel('y');title(['t = ' num2str(t)]);
axis([0, L, -2*A, 2*A]);

% For visualization purposes
MM = getframe;

% Update time variable
t = t + dt;
end

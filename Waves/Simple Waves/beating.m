% Wave beating composed of two waves
% With small difference in frequencies
% and same amplitude

clear;
clc;

% Constants
% All units are SI
% Amplitude of wave
A = 0.2;
% Phase velocities
vp_1 = 5;
vp_2 = 5;
% Phase of waves
phi_1 = pi/2;
phi_2 = 0;
% Frequencies
f_1 = 1;
f_2 = 0.9;
% Wavelengths
lambda_1 = vp_1 / f_1;
lambda_2 = vp_2 / f_2;
% Wavenumbers
k_1 = 2*pi / lambda_1;
k_2 = 2*pi / lambda_2;
% Length of medium and sampling step
L = 200;
dl = 0.5;
% Maximum simulation time and time step
Tmax = 10;
dt = 0.01;
% Form all medium positions
x = 0:dl:L;

% Initialize time
t = 0;

% Simulation start
while (t<=Tmax)
% Compute position with respect to x and t
y=A*cos(k_1*(x-vp_1*t)+phi_1)+A*cos(k_2*(x-vp_2*t)+phi_2);

% Plot results with stem to be visibly more pleasing
plot(x,y);grid;xlabel('x');ylabel('y');title(['t = ',num2str(t),', v_p = ',...
 num2str(0.5*(lambda_1 * f_1 + lambda_2 * f_2)), ', v_g = ' num2str(2*pi*(f_2-f_1)/(k_2-k_1))]);
hold;
% Plot envelope in red
ye=2 * A *cos(0.5*k_1*(x-vp_1*t)-0.5*k_2*(x-vp_2*t)+0.5*(phi_1-phi_2));
plot(x,ye,'--r');
plot(x,-ye,'--r');
hold off
% Fix axis for aesthetic reasons
axis([0, L, -4*A, 4*A]);

% For visualization purposes
MM = getframe;

% Update time variable
t = t + dt;
end

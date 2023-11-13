% Wave beating composed of two waves
% With small difference in frequencies
% and same amplitude

clear;
% clc;

% Constants
% All units are SI
% Number of waves
N = 2 * 4;
% Amplitude of waves
A = 0.2 * ones(N,1);
% Phase velocities
vp = 5 * ones(N,1);
% Max, min frequencies and step
min_f = 0.8;
max_f = 1.2;
f = linspace(min_f,max_f,N)';
% Wavelengths
lambda = vp ./ f;
% Wavenumbers
k = 2*pi ./ lambda;
% Difference between 2 consecutive wavenumbers 
dk = k(2) - k(1);
% Difference between 2 angular frequencies 
domega = 2*pi*(f(2) - f(1));
% Length of medium and sampling step
L = 200;
dl = 0.5;
% Maximum simulation time and time step
Tmax = 50;
dt = 0.01;
% Form all medium positions
x = repmat(0:dl:L,N,1);

% Initialize time
t = 0;

% Simulation start
while (t<=Tmax)
% Compute position with respect to x and t
y=sum(A.*cos(k.*(x-vp.*t)),1);
% Plot results with stem to be visibly more pleasing
% Along with phase, group velocities
plot(0:dl:L,y);grid;xlabel('x');ylabel('y');title(['t = ',num2str(t),', v_p = ',...
 num2str(1/N*sum(lambda.*f)), ', v_g = ' num2str(domega/dk)]);
hold;
% Plot envelope in red
% The equation for the envelope is
% A sin(0.5 N (dk x - domega t) / sin (0.5 (dk x - domega t)
ye=A *sin(0.5*N*(dk*(0:dl:L)-domega*t))./(sin(0.5 *(dk*(0:dl:L)-domega*t)));
plot(0:dl:L,ye,'--r');
plot(0:dl:L,-ye,'--r');
hold off
% Fix axis for aesthetic reasons
axis([0, L, -2*sum(A), 2*sum(A)]);

% For visualization purposes
MM = getframe;

% Update time variable
t = t + dt;
end

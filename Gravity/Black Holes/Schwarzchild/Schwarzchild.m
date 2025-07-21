clear;
clc;

% Parameters
M = 1;                    % Black hole mass
r_init = 50;              % Starting radius
phi_init = -pi/4;         % Initial angle
num_rays = 20;            % Number of rays
b_values = linspace(3.5, 15, num_rays);  % Impact parameters

% ODE options
options = odeset('RelTol',1e-9,'AbsTol',1e-12,'Events',@event_horizon);

% Setup plots and plot event horizon
figure; hold on; axis equal
xlabel('x'); ylabel('y')
title('Light rays bending near Schwarzschild black hole')
xlim([-r_init r_init]); ylim([-r_init r_init])
plot_circle(0,0,2*M,'k') % Event horizon

for i=1:num_rays
    b = b_values(i);

    E = 1;           % Energy per unit mass
    L = b * E;       % Angular momentum per unit mass
    f = 1 - 2*M/r_init;

    % Solve for dt/dlambda from null condition:
    % 0 = -f*(dt/dλ)^2 + (1/f)*(dr/dλ)^2 + r^2*(dφ/dλ)^2
    dphidlambda = b / r_init^2;
    drdlambda = -sqrt(E^2 - f * (L^2)/(r_init^2));
    dtdlambda = 0;%E / f;

    % Initial state vector: [t, r, phi, dt/dλ, dr/dλ, dφ/dλ]
    y0 = [0; r_init; phi_init; dtdlambda; drdlambda; dphidlambda];

    lambda_span = [0 2000];
    [lambda, y] = ode15s(@(lambda,y) geodesic_eqns_full(lambda,y,M), lambda_span, y0, options);

    r = y(:,2);
    phi = y(:,3);

    % Cartesian coords
    [x, y_cart] = pol2cart(phi,r);

    plot(x, y_cart, 'LineWidth', 1.5, 'DisplayName',sprintf('b=%.2f',b));
end

legend();
hold off;

% Full geodesic system including time
function dydlambda = geodesic_eqns_full(~, y, M)
    t = y(1);
    r = y(2);
    phi = y(3);
    dtdl = y(4);
    drdl = y(5);
    dphidl = y(6);

    f = 1 - 2*M/r;

    % Geodesic equations (second derivatives)
    d2tdl2   = - (2*M / (r^2 * f)) * dtdl * drdl;
    d2rdl2   = r * f * dphidl^2 - (M / (r^2 * f)) * drdl^2 + f * M / r^2 * dtdl^2;
    d2phidl2 = - (2 / r) * drdl * dphidl;

    dydlambda = zeros(6,1);
    dydlambda(1) = dtdl;
    dydlambda(2) = drdl;
    dydlambda(3) = dphidl;
    dydlambda(4) = d2tdl2;
    dydlambda(5) = d2rdl2;
    dydlambda(6) = d2phidl2;
end

% Stop at event horizon (r = 2M)
function [value,isterminal,direction] = event_horizon(~, y)
    r = y(2);
    value = r - 2;
    isterminal = 1;
    direction = -1;
end

% Draw event horizon
function plot_circle(xc,yc,R,color)
    theta = linspace(0,2*pi,200);
    x = xc + R*cos(theta);
    y = yc + R*sin(theta);
    plot(x,y,color,'LineWidth',2,'DisplayName','Event Horizon');
end

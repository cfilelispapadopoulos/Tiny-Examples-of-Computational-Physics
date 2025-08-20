clear; clc;

% Black hole parameters
M = 1;                                 % Mass of black hole
a = 0.95;                              % Angular momentum parameter
r0 = 20;                               % Initial radius
phi0 = 0;                              % Initial phi
theta0 = pi/2 - 0.1;                   % Slightly off equator
E = 1;                                 % Photon energy (fixed)
num_rays = 10;
L_values = linspace(2, 10, num_rays);  % Angular momenta
dtMax = 6;                             % Normalized time for visualization

opts = odeset('RelTol',1e-10,'AbsTol',1e-12);

figure; hold on; grid on; axis equal; view(3)
xlabel('x'); ylabel('y'); zlabel('z')
title('3D Photon Trajectories near Kerr BH')


for theta0 = pi/2+linspace(-0.5,0.5,10)
    for i = 1:num_rays
        L = L_values(i);
        p_t = -E;
        p_phi = L;

        r = r0;
        theta = theta0;

        ginv = inv_metric_Kerr_full(r, theta, M, a);
        gtt = ginv(1,1); gtphi = ginv(1,4); gphiphi = ginv(4,4);
        grr = ginv(2,2); gthetatheta = ginv(3,3);

        % Assume p_theta = 0 initially
        p_theta = 0;

        % Solve for pr using null condition
        A = gtt * p_t^2 + 2*gtphi*p_t*p_phi + gphiphi*p_phi^2;
        pr = -sqrt( -A / grr );

        Y0 = [0; r; theta; phi0; pr; p_theta];
        lambda_span = [0 50];

        [~, Y] = ode45(@(lambda,Y) geodesic_eqns_3D(lambda, Y, M, a, p_t, p_phi), lambda_span, Y0, opts);

        r_traj = Y(:,2);
        theta_traj = Y(:,3);
        phi_traj = Y(:,4);

        [x, y, z] = sph2cart(phi_traj, pi/2 - theta_traj, r_traj);

        plot3(x, y, z, 'y', 'LineWidth', 1.2);
    end
end

% Event horizon (approximate sphere)
[Xs,Ys,Zs] = sphere(50);
Rh = M + sqrt(M^2 - a^2);
surf(Rh*Xs, Rh*Ys, Rh*Zs, 'FaceAlpha', 0.1, 'EdgeColor', 'none');

% Ergosphere
theta = linspace(0, pi, 100);
phi = linspace(0, 2*pi, 200);
[thetaGrid, phiGrid] = meshgrid(theta, phi);

% Ergosphere radius as a function of theta
rErgo = M + sqrt(M^2 - a^2 * cos(thetaGrid).^2);

% Convert to Cartesian
x = rErgo .* sin(thetaGrid) .* cos(phiGrid);
y = rErgo .* sin(thetaGrid) .* sin(phiGrid);
z = rErgo .* cos(thetaGrid);

% Plot
surf(x, y, z, 'FaceAlpha', 0.05, 'EdgeColor', 'none');

hold off;

% Grid in spherical coords
r_vals = linspace(2*M, 8*M, 15);  % avoid inside horizon
theta_vals = linspace(0, pi, 15);
phi_vals = linspace(0, 2*pi, 20);

[R, Theta, Phi] = meshgrid(r_vals, theta_vals, phi_vals);

% Metric quantities
Sigma = R.^2 + a^2 .* cos(Theta).^2;
Delta = R.^2 - 2*M.*R + a^2;

g_tphi = - (2*M*a.*R .* sin(Theta).^2) ./ Sigma;
g_phiphi = ((R.^2 + a^2).^2 - a^2 .* Delta .* sin(Theta).^2) .* sin(Theta).^2 ./ Sigma;

% ZAMO angular velocity
omega = - g_tphi ./ g_phiphi;

% Convert to Cartesian for plotting
X = R .* sin(Theta) .* cos(Phi);
Y = R .* sin(Theta) .* sin(Phi);
Z = R .* cos(Theta);

% Tangential unit vector in phi direction
e_phi_x = -sin(Phi);
e_phi_y = cos(Phi);
e_phi_z = zeros(size(Phi));

% Frame-drag velocity field
U = omega .* e_phi_x;
V = omega .* e_phi_y;
W = omega .* e_phi_z;

% Plot
figure;
hold
quiver3(X, Y, Z, U, V, W, 'r');
axis equal;
xlabel('x'); ylabel('y'); zlabel('z');
title('ZAMO Frame-Drag Angular Velocity Field');
grid on;
view(45, 30);

% Event horizon and ergosphere
[shx, shy, shz] = sphere(50);
Rh = M + sqrt(M^2 - a^2);
surf(Rh*shx, Rh*shy, Rh*shz, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'FaceColor', 'k');
surf(x, y, z, 'FaceAlpha', 0.05, 'EdgeColor', 'none');
view(3);

hold off

% Non-uniform coordinate grids (clustered near 0)
Nx = 40; Ny = 40; Nz = 1;
xvals = sinh(linspace(-2, 2, Nx));   % cluster near 0 in x
yvals = sinh(linspace(-2, 2, Ny));   % cluster near 0 in y
zvals = 0;

[xx, yy, zz] = meshgrid(xvals, yvals, zvals);

% Boyer–Lindquist coordinates
r = sqrt(xx.^2 + yy.^2 + zz.^2);
theta = acos(zz ./ r);
phi = atan2(yy, xx);

% Kerr quantities
Sigma = r.^2 + a^2 .* cos(theta).^2;
Delta = r.^2 - 2*M.*r + a^2;
A = (r.^2 + a^2).^2 - a^2 .* Delta .* sin(theta).^2;

% Physical ZAMO angular velocity
omega = 2*M*a.*r ./ A;

% Plot grid lines
figure;

for dt = 1:0.5:dtMax
    clf;
    axis equal; grid on;
    xlabel('x'); ylabel('y'); zlabel('z');
    title('Physically Accurate Frame Dragging Grid (Non-uniform)');
    
    hold on;
    
    % Apply frame dragging for Δt = 1
    phi_def = phi + dt*omega;

    % Convert back to Cartesian
    xd = r .* sin(theta) .* cos(phi_def);
    yd = r .* sin(theta) .* sin(phi_def);
    zd = r .* cos(theta);

    % Draw lines
    for k = 1:size(xx,1)
        for m = 1:size(xx,3)
            plot3(squeeze(xd(k,:,m)), squeeze(yd(k,:,m)), squeeze(zd(k,:,m)), 'b');
        end
    end
    for k = 1:size(xx,2)
        for m = 1:size(xx,3)
            plot3(squeeze(xd(:,k,m)), squeeze(yd(:,k,m)), squeeze(zd(:,k,m)), 'b');
        end
    end
    for k = 1:size(xx,1)
        for m = 1:size(xx,2)
            plot3(squeeze(xd(k,m,:)), squeeze(yd(k,m,:)), squeeze(zd(k,m,:)), 'b');
        end
    end

    % Event horizon and ergosphere
    [shx, shy, shz] = sphere(50);
    Rh = M + sqrt(M^2 - a^2);
    surf(Rh*shx, Rh*shy, Rh*shz, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'FaceColor', 'k');
    surf(x, y, z, 'FaceAlpha', 0.05, 'EdgeColor', 'none');
    hold off;
    view(45,30);
    drawnow;
    
    
    % Smoother visualization
    MM=getframe;
    
    
end

%%%-----------------------
function dYdLambda = geodesic_eqns_3D(~, Y, M, a, p_t, p_phi)
    r = Y(2);
    theta = Y(3);
    phi = Y(4);
    p_r = Y(5);
    p_theta = Y(6);

    ginv = inv_metric_Kerr_full(r, theta, M, a);
    gtt = ginv(1,1); grr = ginv(2,2); gthth = ginv(3,3);
    gtphi = ginv(1,4); gphiphi = ginv(4,4);

    dt = -gtt*p_t - gtphi*p_phi;
    dr = grr * p_r;
    dtheta = gthth * p_theta;
    dphi = -gtphi*p_t - gphiphi*p_phi;

    % Finite diff for derivatives wrt r and theta
    delta = 1e-5;
    ginv_r_plus = inv_metric_Kerr_full(r + delta, theta, M, a);
    ginv_r_minus = inv_metric_Kerr_full(r - delta, theta, M, a);
    dginv_dr = (ginv_r_plus - ginv_r_minus) / (2*delta);

    ginv_th_plus = inv_metric_Kerr_full(r, theta + delta, M, a);
    ginv_th_minus = inv_metric_Kerr_full(r, theta - delta, M, a);
    dginv_dth = (ginv_th_plus - ginv_th_minus) / (2*delta);

    ps = [p_t; p_r; p_theta; p_phi];
    dp_r = -0.5 * sum(sum(dginv_dr .* (ps * ps')));
    dp_theta = -0.5 * sum(sum(dginv_dth .* (ps * ps')));

    dYdLambda = [dt; dr; dtheta; dphi; dp_r; dp_theta];
end

%%%-----------------------
function ginv = inv_metric_Kerr_full(r, theta, M, a)
    Sigma = r^2 + a^2 * cos(theta)^2;
    Delta = r^2 - 2*M*r + a^2;
    A = (r^2 + a^2)^2 - a^2 * Delta * sin(theta)^2;

    ginv = zeros(4,4);

    ginv(1,1) = -(A / (Sigma * Delta));
    ginv(1,4) = -2*M*a*r / (Sigma * Delta);
    ginv(4,1) = ginv(1,4);
    ginv(2,2) = Delta / Sigma;
    ginv(3,3) = 1 / Sigma;
    ginv(4,4) = (Sigma - 2*M*r) / (Sigma * Delta * sin(theta)^2);
end

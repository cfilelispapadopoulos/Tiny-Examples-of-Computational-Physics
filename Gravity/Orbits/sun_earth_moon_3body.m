clear;
% clc;

% Constants
m = [1.989e30;5.972e24;7.34767309e22];
Tmax = 365.242199*24*3600;
G = 6.67428e-11;

% Stricter tolerance
options = odeset('RelTol',1e-8,'AbsTol',1e-10);

% Initial conditions
y0 = zeros(12,1);

% x - Positions
y0(1) = 0; y0(2) = 147098074e3; y0(3) = y0(2) + 384400e3;
% y - Positions
y0(4:6) = 0;
% x - Velocities
y0(7:9) = 0;
% y - Velocities
y0(10) = 0; y0(11) = 30.3e3; y0(12) = y0(11) + 1022;

% Solve with ode45
sol = ode45(@(t,p) F(t,p,G,m),[0 Tmax],y0,options);

% Plot trajectories
plot(sol.y(1,:),sol.y(4,:),'o');hold
plot(sol.y(2,:),sol.y(5,:));
plot(sol.y(3,:),sol.y(6,:));
xlabel('x (km)');ylabel('y (km)');
legend('Sun','Earth','Moon');hold off
figure
plot(sol.x / (24*3600),sqrt((sol.y(1,:)-sol.y(2,:)).^2+(sol.y(4,:)-sol.y(5,:)).^2) / 1000);hold
plot(sol.x / (24*3600),sqrt((sol.y(3,:)).^2+(sol.y(6,:)).^2) / 1000);
xlabel('Time (Days)');ylabel('Distance (km)')
figure
plot(sol.x / (24*3600),sqrt((sol.y(2,:)-sol.y(3,:)).^2+(sol.y(5,:)-sol.y(6,:)).^2) / 1000);hold
xlabel('Time (Days)');ylabel('Distance (km)')

% System of ODEs
function y = F(t,y0,G,m)
    % Initialize output
    y = zeros(12,1);
    % x-Positions
    y(1:3) = y0(7:9);
    % y-Positions
    y(4:6) = y0(10:12);
    % Distances
    dx = y0(1:3)-y0(1:3)';
    dy = y0(4:6)-y0(4:6)';
    d = sqrt(dx.^2+dy.^2);
    % x-Velocities
    y(7) = -G * m(2) / d(1,2)^3 * dx(1,2)...
        - G *m(3) / d(1,3)^3 * dx(1,3);
    y(8) =  G * m(1) / d(1,2)^3 * dx(1,2)...
        - G *m(3) / d(2,3)^3 * dx(2,3);
    y(9) =  G * m(1) / d(1,3)^3 * dx(1,3)...
        + G *m(2) / d(2,3)^3 * dx(2,3);
    % y-Velocities
    y(10) = -G * m(2) / d(1,2)^3 * dy(1,2)...
        - G *m(3) / d(1,3)^3 * dy(1,3);
    y(11) =  G * m(1) / d(1,2)^3 * dy(1,2)...
        - G *m(3) / d(2,3)^3 * dy(2,3);
    y(12) =  G * m(1) / d(1,3)^3 * dy(1,3)...
        + G *m(2) / d(2,3)^3 * dy(2,3);
end
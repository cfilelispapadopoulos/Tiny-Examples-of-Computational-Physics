clear;
clc;

% Constants
G = 6.67428e-11;
M = 1.989e30; % Mass of the sun
m = 5e2; % Mass of the satelite
Tmax = 1.5*365.242199*24*3600; %Year in seconds
r0 = 147098074e3; % Initial distance of the lower orbit
v0 = 30.3e3; % Initial velocity of the lower orbit
% Number of trajectories to check in order
% to find the correct one
Nst = 5000;

% Percentage of the trajectory where the
% starting and ending poitns lie
per1 = 0.2;
per2 = 0.7;

% Get and plot the 2 trajectories
hold
plot(0,0,'ok');
[t1,r1,th1] = getTrajectory(G,M,m,Tmax,r0,0,v0);
plot(r1.*cos(th1),r1.*sin(th1),'k');

[t2,r2,th2] = getTrajectory(G,M,m,Tmax,1.2*r0,0,0.95*v0);
plot(r2.*cos(th2),r2.*sin(th2),'k');

% Convert percentage to index
ind1 = min(max(1,round(per1*length(t1))),round(per1*length(t1)));
ind2 = min(max(1,round(per2*length(t2))),round(per2*length(t2)));

% Get velocities
v1 = (th1(ind1+1)-th1(ind1-1))/(t1(ind1+1)-t1(ind1-1))*r1(ind1);
v2 = (th2(ind2+1)-th2(ind2-1))/(t2(ind2+1)-t2(ind2-1))*r2(ind2);

% Plot beginning and ending of the transfer
plot(r1(ind1)*cos(th1(ind1)),r1(ind1)*sin(th1(ind1)),'or');
plot(r2(ind2)*cos(th2(ind2)),r2(ind2)*sin(th2(ind2)),'or');

% Shooting method to find the Hohmann transfer orbit
[t,r,th,dvs] = getTransfer(r1(ind1),th1(ind1),v1,r2(ind2),th2(ind2),v2,G,m,M,Tmax,Nst);
plot(r.*cos(th),r.*sin(th),'r--');
xlabel('x (m)');ylabel('y (m)');
hold off

% Get transfer orbit by the shooting method
function [tt,rr,thh,dvs] = getTransfer(r1,th1,v1,r2,th2,v2,G,m,M,Tmax,Nst)
    % Compute increments of velocity
    dv = linspace(0,0.15*v1,Nst);
    mind = inf;
    for i=1:Nst
        % Compute Angular Momentum
        L = m * (v1+dv(i)) * r1;

        % Initial Condition
        y0 = [r1;0;th1];

        % Stricter tolerance
        options=odeset('RelTol',1e-10,'AbsTol',1e-10);

        % Solve with ode45
        sol = ode45(@(t,p) F(t,p,G,m,M,L),[0 Tmax],y0,options);

        % Extract vectors from structure
        t = sol.x;
        r = sol.y(1,:);
        th = sol.y(3,:);
        
        % Compute distances from ending point
        d = sqrt((r2*cos(th2)-r'.*cos(th')).^2+(r2*sin(th2)-r'.*sin(th')).^2);
        % Get minimum
        [~,indm] = min(d);
        % Check in new orbit is closer than the older one and
        % safeguard against the ending that causes error
        if indm >= length(d)
            continue;
        end
        if (d(indm) < mind)
            mind = d(indm);
            % Get DV impulses
            dvs(1) = dv(i);
            dvs(2) = v2-(th(indm+1)-th(indm-1))/(t(indm+1)-t(indm-1))*r(1,indm);
            % Retain the trajectory to be returned
            tt = t(1:indm);
            rr = r(1:indm);
            thh = th(1:indm);
        end
    end
    
    
end


% Get trajectory
function [t,r,th] = getTrajectory(G,M,m,Tmax,r0,th0,v0)
    % Compute Angular Momentum
    L = m * v0 * r0;

    % Initial Condition
    y0 = [r0;0;th0];

    % Stricter tolerance
    options=odeset('RelTol',1e-10,'AbsTol',1e-10);

    % Solve with ode45
    sol = ode45(@(t,p) F(t,p,G,m,M,L),[0 Tmax],y0,options);

    % Plot results
    t = sol.x;
    r = sol.y(1,:);
    th = sol.y(3,:);
end

% Function
function y = F(t,y0,G,m,M,L)
    % Initialize output
    y = zeros(3,1);
    % Equation r\dot = v_r
    y(1) = y0(2);
    % Equation v_r\dot = L^2 / (m^2 r^3) - G M / r^2
    y(2) =  L^2 ./ (m^2 * y0(1)^3) - G * M ./ y0(1).^2;
    % Equation \phi\dot = L / (m r ^2)
    y(3) = L ./ (m * y0(1).^2);
end

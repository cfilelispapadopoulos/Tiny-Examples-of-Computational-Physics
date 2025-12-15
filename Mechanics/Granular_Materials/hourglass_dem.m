% 2D DEM hourglass,
clear; close all; rng(1);

%% ---------------- Simulation Parameters ----------------
N = 600;                    
Rmean = 0.003;
Rstd  = 0.001;
boxW  = 0.3;
boxH  = 0.6;
neckW = 0.04;
neckY = 0.05;

g = 9.81;
rho = 2500;
kn = 1e5;     
kt = 2e4;     
gamma_n = 50;  
mu = 0.4;     
mu_roll = 0.02;

CFL = 0.1;      
t_end = 1.7;
plotEvery = 20;

relaxSteps = 1000;
gamma_relax = 200;

%% ---------------- Hourglass Geometry ----------------
neckHalf = neckW/2;
topY =  boxH/2;
botY = -boxH/2;

% Vertical neck region
neckTop = neckY;          % top of neck
neckBottom = neckY - 0.1; % bottom of vertical neck segment

% Walls: upper funnel
walls = [
    -boxW/2, topY,  boxW/2, topY;       % top horizontal
    -boxW/2, topY, -neckHalf, neckTop;  % left slope
     boxW/2, topY,  neckHalf, neckTop;  % right slope
];

% Vertical neck walls (middle portion)
walls = [walls;
    -neckHalf, neckTop, -neckHalf, neckBottom;  % left vertical neck
     neckHalf, neckTop,  neckHalf, neckBottom]; % right vertical neck

% Lower funnel walls (angled)
walls = [walls;
    -neckHalf, neckBottom, -boxW/2, botY;  % left lower funnel
     neckHalf, neckBottom,  boxW/2, botY]; % right lower funnel

% Bottom horizontal wall
walls = [walls;
    -boxW/2, botY, boxW/2, botY];

polyX = [-boxW/2, boxW/2, neckHalf, -neckHalf];
polyY = [topY, topY, neckTop, neckTop];

%% ---------------- Particle Initialization (log-normal) ----------------
Rs = lognrnd(log(Rmean)-0.5*log(1+(Rstd/Rmean)^2), sqrt(log(1+(Rstd/Rmean)^2)), N,1);
rmax = max(Rs);

pos = zeros(N,2);
placed = 0;
while placed < N
    x = (rand-0.5)*boxW;
    y = neckTop + rand*(topY-neckTop);
    if inpolygon(x,y,polyX,polyY)
        placed = placed + 1;
        pos(placed,:) = [x y];
    end
end
vel = zeros(N,2);
omega = zeros(N,1);
m = rho*pi.*Rs.^2;
I = 0.5*m.*Rs.^2;

%% ---------------- Visualization ----------------
fig = figure('Color','w','Position',[100 100 600 900]);
ax = axes('Parent',fig); hold on; axis equal; box on;
xlim([-boxW/2 boxW/2]); ylim([botY topY]);
title('DEM Hourglass');

for i=1:size(walls,1)
    plot([walls(i,1) walls(i,3)],[walls(i,2) walls(i,4)],'k-','LineWidth',2);
end

cmap = jet(256);
Rs_norm = (Rs - min(Rs)) / (max(Rs) - min(Rs));
colorIdx = round(Rs_norm*255) + 1;
colors = cmap(colorIdx,:);
sc = scatter(pos(:,1), pos(:,2), (Rs*2000).^2, colors, 'filled');
drawnow;

%% ---------------- Compute stable dt ----------------
dt = CFL * min(sqrt(m./kn));

%% ---------------- Initial Relaxation ----------------
force = computeForcesHertzFull(pos, vel, omega, Rs, m, I, walls, kn, kt, gamma_relax, mu, mu_roll, g);
acc = force(:,1:2)./m;
alpha = force(:,3)./I;

for step=1:relaxSteps
    pos = pos + vel*dt + 0.5*acc*dt^2;
    omega = omega + 0.5*alpha*dt;
    force_new = computeForcesHertzFull(pos, vel, omega, Rs, m, I, walls, kn, kt, gamma_relax, mu, mu_roll, g);
    acc_new = force_new(:,1:2)./m;
    alpha_new = force_new(:,3)./I;
    
    vel = vel + 0.5*(acc + acc_new)*dt;
    omega = omega + 0.5*(alpha + alpha_new)*dt;
    
    acc = acc_new; alpha = alpha_new; force = force_new;
end

%% ---------------- DEM Simulation ----------------
t = 0; step = 0;
while t < t_end
    step = step + 1;
    
    pos = pos + vel*dt + 0.5*acc*dt^2;
    omega = omega + 0.5*alpha*dt;
    
    force_new = computeForcesHertzFull(pos, vel, omega, Rs, m, I, walls, kn, kt, gamma_n, mu, mu_roll, g);
    acc_new = force_new(:,1:2)./m;
    alpha_new = force_new(:,3)./I;
    
    dt = CFL * min(sqrt(m./kn));
    
    vel = vel + 0.5*(acc + acc_new)*dt;
    omega = omega + 0.5*(alpha + alpha_new)*dt;
    
    acc = acc_new; alpha = alpha_new; force = force_new;
    
    t = t + dt;
    
    if mod(step,plotEvery)==0
        set(sc,'XData',pos(:,1),'YData',pos(:,2));
        drawnow limitrate;
    end
end

disp('Simulation complete.');

%% ---------------- Function: Hertzian + Rolling Friction (Full Geometry) ----------------
function F = computeForcesHertzFull(pos, vel, omega, R, m, I, walls, kn, kt, gamma_n, mu, mu_roll, g)
    N = size(pos,1);
    F = zeros(N,3);
    
    F(:,2) = -m*g;
    
    maxR = max(R);
    Mdl = KDTreeSearcher(pos);
    neighbors = rangesearch(Mdl, pos, 2*maxR);
    
    %% Particle–particle
    for i=1:N
        ni = neighbors{i};
        for j = ni
            if j<=i, continue; end
            rij = pos(j,:) - pos(i,:);
            dist = norm(rij);
            overlap = R(i)+R(j)-dist;
            if overlap>0
                n = rij/dist;
                dv = vel(j,:) - vel(i,:);
                vn = dot(dv,n);
                fn_mag = kn*overlap^(1.5) - gamma_n*vn;
                fn_mag(fn_mag<0)=0;
                fn = fn_mag*n;
                
                tdir = [-n(2) n(1)];
                vt = dot(dv,tdir) + R(i)*omega(i) + R(j)*omega(j);
                ft_mag = -kt*vt;
                if abs(ft_mag) > mu*fn_mag, ft_mag = -mu*fn_mag*sign(vt); end
                ft = ft_mag*tdir;
                
                torque_i = -mu_roll*fn_mag*R(i)*sign(omega(i));
                torque_j = -mu_roll*fn_mag*R(j)*sign(omega(j));
                
                F(i,:) = F(i,:) + [-(fn+ft), torque_i];
                F(j,:) = F(j,:) + [(fn+ft), torque_j];
            end
        end
    end
    
    %% Particle–wall
    for wi=1:size(walls,1)
        x1 = walls(wi,1:2); x2 = walls(wi,3:4);
        seg = x2 - x1; segL2 = dot(seg,seg);
        for i=1:N
            v = pos(i,:) - x1;
            tproj = max(0,min(1,dot(v,seg)/segL2));
            closest = x1 + tproj*seg;
            
            rij = pos(i,:) - closest;
            dist = norm(rij);
            overlap = R(i) - dist;
            if overlap>0
                if dist==0, n=[0 1]; else n=rij/dist; end
                vn = dot(vel(i,:),n);
                fn_mag = kn*overlap - gamma_n*vn;
                if fn_mag<0, fn_mag=0; end
                fn = fn_mag*n;
                
                tdir = [-n(2) n(1)];
                vt = dot(vel(i,:),tdir);
                ft_mag = -kt*vt;
                if abs(ft_mag)>mu*fn_mag, ft_mag = -mu*fn_mag*sign(vt); end
                ft = ft_mag*tdir;
                
                F(i,1:2) = F(i,1:2) + fn + ft;
                F(i,3) = F(i,3) - mu_roll*fn_mag*R(i)*sign(omega(i));
                
                % Project particle outside wall slightly
                pos(i,:) = closest + n*(R(i)+1e-6);
            end
        end
    end
end

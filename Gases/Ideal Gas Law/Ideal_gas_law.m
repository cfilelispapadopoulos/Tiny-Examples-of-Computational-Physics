% Simulation of oxygen filled piston with two compartments in 3D 
% using Ideal Gas Law and under prescribed temperatures


clear;
clc;
% All units are in SI
% Width of piston
w = 0.3;
% Height of piston
h = 0.1;
% Length of piston
l = 0.5;
% Mass of piston
Mb = 0.1;
% Initial position of the piston
p = 0.4;
pin = p;
% Safeguard against choosing an initial position that towards the very end
% of the available length of the piston
p = min(p, 0.75 * l);

% Spring coeff
spr_coeff = 3;

% Maximum simulation time
T_max = 10;
% Time step
t_step = 0.005;

% Boltzmann's constant (m^2 kg s^-2 K^-1)
kB = 1.380649e-23;
% Ideal Gas Constant (J/(K mol))
R = 8.314;
% Avogadro's number
NA = 6.0221408e+23;
% Number of moles
n = [2 2];
% Number of particles (for visualization)
N = [100 100];
% Mass of each particle (in simulation)
% 1 mol of O_2 is 32g (molar mass)
mmol = 32e-3;
M = n * mmol ./ N;
% Pressure in Pa
T = [300 300];
% Temperature of Gas in K from ideal gas law
P = (n .* R .* T) ./ [p * w * h, (l - p) * w * h];

% Create initial positions and velocities
[pp,vp] = create_initial_pos_vel(N, l, w, h, p, kB, T, mmol, NA);

% Accelerations of particles
ap = zeros(size(pp));

% Velocity of wall
vw = 0;
% Acceleration of wall
aw = 0;

% Create a figure
figure
% Loop over
% Initialize time
time = 0;
while time < T_max
    % Clear previous plot
    clf
    % Compute adaptive time step
    t_step=min(0.001,0.01/max(max((abs(vp)))));
    % Hold the plots in order to simultaneously present the geometry
    % and particles
    hold
    % Cleate geometry
    create_geometry(l, w, h, p);
    
    % Leapfrog KDK (particles)
    % Kick (1/2)
    vp = vp + ap * t_step / 2.0; 
    % Drift
    pp = pp + vp * t_step;
    % Kick (1/2)
    vp = vp + ap * t_step / 2.0;

    % Leapfrog KDK (wall)
    % Kick (1/2)
    vw = vw + aw * t_step / 2.0; 
    % Drift
    p = p + vw * t_step;
 
    % Kick (1/2)
    vw = vw + aw * t_step / 2.0;

    % Check position with respect to the box
    [pp, vp] = checkPos(pp, vp, l, w, h, p);
    
    % Recompute pressure
    P = (n .* R .* T) ./ [p * w * h, (l - p) * w * h];
    % Compute the acceleration due to pressure
    aw = barrierAcc(P, h, w, Mb);
    disp(aw)
    % Plot particles
    scatter3(pp(:,1), pp(:,2), pp(:,3), [],sqrt(sum(vp.^2,2))./max(sqrt(sum(vp.^2,2))),'filled');colormap(autumn);

    % Update time
    time = time + t_step;
    % Added for visualization reasons
    MM=getframe;
    % Remove hold to move to next frame
    hold off
end

% Function that creates the geometry of piston and spring
% Inputs:
%   w_  : Width of piston
%   h_  : Height of piston
%   t_  : Thickness of piston
%   p_  : Position of piston
function create_geometry(l_, w_, h_, p_)
    % Parallelepiped 
    coord = [...
        0 0 0;
        p_ 0 0;
        p_ w_ 0;
        0 w_ 0;
        0 0 h_;
        p_ 0 h_;
        p_ w_ h_;
        0 w_ h_;
        p_ 0 0;
        l_ 0 0;
        l_ w_ 0;
        p_ w_ 0;
        p_ 0 h_;
        l_ 0 h_;
        l_ w_ h_;
        p_ w_ h_;        
        ];
    
    idx = [4 8 5 1 4; 1 5 6 2 1; 2 6 7 3 2; 3 7 8 4 3; 5 8 7 6 5; 1 4 3 2 1;]';
    xc = coord(:,1); yc = coord(:,2); zc = coord(:,3);
    
    patch(xc(idx), yc(idx), zc(idx), 'r', 'facealpha', 0.05);
    patch(xc(idx+8), yc(idx+8), zc(idx+8), 'r', 'facealpha', 0.05);
    view(3);
end

% Function that creates the initial positions using uniform random numbers
% and velocities following Maxwell - Boltzmann distribution
% Inputs:
%   N_     : Number of particles
%   w_     : Width of the piston
%   h_     : Height of the piston
%   l_     : Length of piston (towards x axis)
%   p_     : Position of the piston with respect to reference point
%   kB_    : Boltzmann constant
%   T_     : Temperature (K)
%   mmol_  : Molar mass of O_2
%   NA_    : Avogadro's number
% Outputs:
%   pos_   : Position of all particles in both sides
%   vel_   : Velocities of particles in both sides
function [pos_, vel_] = create_initial_pos_vel(N_, l_, w_, h_, p_, kB_, T_, mmol_, NA_)
    % Create position of particles in the cylinder
    pos_ = rand(N_(1),3).*[p_, w_, h_];
    pos_ = [pos_; p_*[1 0 0]+rand(N_(1),3).*[l_ - p_, w_, h_]];
    v = sqrt((kB_ * T_) / (mmol_ / NA_));
    vel_ = v(1) * randn(N_(1), 3);
    vel_ = [vel_; v(2) * randn(N_(2), 3)];
end

% Function that checks position of particles 
% Inputs:
%   pos_:   Position of all particles
%   cel_:   Velocities of all particles
%   l_  :   Length of piston
%   w_  :   Weight of piston
%   h_  :   Height of piston
%   p_  :   Position of piston relative to refernce (0,0)
% Outputs:
%   pos_:   New positions after collision
%   vel_:   New velocities after collision
function [pos_, vel_] = checkPos(pos_, vel_, l_, w_, h_, p_)
    N = size(pos_,1)/2;
    
    % Check and reflect (elastic collisions)
    % Side 1 Half bodies
    ind = (pos_(1:N,1) <= 0);
    pos_(ind,1) = - pos_(ind,1);
    vel_(ind,1) = - vel_(ind,1);
    
    % Barrier
    % Left
    ind = (pos_(1:N,1) >= p_);
    pos_(ind,1) = 2 * p_ - pos_(ind,1);
    vel_(ind,1) = - vel_(ind,1);
    
    % Right
    ind = find(pos_(N+1:end,1) <= p_);
    if ~isempty(ind)
        pos_(N+ind,1) = 2 * p_ - pos_(N+ind,1);
        vel_(N+ind,1) = - vel_(N+ind,1);
    end
    
    % Side 2 Half bodies
    ind = find(pos_(N+1:end,1) >= l_);
    if ~isempty(ind)
        pos_(N+ind,1) = 2 * l_ - pos_(N+ind,1);
        vel_(N+ind,1) = - vel_(N+ind,1);
    end
    
    % Side 3
    ind = (pos_(:,2) >= w_);
    pos_(ind,2) = 2 * w_ - pos_(ind,2);
    vel_(ind,2) = - vel_(ind,2);
    
    % Side 4
    ind = (pos_(:,2) <= 0);
    pos_(ind,2) = - pos_(ind,2);
    vel_(ind,2) = - vel_(ind,2);

    % Side 5
    ind = (pos_(:,3) >= h_);
    pos_(ind,3) = 2 * h_ - pos_(ind,3);
    vel_(ind,3) = - vel_(ind,3);
    
    % Side 6
    ind = (pos_(:,3) <= 0);
    pos_(ind,3) = - pos_(ind,3);
    vel_(ind,3) = - vel_(ind,3);        
end

% Function to compute acceleration of the barrier due to opposing forces
% Inputs:
%   P         : Pressure in each compartment
%   h_        : Height of the box
%   w_        : Width of the box
%   Mb_       : Mass of the barier
% Outputs:
%   aw_       : Acceleration of the piston
function [aw_] = barrierAcc(P, h_, w_, Mb)
    aw_ = (P(1)-P(2)) * h_ * w_ / Mb;
end
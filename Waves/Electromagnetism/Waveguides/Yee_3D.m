%% 3D FDTD Waveguide with stable CPML on x boundaries and PEC walls
clc; clear;

%% --- Physical constants ---
c0 = 3e8; eps0 = 8.854e-12; mu0 = 4*pi*1e-7;

%% --- Grid parameters ---
Nx = 100; Ny = 80; Nz = 80;
dx = 1e-3; dy = 1e-3; dz = 1e-3;
dt = 0.6/(c0*sqrt(1/dx^2 + 1/dy^2 + 1/dz^2));  % slightly smaller CFL for stability
Nt = 1000;

%% --- CPML parameters (x only) ---
Npml = 20;
sigma_max = 10;    % moderate damping
kappa_max = 3;     % typical value
alpha_max = 0.05;  % small
m = 4;             % polynomial order

sigma_x = zeros(Nx,1); kappa_x = ones(Nx,1); alpha_x = zeros(Nx,1);
for i = 1:Npml
    x_norm = (Npml-i+0.5)/Npml;
    sigma_val = sigma_max * x_norm^m;         % smooth quartic grading
    kappa_val = 1 + (kappa_max-1)*x_norm^m;
    alpha_val = alpha_max*(1-x_norm);
    
    sigma_x(i) = sigma_val; sigma_x(end-i+1) = sigma_val;
    kappa_x(i) = kappa_val; kappa_x(end-i+1) = kappa_val;
    alpha_x(i) = alpha_val; alpha_x(end-i+1) = alpha_val;
end

% Precompute CPML coefficients
bx = exp(-(sigma_x./kappa_x + alpha_x)*dt/eps0);
bx = reshape(bx,[Nx,1,1]);
kappa_x = reshape(kappa_x,[Nx,1,1]);

%% --- Initialize fields ---
Ex = zeros(Nx,Ny,Nz); Ey = zeros(Nx,Ny,Nz); Ez = zeros(Nx,Ny,Nz);
Hx = zeros(Nx,Ny,Nz); Hy = zeros(Nx,Ny,Nz); Hz = zeros(Nx,Ny,Nz);

Exy = zeros(Nx,Ny,Nz); Exz = zeros(Nx,Ny,Nz);
Eyx = zeros(Nx,Ny,Nz); Eyz = zeros(Nx,Ny,Nz);
Ezx = zeros(Nx,Ny,Nz); Ezy = zeros(Nx,Ny,Nz);

Psi_Exy = zeros(Nx,Ny,Nz); Psi_Exz = zeros(Nx,Ny,Nz);
Psi_Eyx = zeros(Nx,Ny,Nz); Psi_Eyz = zeros(Nx,Ny,Nz);
Psi_Ezx = zeros(Nx,Ny,Nz); Psi_Ezy = zeros(Nx,Ny,Nz);

%% --- PEC obstacle ---
obs_size = 10;
obs_x = round(Nx/2-obs_size/2):round(Nx/2+obs_size/2);
obs_y = round(Ny/2-obs_size/2):round(Ny/2+obs_size/2);
obs_z = round(Nz/2-obs_size/2):round(Nz/2+obs_size/2);

%% --- Source parameters ---
t0 = 35; spread = 10;
source_pos = [30, round(Ny/2), round(Nz/2)];

%% --- Energy History ---
Uhist = zeros(Nt,1);

%% --- Main FDTD loop ---
figure;
for n = 1:Nt
    %% --- Update H fields ---
    Hx(1:end-1,1:end-1,1:end-1) = Hx(1:end-1,1:end-1,1:end-1) ...
        - dt/(mu0*dy)*(Ez(1:end-1,2:end,1:end-1)-Ez(1:end-1,1:end-1,1:end-1)) ...
        + dt/(mu0*dz)*(Ey(1:end-1,1:end-1,2:end)-Ey(1:end-1,1:end-1,1:end-1));
    
    Hy(1:end-1,1:end-1,1:end-1) = Hy(1:end-1,1:end-1,1:end-1) ...
        - dt/(mu0*dz)*(Ex(1:end-1,1:end-1,2:end)-Ex(1:end-1,1:end-1,1:end-1)) ...
        + dt/(mu0*dx)*(Ez(2:end,1:end-1,1:end-1)-Ez(1:end-1,1:end-1,1:end-1));
    
    Hz(1:end-1,1:end-1,1:end-1) = Hz(1:end-1,1:end-1,1:end-1) ...
        - dt/(mu0*dx)*(Ey(2:end,1:end-1,1:end-1)-Ey(1:end-1,1:end-1,1:end-1)) ...
        + dt/(mu0*dy)*(Ex(1:end-1,2:end,1:end-1)-Ex(1:end-1,1:end-1,1:end-1));
    
    %% --- Update E fields with stable CPML ---
    % Ex
    dHz_dy = (Hz(:,2:Ny,2:Nz-1) - Hz(:,1:Ny-1,2:Nz-1))/dy;
    dHy_dz = (Hy(:,2:Ny-1,2:Nz) - Hy(:,2:Ny-1,1:Nz-1))/dz;
    
    % update split accumulators (standard Yee integration part)
    Exy(:,2:Ny,2:Nz-1) = Exy(:,2:Ny,2:Nz-1) + (dt/eps0) * dHz_dy;
    Exz(:,2:Ny-1,2:Nz) = Exz(:,2:Ny-1,2:Nz) - (dt/eps0) * dHy_dz;
    
    Psi_Exy(:,2:Ny,2:Nz-1) = bx(:,1,1).*Psi_Exy(:,2:Ny,2:Nz-1) + dt/eps0*dHz_dy./kappa_x(:,1,1);
    Psi_Exz(:,2:Ny-1,2:Nz) = bx(:,1,1).*Psi_Exz(:,2:Ny-1,2:Nz) - dt/eps0*dHy_dz./kappa_x(:,1,1);
    
    Ex(:,2:Ny-1,2:Nz-1) = Exy(:,2:Ny-1,2:Nz-1) + Exz(:,2:Ny-1,2:Nz-1) + Psi_Exy(:,2:Ny-1,2:Nz-1) + Psi_Exz(:,2:Ny-1,2:Nz-1);
    
    % Ey
    dHz_dx = (Hz(2:Nx,2:Ny-1,2:Nz-1) - Hz(1:Nx-1,2:Ny-1,2:Nz-1))/dx;
    dHx_dz = (Hx(2:Nx-1,2:Ny-1,2:Nz) - Hx(2:Nx-1,2:Ny-1,1:Nz-1))/dz;
    
    % update split accumulators (standard Yee integration part)
    Eyx(2:Nx,2:Ny-1,2:Nz-1) = Eyx(2:Nx,2:Ny-1,2:Nz-1) - (dt/eps0) * dHz_dx;
    Eyz(2:Nx-1,2:Ny-1,2:Nz) = Eyz(2:Nx-1,2:Ny-1,2:Nz) + (dt/eps0) * dHx_dz;
    
    Psi_Eyx(2:Nx,2:Ny-1,2:Nz-1) = bx(2:Nx,1,1).*Psi_Eyx(2:Nx,2:Ny-1,2:Nz-1) - dt/eps0*dHz_dx./kappa_x(2:Nx,1,1);
    Psi_Eyz(2:Nx-1,2:Ny-1,2:Nz) = bx(2:Nx-1,1,1).*Psi_Eyz(2:Nx-1,2:Ny-1,2:Nz) + dt/eps0*dHx_dz./kappa_x(2:Nx-1,1,1);
    Ey(2:Nx-1,2:Ny-1,2:Nz-1) = Eyx(2:Nx-1,2:Ny-1,2:Nz-1) + Eyz(2:Nx-1,2:Ny-1,2:Nz-1) + Psi_Eyx(2:Nx-1,2:Ny-1,2:Nz-1) + Psi_Eyz(2:Nx-1,2:Ny-1,2:Nz-1);
    
    % Ez
    dHy_dx = (Hy(2:Nx,2:Ny-1,2:Nz-1) - Hy(1:Nx-1,2:Ny-1,2:Nz-1))/dx;
    dHx_dy = (Hx(2:Nx-1,2:Ny,2:Nz-1) - Hx(2:Nx-1,1:Ny-1,2:Nz-1))/dy;
    
    % update split accumulators (standard Yee integration part)
    Ezx(2:Nx,2:Ny-1,2:Nz-1) = Ezx(2:Nx,2:Ny-1,2:Nz-1) + (dt/eps0) * dHy_dx;
    Ezy(2:Nx-1,2:Ny,2:Nz-1) = Ezy(2:Nx-1,2:Ny,2:Nz-1) - (dt/eps0) * dHx_dy;    
    
    Psi_Ezx(2:Nx,2:Ny-1,2:Nz-1) = bx(2:Nx,1,1).*Psi_Ezx(2:Nx,2:Ny-1,2:Nz-1) + dt/eps0*dHy_dx./kappa_x(2:Nx,1,1);
    Psi_Ezy(2:Nx-1,2:Ny,2:Nz-1) = bx(2:Nx-1,1,1).*Psi_Ezy(2:Nx-1,2:Ny,2:Nz-1) - dt/eps0*dHx_dy./kappa_x(2:Nx-1,1,1);
    Ez(2:Nx-1,2:Ny-1,2:Nz-1) = Ezx(2:Nx-1,2:Ny-1,2:Nz-1) + Ezy(2:Nx-1,2:Ny-1,2:Nz-1) + Psi_Ezx(2:Nx-1,2:Ny-1,2:Nz-1) + Psi_Ezy(2:Nx-1,2:Ny-1,2:Nz-1);
    
    %% --- Apply Gaussian source ---
    Ez(source_pos(1),source_pos(2),source_pos(3)) = Ez(source_pos(1),source_pos(2),source_pos(3)) + exp(-((n-t0)/spread)^2);
    
    %% --- Apply PEC obstacle ---
    Ex(obs_x,obs_y,obs_z) = 0;
    Ey(obs_x,obs_y,obs_z) = 0;
    Ez(obs_x,obs_y,obs_z) = 0;
    
    %% --- Apply PEC waveguide walls ---
    Ex(:,1,:) = 0; Ex(:,Ny,:) = 0;
    Ey(:,1,:) = 0; Ey(:,Ny,:) = 0;
    Ez(:,1,:) = 0; Ez(:,Ny,:) = 0;
    Ex(:,:,1) = 0; Ex(:,:,Nz) = 0;
    Ey(:,:,1) = 0; Ey(:,:,Nz) = 0;
    Ez(:,:,1) = 0; Ez(:,:,Nz) = 0;
    
    %% --- Visualization ---
    if mod(n,5)==0
        mid_plane = round(Nz/2);
        imagesc(squeeze(Ez(:,:,mid_plane))');
        colorbar; caxis([-0.0001 0.0001]);
        xlabel('X'); ylabel('Y');
        title(['Ez (V/m) at middle plane, t = ',num2str(n)]);
        axis equal tight; drawnow;
    end
    
    % Compute instantaneous energy at time step n
    u = 0.5*( eps0*(Ex.^2 + Ey.^2 + Ez.^2) + mu0*(Hx.^2 + Hy.^2 + Hz.^2) );
    U_total = sum(u(:)) * dx * dy * dz;
    Uhist(n) = U_total;
    
    fprintf('Total EM energy at step %d: %e J\n', n, U_total);    
end

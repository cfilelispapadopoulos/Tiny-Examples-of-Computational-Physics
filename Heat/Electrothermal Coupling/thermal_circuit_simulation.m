clear;
clc;

% Size of substrate
Lx = 0.2;
Ly = 0.2;

% Number of grid points per dimension
nx = 100;
ny = 100;

% Maximum number of iterations 
Nmax = 10;

% Convergence tolerance
tol = 1e-3;

% Temperature coefficient (per K)
alpha = 3e-3; % (Steel)

% Heat transfer coefficient
h = 5;
d = 3e-3; % Thickness of substrate
h_vol = 2 * h / d;

% Conductivity of stainless steel at room temperature
sigma0 = 1.45e6;

% Reference and ambient temperature (K)
T0 = 293;
Tinf = T0;

% Source Voltage (V)
V0 = 1;

% Mesh sizes
hx = Lx / nx;
hy = Ly / ny;

% Materials
materials(1).name = 'substrate';
materials(1).sigma = 1e-6;
materials(1).kappa = 0.3; % (W/(m K))

materials(2).name = 'wire';
materials(2).sigma = 5.8e7;
materials(2).kappa = 400; % (W/(m K))
materials(2).width = 1;

materials(3).name = 'resistor';
materials(3).sigma = @(T) sigma0 ./ (1 + alpha * (T - T0));
materials(3).kappa = 20; % (W/(m K))
materials(3).width = 2;
materials(3).length = 3;

% Matrices for cell sigmas and kappas
sigmas = zeros(nx+1,ny+1);
kappas = zeros(nx+1,ny+1);

% Substrate
sigmas(:,:) = materials(1).sigma;
kappas(:,:) = materials(1).kappa;

% Wiring
sigmas(nx/2-materials(2).width:nx/2+materials(2).width,1:ny/2) = materials(2).sigma;
kappas(nx/2-materials(2).width:nx/2+materials(2).width,1:ny/2) = materials(2).kappa;
sigmas(nx/2-nx/4-1:nx/2+nx/4+1,ny/2-materials(2).width:ny/2+materials(2).width) = materials(2).sigma;
kappas(nx/2-nx/4-1:nx/2+nx/4+1,ny/2-materials(2).width:ny/2+materials(2).width) = materials(2).kappa;
sigmas(nx/2-nx/4-1-materials(2).width+1:nx/2-nx/4-1+materials(2).width+1,ny/2+1-materials(2).width:ny-ny/10+materials(2).width) = materials(2).sigma;
kappas(nx/2-nx/4-1-materials(2).width+1:nx/2-nx/4-1+materials(2).width+1,ny/2+1-materials(2).width:ny-ny/10+materials(2).width) = materials(2).kappa;
sigmas(nx/2+nx/4+1-materials(2).width-1:nx/2+nx/4+1+materials(2).width-1,ny/2+1-materials(2).width:ny-ny/10+materials(2).width) = materials(2).sigma;
kappas(nx/2+nx/4+1-materials(2).width-1:nx/2+nx/4+1+materials(2).width-1,ny/2+1-materials(2).width:ny-ny/10+materials(2).width) = materials(2).kappa;
sigmas(nx/2-nx/4-1-materials(2).width+1:nx/2+nx/4+1+materials(2).width-1,ny-ny/10-materials(2).width:ny-ny/10+materials(2).width) = materials(2).sigma;
kappas(nx/2-nx/4-1-materials(2).width+1:nx/2+nx/4+1+materials(2).width-1,ny-ny/10-materials(2).width:ny-ny/10+materials(2).width) = materials(2).kappa;
sigmas(nx/2-materials(2).width:nx/2+materials(2).width,ny-ny/10+materials(2).width+1:end) = materials(2).sigma;
kappas(nx/2-materials(2).width:nx/2+materials(2).width,ny-ny/10+materials(2).width+1:end) = materials(2).kappa;

% Set current temp to ambient
T = T0*ones(nx+1,ny+1);

% Get max temp
maxT = max(T(:));

% Iterate towards equilibrium
for iter = 1:Nmax

    % Resistors (Kappas can be set outside of the loop for performance)
    sigmas(nx/2-materials(3).width:nx/2+materials(3).width,ny/4-materials(3).length:ny/4+materials(3).length) = materials(3).sigma(T(nx/2-materials(3).width:nx/2+materials(3).width,ny/4-materials(3).length:ny/4+materials(3).length));
    kappas(nx/2-materials(3).width:nx/2+materials(3).width,ny/4-materials(3).length:ny/4+materials(3).length) = materials(3).kappa;
    sigmas(nx/2-nx/4-1-materials(3).width+1:nx/2-nx/4-1+materials(3).width+1,3*ny/4-materials(3).length:3*ny/4+materials(3).length) = materials(3).sigma(T(nx/2-nx/4-1-materials(3).width+1:nx/2-nx/4-1+materials(3).width+1,3*ny/4-materials(3).length:3*ny/4+materials(3).length));
    kappas(nx/2-nx/4-1-materials(3).width+1:nx/2-nx/4-1+materials(3).width+1,3*ny/4-materials(3).length:3*ny/4+materials(3).length) = materials(3).kappa;
    sigmas(nx/2+nx/4+1-materials(3).width-1:nx/2+nx/4+1+materials(3).width-1,3*ny/4-materials(3).length:3*ny/4+materials(3).length) = materials(3).sigma(T(nx/2+nx/4+1-materials(3).width-1:nx/2+nx/4+1+materials(3).width-1,3*ny/4-materials(3).length:3*ny/4+materials(3).length));
    kappas(nx/2+nx/4+1-materials(3).width-1:nx/2+nx/4+1+materials(3).width-1,3*ny/4-materials(3).length:3*ny/4+materials(3).length) = materials(3).kappa;
    
    % Form coefficient matrix and Identity matrix
    A = form_coeff_matrix(sigmas,hx,hy);
    I = speye((nx+1)*(ny+1));

    % Apply Dirichlet Lines for the known Voltages
    A(nx/2-materials(2).width:nx/2+materials(2).width,:) = I(nx/2-materials(2).width:nx/2+materials(2).width,:);
    A((nx+1)*ny+(nx/2-materials(2).width:nx/2+materials(2).width),:) = I((nx+1)*ny+(nx/2-materials(2).width:nx/2+materials(2).width),:);

    % Right hand size
    Vrhs = zeros(nx+1,ny+1);
    Vrhs(nx/2-materials(2).width:nx/2+materials(2).width) = V0;

    % Solve linear system for potential
    V = reshape(A\Vrhs(:),nx+1,ny+1);
    
    % Form coefficient matrix for temperature
    B = form_coeff_matrix(kappas,hx,hy) + spdiags(h_vol * ones((nx+1) * (ny+1),1), 0, (nx+1) * (ny+1), (nx+1) * (ny+1));
    
    % Compute the gradient required for rhs of the thermal equation
    [Vx,Vy] = gradient(V,hx,hy);
    
    % Heat source
    Q = sigmas .* (Vx.^2 + Vy.^2) * hx * hy;     % Joule heating in W/m^3
    Q = Q(:);                                    % Flatten for matrix solve

    % Add Robin boundary conditions to thermal matrix
    % Copy thermal matrix B and RHS
    Trhs = Q + h_vol * Tinf;

    % Apply Robin BC on all four sides:
    for i = 1:(nx+1)
        % Top (j = ny+1)
        row = sub2ind([nx+1, ny+1], i, ny+1);
        B(row, row) = B(row, row) + h * hx / kappas(i, ny+1);
        Trhs(row) = Trhs(row) + h * hx * Tinf / kappas(i, ny+1);

        % Bottom (j = 1)
        row = sub2ind([nx+1, ny+1], i, 1);
        B(row, row) = B(row, row) + h * hx / kappas(i,1);
        Trhs(row) = Trhs(row) + h * hx * Tinf / kappas(i,1);
    end

    for j = 1:(ny+1)
        % Left (i = 1)
        row = sub2ind([nx+1, ny+1], 1, j);
        B(row, row) = B(row, row) + h * hy / kappas(1,j);
        Trhs(row) = Trhs(row) + h * hy * Tinf / kappas(1,j);

        % Right (i = nx+1)
        row = sub2ind([nx+1, ny+1], nx+1, j);
        B(row, row) = B(row, row) + h * hy / kappas(nx+1,j);
        Trhs(row) = Trhs(row) + h * hy * Tinf / kappas(nx+1,j);
    end

    % Solve temperature
    T = reshape(B \ Trhs, nx+1, ny+1);
    
    % Check change in max temp
    curr_maxT = max(T(:));
    if abs(curr_maxT-maxT) < tol
        disp(['Converged in ' num2str(iter) ' iterations'])
        break;
    end
    maxT = curr_maxT;
end

% Compute coordinates of the grid points
[xc,yc] = ndgrid(linspace(0,Lx,nx+1),linspace(0,Ly,ny+1));

% Plot voltage
figure;
contourf(xc, yc, V, 30, 'LineColor', 'none');xlabel('x');ylabel('y');title('V (V)');
colorbar;

% Compute voltage gradients
[Vy, Vx] = gradient(V, hy, hx);

% Compute current density components
Jx = -sigmas .* Vx;
Jy = -sigmas .* Vy;
J = sqrt(Jx.^2 + Jy.^2);

% Plot current density
figure;
contourf(xc, yc, J, 30, 'LineColor', 'none');xlabel('x');ylabel('y');title('|J| (A/m^2)');
colorbar;

% Optional: quiver plot of current vectors
figure;
quiver(xc, yc, Jx, Jy);
title('Current density vector field');
xlabel('x'); ylabel('y');
axis equal tight;

% Plot temp
figure;
contourf(xc, yc, T, 30, 'LineColor', 'none');xlabel('x');ylabel('y');title('T (K)');
colorbar;

function A = form_coeff_matrix(sigma_mat, dx, dy)
    [Nx, Ny] = size(sigma_mat);
    N = Nx * Ny;

    % Harmonic averaging of conductivities at cell interfaces

    % x-direction (horizontal interfaces between rows)
    sigma_xp = 2 ./ (1./sigma_mat(1:end-1,:) + 1./sigma_mat(2:end,:));  % interfaces in x
    sx_p = zeros(Nx, Ny); sx_p(1:end-1,:) = sigma_xp;
    sx_m = zeros(Nx, Ny); sx_m(2:end,:)   = sigma_xp;

    % y-direction (vertical interfaces between columns)
    sigma_yp = 2 ./ (1./sigma_mat(:,1:end-1) + 1./sigma_mat(:,2:end));  % interfaces in y
    sy_p = zeros(Nx, Ny); sy_p(:,1:end-1) = sigma_yp;
    sy_m = zeros(Nx, Ny); sy_m(:,2:end)   = sigma_yp;

    % Main diagonal
    main_diag = (sx_p + sx_m)/dx^2 + (sy_p + sy_m)/dy^2;

    % Get linear indices
    idx = reshape(1:N, Nx, Ny);
    I = []; J = []; S = [];

    % Center
    I = [I; idx(:)];
    J = [J; idx(:)];
    S = [S; main_diag(:)];

    % x-direction neighbors
    I = [I; reshape(idx(2:end,:),[],1)];
    J = [J; reshape(idx(1:end-1,:),[],1)];
    S = [S; -reshape(sx_m(2:end,:),[],1)/dx^2];

    I = [I; reshape(idx(1:end-1,:),[],1)];
    J = [J; reshape(idx(2:end,:),[],1)];
    S = [S; -reshape(sx_p(1:end-1,:),[],1)/dx^2];

    % y-direction neighbors
    I = [I; reshape(idx(:,2:end),[],1)];
    J = [J; reshape(idx(:,1:end-1),[],1)];
    S = [S; -reshape(sy_m(:,2:end),[],1)/dy^2];

    I = [I; reshape(idx(:,1:end-1),[],1)];
    J = [J; reshape(idx(:,2:end),[],1)];
    S = [S; -reshape(sy_p(:,1:end-1),[],1)/dy^2];

    % Build sparse matrix
    A = sparse(I, J, S, N, N);
end

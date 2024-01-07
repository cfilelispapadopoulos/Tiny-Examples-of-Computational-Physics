## Incompressible Navier - Stokes in Vorticity - Velocity formulation
Let us consider the Incompresible Navier - Stokes equation:
$$\frac{\vartheta \vec{u}}{\vartheta t}=-\nabla P + \frac{1}{Re} \Delta \vec{u} - (\vec{u} \cdot \nabla)  \vec{u},  (1)$$
$$\nabla \cdot u = 0,  (2)$$
where $P$ is the pressure, $Re$ is the [Reynolds number](https://en.wikipedia.org/wiki/Reynolds_number) and $\vec{u} = (u_1,u_2)$ is the velocity vector with components $u_1(x,y,t)$ and $u_2(x,y,t)$ across the $x\in[-\ell_x,\ell_x]$ and $y\in[-\ell_y,\ell_y]$ directions respectively. In this particular example we will consider the domain to be double periodic.
Let us define the vorticity $\vec{w}=\nabla \times \vec{u}$ which is a vector quantity which is normal to the the plane of the flow. By taking the curl operator $(\nabla \times )$ of the conservation of momentum equation (1):
$$\nabla \times\frac{\vartheta \vec{u}}{\vartheta t}=-\nabla \times\nabla P + \nabla \times \left( \frac{1}{Re} \Delta \vec{u} \right) - \nabla \times \left( (\vec{u} \cdot \nabla)  \vec{u} \right)$$
and taking into account that:
$$-\nabla \times\nabla P = 0,$$
since the curl of a gradient of a scalar field leads to the zero vector field $\vec{0}$. Moreover, the term:
$$\nabla \times \left( \frac{1}{Re} \Delta \vec{u} \right) = \frac{1}{Re} \left( \nabla \times \Delta \vec{u} \right)=\frac{1}{Re} \Delta \left( \nabla \times \vec{u} \right)=\frac{1}{Re} \Delta \vec{w},$$
since $\nabla \times$ and $\Delta$ operators commute. Similarly:
$$- \nabla \times \left( (\vec{u} \cdot \nabla)  \vec{u} \right)= -(\vec{u} \cdot \nabla) \vec{w}$$
and
$$\nabla \times\frac{\vartheta \vec{u}}{\vartheta t}=\frac{\vartheta \vec{w}}{\vartheta t}.$$
Thus, equation (1) takes the form:
$$\frac{\vartheta \vec{w}}{\vartheta t}=\frac{1}{Re} \Delta \vec{w}-(\vec{u} \cdot \nabla) \vec{w},$$
which, in practice, expresses the conservation of angular momentum. This formulation is advantageous, since it implicitly satisfies the incompressibility (conservation of mass or equation (2)). We can compute the velocities through the stream function $\psi$:
$$\vec{u} = \nabla^\perp \psi = \left( \frac{\vartheta \psi}{\vartheta y},-\frac{\vartheta \psi}{\vartheta x} \right)=(u_1,u_2)$$
and taking the second derivative of the stream function $\psi$ we obtain:
$$\nabla^2 \psi = \frac{\vartheta^2 \psi}{\vartheta x^2}+\frac{\vartheta^2 \psi}{\vartheta y^2}=-\frac{\vartheta u_2}{\vartheta x}+\frac{\vartheta u_1}{\vartheta y}=-w$$
where $w$ is scalar such that $\vec{w} = w\hat{z}$, since the velocity vector is two dimensional. Using the stream function we can prove that incompressibility is satfisfied implicitly:
$$\nabla \cdot \vec{u} = \frac{\vartheta u_1}{\vartheta x}+\frac{\vartheta u_2}{\vartheta y}=\frac{\vartheta^2 \psi}{\vartheta x \vartheta y}-\frac{\vartheta^2 \psi}{\vartheta y \vartheta x}=0.$$
Solving this PDE can be easily performed using a pseudospectral approach for space discretization and a method such as Crank - Nicolson / Adams - Bashforth 2nd order (CNAB2) similar to the one followed for solving the Kuramoto - Sivashinsky or the Nikolaevskiy (Combustion folder). Let us denote the [Fourier transform](https://en.wikipedia.org/wiki/Fourier_transform) operator as $\mathcal{F}( )$, the inverse Fourier transform operator as $\mathcal{F}^{-1}( )$ and the following transformed quantities:
$$U_1 = \mathcal{F}(u_1)$$
$$U_2 = \mathcal{F}(u_2)$$
$$W = \mathcal{F}(w)$$
$$\mathcal{F} \left(\frac{\vartheta f }{\vartheta x}\right)=i k_x F$$
$$\mathcal{F} \left(\frac{\vartheta f }{\vartheta y}\right)=i k_y F$$
$$\mathcal{F} \left(\frac{\vartheta^2 f }{\vartheta x^2}\right)=- k_x^2 F$$
$$\mathcal{F} \left(\frac{\vartheta^2 f }{\vartheta y^2}\right)=-k_y^2 F$$
where $k_x=\frac{2\pi}{\ell_x} \xi_x,k_y=\frac{2\pi}{\ell_y}\xi_y$, with $\xi_x = -N_x/2+1,-N_x/2+2,...,N_x/2,\xi_y = -N_y/2+1,-N_y/2+2,...,N_y/2$, are the wavenumbers across $x$ and $y$ directions. The points per dimension are denoted as $N_x$ and $N_y$. Thus, the Vorticity - Velocity equation, takes the form:
$$\frac{\vartheta W}{\vartheta t}=-\frac{1}{Re} (k_x^2+k_y^2) W-\mathcal{F}((\vec{u}\cdot \nabla) w),$$
applying the Crank - Nicolson scheme in time and the Adams - Bashforth scheme for the nonlinear part we obtain:
$$\frac{W^{t+1}-W^t}{\delta t}=-\frac{k_x^2+k_y^2}{2Re}(W^{t+1}+W^t)-1.5\mathcal{F}(u_1^{t}\mathcal{F}^{-1}(ik_x W^t)+u_2^{t}\mathcal{F}^{-1}(ik_y W^t))+0.5\mathcal{F}(u_1^{t-1}\mathcal{F}^{-1}(ik_x W^{t-1})+u_2^{t-1}\mathcal{F}^{-1}(ik_y W^{t-1})),$$
or equivalently:
$$\left(\frac{1}{\delta t}+\frac{k_x^2+k_y^2}{2Re}\right)W^{t+1}=\left(\frac{1}{\delta t}-\frac{k_x^2+k_y^2}{2Re}\right)W^{t}+1.5\mathcal{N}(\vec{u}^t,W^t)-0.5\mathcal{N}(\vec{u}^{t-1},W^{t-1}), t=0,1,2,...   (3)$$
with:
$$\mathcal{N}(\vec{u}^t,W^t)=-\mathcal{F}(u_1^{t}\mathcal{F}^{-1}(ik_x W^t)+u_2^{t}\mathcal{F}^{-1}(ik_y W^t)).$$
In order to start the iteration an initial condition $w^0$ is required. Using this initial condition the initial velocities can be computed from the stream function solution at $t=0$:
$$\Psi^0 = (k_x^2+k_y^2)^{-1} W^0=(k_x^2+k_y^2)^{-1} \mathcal{F}(w^0),$$
$$U_1^0 = i k_y \Psi^0,$$
$$U_2^0 = -i k_x \Psi^0.$$
From equation (3) is evident that the velocity and vorticity are required at time $t=-1$, thus for the first time step we set $\vec{u}^{-1}=\vec{u}^0$ and $w^{-1}=w^{0}$ $(W^{-1}=W^{0})$. For this example we can set the initial condition using a uniform random number generator in the interval $[-1,1]$, while the size of the domain was set to $\ell_x=1,\ell_y=1$. The Reynolds number was set to $Re=1000$, the maximum time for the simulation was to $T_{max}=1000$ and the timestep to $\delta t = 0.1$. For $N_x=N_y=128$, the solution is of the following form:
![nsvv](https://github.com/cfilelispapadopoulos/Tiny-Examples-of-Computational-Physics/assets/137081674/72edff28-f584-4aef-b515-3cdd4558d742)

One major issue arising in this example, which often arises when the solution of the Poisson equation with a pseudospectral method is considered, is division by zero due to the zero wavenumbers $(k_x,k_y)=(0,0)$. This can be mitigated by explicitly setting the value of the constant term (k_x,k_y) to $1$, during the formation of the operator matrx $-(k_x^2+k_y^2)$. Another option is using [Cauchy's integral formula](https://en.wikipedia.org/wiki/Cauchy%27s_integral_formula):
$$f(\mathcal{L})=\frac{1}{2\pi i}\oint_\Gamma f(z)(z I-\mathcal{L})^{-1}dz,$$
where $I$ denotes the identity matrix and the path $\Gamma$ can be chosen arbitrarily as long as it encloses all the eigenvalues of $\mathcal{L}$. These integrals can be computed numerically using techniques such as the trapezoidal rule, [which converges exponentially](https://epubs.siam.org/doi/book/10.1137/1.9780898719598). Due to the choice of Fourier transform for space discretization the the operator $\mathcal{L}$ is diagonal, thus the integrals can be performed individually as:
$$f(\mathcal{L_{j,j}})=\frac{1}{2\pi i}\oint_{\Gamma_j} f(z)(z-\mathcal{L_{j,j}})^{-1}dz.$$
The path $\Gamma_j$ is chosen as the unit circle centered at a diagonal element $\mathcal{L_{j,j}}$ (integrand pole), $\Gamma_j=( \mathcal{L_{j,j}}+e^{2\pi i z},0 < z \leq 1 )$. Thus, the path integrals take the form:
$$f(\mathcal{L_{j,j}})=\int_0^1 f \left( \mathcal{L_{j,j}}+e^{2\pi i z} \right)dz \approx \frac{1}{M} \sum_{\ell=1}^M f\left( \mathcal{L_{j,j}}+e^{2\pi i \ell / M} \right),$$
where $\mathcal{L_{j,j}}$ in the case of the Poisson operator takes the form $\mathcal{L_{k_x,k_y}}=-(k_x^2+k_y^2).$ In practice, $M=16,32,64$ points are sufficient to achieve machine precision.

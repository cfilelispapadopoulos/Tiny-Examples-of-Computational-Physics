## Incompressible Navier Stokes Equation and the Kelvin - Helmholtz Instability

The [Kelvin - Helmholtz instability](https://en.wikipedia.org/wiki/Kelvin%E2%80%93Helmholtz_instability) arises due to velocity differences between two fluids or layers of the same fluid. This type of instability is visible in the atmospheres of planets like Jupiter or the Sun. In order to for this instability to arise we have to model a viscous fluid flow model based on the Incompresible Navier - Stokes equation:
$$\frac{\vartheta \vec{u}}{\vartheta t}=-\nabla P + \frac{1}{Re} \Delta \vec{u} - (\vec{u} \cdot \nabla)  \vec{u},  (1)$$
$$\nabla \cdot u = 0,  (2)$$
where $P$ is the pressure, $Re$ is the [Reynolds number](https://en.wikipedia.org/wiki/Reynolds_number) and $\vec{u} = (u_1,u_2)$ is the velocity vector with components $u_1(x,y,t)$ and $u_2(x,y,t)$ across the $x\in[0,\ell_x]$ and $y\in[0,\ell_y]$ directions respectively. The first equation is the conservation of momentum, while the second, acting as a constraint, is the concervation of mass. The constraint is derived from the mass continuity equation (Lagrangian description):
$$\frac{D\rho}{Dt}+\rho (\nabla \vec{u})=0,$$
by considering that along the flow line the density is constant in an incompressible fluid:
$$\frac{D\rho}{Dt}=0 \rightarrow \nabla \cdot u = 0$$
The operator $\frac{D}{Dt}$ is the [material or sustantive derivative](https://en.wikipedia.org/wiki/Material_derivative) and is defined as follows:
$$\frac{D}{Dt}=\frac{\vartheta}{\vartheta t}+\vec{u}\cdot \nabla$$
defined for any tensor field depending on position and time coordinates only. In order to solve the Incompressible Navier - Stokes equation we will adopt the approach presented by [Novak](https://www.equalsharepress.com/media/NMFSC.pdf). Initially, an Implicit - Explicit (IMEX) scheme, namely [Crank - Nicolson / Adams - Bashforth (CNAB)](https://www.sciencedirect.com/science/article/pii/S016892741730226X) approach is applied to the conservation of momentum equation:
$$\frac{\vartheta \vec{u}}{\vartheta{t}}=L(\vec{u})+N(\vec{u}),$$
where $L$ is the linear part and $N$ is the nonlinear part. The linear part is discretized implicitly (using Crank - Nicolson) with respect to time, while the nonlinear part is discretized explicitly (using Adams - Bashforth).
$$\frac{\vec{u}^{t+1}-\vec{u}^{t}}{\Delta t}=-\frac{1}{2}(\nabla P^{t+1}+\nabla P^t)+\frac{1}{2Re}(\Delta \vec{u}^{t+1} +\Delta \vec{u}^t)-\bigg(\frac{3}{2}N(\vec{u}^{t})-\frac{1}{2}N(\vec{u}^{t-1})\bigg), t=0,1,2,...$$
where $\Delta t = \frac{T}{n_t}$ with $T$ denoting the maximum time of simulation and $n_t$ the number of time steps. The quantities $N(\vec{u}^t)$ and $N(\vec{u}^{t-1})$ are denoted as $N^{t}$ and $N^{t-1}$, respectively.
The biggest issue with the above formulation is that the pressure at time step $t+1$ is not known explicitly. In order to mitigate this issue we can split the discretized conservation of momentum into two equations, using an intermediate solution denoted as $\vec{u}^\star$:
$$\frac{\vec{u}^\star-\vec{u}^{t}}{\Delta t}=-\frac{1}{2}\nabla P^{t}+\frac{1}{2Re}(\Delta\vec{u}^\star +\Delta \vec{u}^t)-\frac{3}{2}N^t+\frac{1}{2}N^{t-1},  (3)$$
$$\frac{\vec{u}^{t+1}-\vec{u}^\star}{\Delta t}=-\frac{1}{2}\nabla P^{t+1}+\frac{1}{2Re}(\Delta\vec{u}^{t+1} - \Delta \vec{u}^\star).  (4)$$
Substituting the operators with their discrete analogue, formed either by using a [Spectral approach](https://www.equalsharepress.com/media/NMFSC.pdf), namely FFT, or either [Finite Difference or Finite Element method](https://uk.mathworks.com/academia/books/computational-science-and-engineering-strang.html), we obtain: 
$$\big(\Delta t^{-1} - \frac{1}{2Re} D^2\big) \vec{u}^\star = -\frac{1}{2} D P^t+\big(\Delta t^{-1} + \frac{1}{2Re} D.^2\big)\vec{u}^t-\frac{3}{2}N^t+\frac{1}{2}N^{t-1},$$
$$\big(\Delta t^{-1} - \frac{1}{2Re} D^2\big)(\vec{u}^{t+1}-\vec{u}^\star)=-\frac{1}{2} D P^{t+1},$$
where the discrete operators are defined as $D\approx \nabla$ and $D^2 \approx \Delta$. It should be noted that a different symbol should be chosen for the discete velocities and pressures, however this would introduce unnecessary complexity, thus are left as they were. Let us also denote the two matrices in the above equations as follows:
$$A=\big(\Delta t^{-1} - \frac{1}{2Re} D^2\big),$$
$$B=\big(\Delta t^{-1} + \frac{1}{2Re} D^2\big).$$
Thus, equations (3) and (4) become:
$$A \vec{u}^\star = -\frac{1}{2} D P^t+B\vec{u}^t-\frac{3}{2}N^t+\frac{1}{2}N^{t-1}, (5)$$
$$A(\vec{u}^{t+1}-\vec{u}^\star)=-\frac{1}{2} D P^{t+1}.$$
The last equation for $t \rightarrow t-1$ becomes:
$$A(\vec{u}^{t}-\vec{u}^{\star-1})=-\frac{1}{2} D P^{t},$$
where $\vec{u}^{\star-1}$ denotes the intermediate solution at a previous time step. We can substitute $D P^{t}$ in equation (5):
$$A \vec{u}^\star = A(\vec{u}^{t}-\vec{u}^{\star-1})+B\vec{u}^t-\frac{3}{2}N^t+\frac{1}{2}N^{t-1},$$
or
$$\vec{u}^\star = \vec{u}^{t}-\vec{u}^{\star-1}+A^{-1}\big(B\vec{u}^t-\frac{3}{2}N^t+\frac{1}{2}N^{t-1}\big) (6).$$
Then solution at the next time step can be computed using equation (2). This equation requires that the divergence of the solution should be equal to zero. We can achieve that by a projection of the form,
$$\vec{u}^{t+1}=\vec{u}^\star- \nabla \Delta^{-1} (\nabla \cdot \vec{u}^\star). (7)$$
which is often referred in literature as [Chorin's projection method](https://en.wikipedia.org/wiki/Projection_method_(fluid_dynamics)). Technically, this projection removes the part of $\vec{u}^\star$ which is colinear to $\nabla$, thus rendering the solution $\vec{u}^{t+1}$ perpendicular to $\nabla$.
For this particular example the spatial discretization is going to be performed using a spectral approach through Fast Fourier Transform (FFT). Thus, the discrete differential operators are: $D_x\approx\mathcal{F}\left(\frac{\vartheta}{\vartheta x}\right)=i k_x$, $D_y\approx\mathcal{F}\left(\frac{\vartheta}{\vartheta y}\right)=i k_y$ and $D^2=(i k_x,i k_y)\cdot(i k_x,i k_y) = -k_x^2-k_y^2$, where $\mathcal{F}$ denotes the Fast Fourier Transform (FFT). Thus, equations (6) and (7) become:
$$U^\star_x = U^{t}_x-U^{\star-1}_x+A^{-1} \odot \big(B \odot U^t_x-\frac{3}{2}\mathcal{N_x}^t+\frac{1}{2}\mathcal{N_x}^{t-1}\big)$$
$$U^\star_y = U^{t}_y-U^{\star-1}_y+A^{-1} \odot \big(B \odot U^t_y-\frac{3}{2}\mathcal{N_x}^t+\frac{1}{2}\mathcal{N_y}^{t-1}\big)$$
with
$$U_x = \mathcal{F}(u_x), U_y = \mathcal{F}(u_y)$$
and
$$\mathcal{N}_x^t=\mathcal{F} \left( \mathcal{F}^{-1}(U_x^t) \odot \mathcal{F}^{-1}(D_x \odot U_x^t) + \mathcal{F}^{-1}(U_y^t) \odot \mathcal{F}^{-1}(D_y \odot U_x^t) \right),$$
$$\mathcal{N}_y^t=\mathcal{F} \left( \mathcal{F}^{-1}(U_y^t) \odot \mathcal{F}^{-1}(D_x \odot U_y^t) + \mathcal{F}^{-1}(U_x^t) \odot \mathcal{F}^{-1}(D_y \odot U_y^t) \right),$$
$$U^{t+1}_x=U^\star_x- D_x \odot D^{-2} \odot (D_x \odot U^\star_x+D_y \odot U^\star_y),$$
$$U^{t+1}_y=U^\star_y- D_y \odot D^{-2} \odot (D_x \odot U^\star_x+D_y \odot U^\star_y).$$
The frequencies corresponding to the two spatial dimensions are: $k_x=\frac{2\pi}{\ell_x}j_x$ and $k_y=\frac{2\pi}{\ell_y}j_y$ with $j_x=-n_x/2+1,...,n_x/2$ and $j_y=-n_y/2+1,...,n_y/2$, respectively. The number of intervals for the discretization is denoted as $n_x$ and $n_y$ for the $x$ and $y$ directions, respectively. It should be mentioned that due to the diagonal nature of the discrete derivative operators the inversion of matrix $A$ is performed by element.

# Avoiding divisions by zero
To compute the velocities accross $x$ and $y$ we have to apply at some point the operator $D^{-2}=(-k_x^2-k_y^2)^{-1}$, which for very small $k_x,k_y\approx 0$ caused by cancellations and round off errors or even the fact that we include zero frequencies, will cause the operator to be represented incorrectly. In order to fix that issue we can use the generalized [Cauchy's integral formula](https://en.wikipedia.org/wiki/Cauchy%27s_integral_formula), similarly to the case in (Fluid Dynamics\Phase Separation):
$$f(\mathcal{L_{j,j}})=\frac{1}{2\pi i}\oint_{\Gamma_j} f(z)(z-\mathcal{L_{j,j}})^{-1}dz=\int_0^1 f \left( \mathcal{L_{j,j}}+e^{2\pi i z} \right)dz \approx \frac{1}{M} \sum_{\ell=1}^M f\left( \mathcal{L_{j,j}}+e^{2\pi i \ell / M} \right),$$
with $\Gamma_j=( \mathcal{L_{j,j}}+e^{2\pi i z},0 < z \leq 1 )$ denoting the unit circle centered at each point in the operator $D^2$. The new operator resulting from the integration, denoted $D'^2$, can be inverted without issue, thus we can substitute it to the two formulas for computing the new velocities. It should be noted that in practive 16, 32 or 64 are sufficient to achieve accuracy close to machine precision, since the sum [converges exponentially to the value of the integral](https://people.maths.ox.ac.uk/trefethen/publication/PDF/2005_111.pdf).

# Visualization
Visualization of a flow can be performed using the vectors or the pressure, but the produced figures will be busy since for each point of the grid a vector is required. One of the most popular approaches is the use of the stream lines, which are lines that a particle would follow for a given velocity field at a specific time. Of course a change in velocities in time would result to a change in these lines. In our case, where the flow is incompressible, which leads to divergence-free velocity fields (solenoidal field), thus the stream lines are closed. In practice these lines are just contours of a stream function which can be defined as:
$$\nabla^\perp Q = \left(-\frac{\vartheta Q}{\vartheta y},\frac{\vartheta Q}{\vartheta x} \right) = \vec{u}=(u_1,u_2)$$
or
$$\nabla \times \nabla^\perp Q = \nabla \times \vec{u}=-\Delta Q=-\frac{\vartheta^2 Q}{\vartheta y^2}-\frac{\vartheta^2 Q}{\vartheta x^2}=\frac{\vartheta u_1}{\vartheta y}-\frac{\vartheta u_2}{\vartheta x}$$
or in discrete form:
$$Q = \mathcal{F}^{-1}(-D'^{-2}\odot(D_y \odot U_1-D_x \odot U_2)).$$
The existence of the flow function is guaranteed, since the flow is incompressible, and it satisfies the conservation of mass (aka incompressibility condition or continuity equation).
Another way to visualize the flow is by evolving a quantity (field or tracer) $Q$ in time using an advection equation of the form:
$$\frac{\vartheta Q}{\vartheta t}+\vec{u}\cdot \nabla Q=0,$$
using a scheme like the Lax-Wendroff method. The quantity $\vec{u}$ is the velocity field obtained by the solution of the Incompressible Navier - Stokes. This approach has been followed by [Novak](https://www.equalsharepress.com/media/NMFSC.pdf). The quantity $Q$ acts like a tracer, it evolves following the flow, and it can be used to visualize the flow. The quantity $Q$ denotes the distribution of these "particles" that we want to monitor. After application of Lax - Wendroff we have (quantities are in physical space):
$$Q^{n+1}=Q^n-\left(T_2 \left(Q^n,\frac{\Delta t}{h_y}u_2^n \right)+T_1 \left(Q^n,\frac{\Delta t}{h_x}u_1^n \right) \right)$$
with $h_x=\ell_x/n_x$, $h_y=\ell_y / n_y$ ($n_x, n_y$ denote the number of intervals accross each dimension) and
$$T_1(Q,s)=s \odot (Q_{:,j}-Q_{:,j+1})-\frac{1}{2} s\odot(1-s)\odot((Q_{:,j}-Q_{:,j+1})+(Q_{:,j}-Q_{:,j-1})),$$
$$T_2(Q,s)=s \odot (Q_{i,:}-Q_{i+1,:})-\frac{1}{2} s\odot(1-s)\odot((Q_{i,:}-Q_{i+1,:})+(Q_{i,:}-Q_{i-1,:})).$$
The differences between the lines and columns of $Q$ can be easily handle due to the doubly periodic domain. The initial conditions chosen for the 2D flow are adopted by [Novak](https://www.equalsharepress.com/media/NMFSC.pdf) ($\ell_x=2,\ell_y=2$):
$$u_1(x,y)=(1+0.5sin(\ell_x \pi x))(0.5+0.5tanh(10-20|1-2y/(\ell_y)|)),$$
$$u_2(x,y)=0.$$
Initially, the velocities across the $y$ axis are set to zero. Thus, only the velocities in the $x$ axis are considered. Thus, we can choose the initial values for the tracer particles equal to $Q^0(x,y)=0.5+0.5 tanh(10-20|1-2y/(\ell_y)|)$ in order to separate them (in value and in color) only with respect to their position in the $y$ axis. This choice creates lines, initially, only towards the $y$ direction. Solution of the Navier - Stokes equation modeling the Kelvin - Helmholtz instability in 2D at some time step (with $Re=1000$, $n_t=3000$, $n_x=256$, $n_y=256$):
![ns](https://github.com/cfilelispapadopoulos/Tiny-Examples-of-Computational-Physics/assets/137081674/3198a731-6175-426d-856e-334f9bdbd89f)

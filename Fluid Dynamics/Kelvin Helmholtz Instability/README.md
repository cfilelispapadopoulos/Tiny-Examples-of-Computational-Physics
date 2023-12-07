## Incompressible Navier Stokes Equation and the Kelvin - Helmholtz Instability

The [Kelvin - Helmholtz instability](https://en.wikipedia.org/wiki/Kelvin%E2%80%93Helmholtz_instability) arises due to velocity differences between two fluids or layers of the same fluid. This type of instability is visible in the atmospheres of planets like Jupiter or the Sun. In order to for this instability to arise we have to model a viscous fluid flow model based on the Incompresible Navier - Stokes equation:
$$\frac{\vartheta \vec{u}}{\vartheta t}=-\nabla P + \frac{1}{Re} \Delta \vec{u} - \vec{u} \cdot \nabla  \vec{u},  (1)$$
$$\nabla \cdot u = 0,  (2)$$
where $P$ is the pressure, $Re$ is the [Reynolds number](https://en.wikipedia.org/wiki/Reynolds_number) and $\vec{u} = (u_1,u_2)$ is the velocity vector with componenets $u_1(x,y,t)$ and $u_2(x,y,t)$ across the $x$ and $y$ directions respectively. The first equation is the conservation of momentum, while the second, acting as a constraint, is the concervation of mass. The constraint is derived from the mass continuity equation:
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
$$\frac{\vec{u}^\star-\vec{u}^{t}}{\Delta t}=-\frac{1}{2}\nabla P^{t}+\frac{1}{2Re}(\Delta\vec{u}^\star +\Delta \vec{u}^t)-\frac{3}{2}H^t+\frac{1}{2}H^{t-1},  (3)$$
$$\frac{\vec{u}^{t+1}-\vec{u}^\star}{\Delta t}=-\frac{1}{2}\nabla P^{t+1}+\frac{1}{2Re}(\Delta\vec{u}^{t+1} - \Delta \vec{u}^\star).  (4)$$
Substituting the operators with their discrete analogue, formed either by using a [Spectral approach](https://www.equalsharepress.com/media/NMFSC.pdf), namely FFT, or either [Finite Difference or Finite Element method](https://uk.mathworks.com/academia/books/computational-science-and-engineering-strang.html), we obtain: 
$$\big(\Delta t^{-1} - \frac{1}{2Re} D^2\big) \vec{u}^\star = -\frac{1}{2} D P^t+\big(\Delta t^{-1} + \frac{1}{2Re} D.^2\big)\vec{u}^t-\frac{3}{2}H^t+\frac{1}{2}H^{t-1},$$
$$\big(\Delta t^{-1} - \frac{1}{2Re} D^2\big)(\vec{u}^{t+1}-\vec{u}^\star)=-\frac{1}{2} D P^{t+1},$$
where the discrete operators are defined as $D\approx \nabla$ and $D^2 \approx \Delta$. It should be noted that a different symbol should be chosen for the discete velocities and pressures, however this would introduce unnecessary complexity, thus are left as they were. Let us also denote the two matrices in the above equations as follows:
$$A=\big(\Delta t^{-1} - \frac{1}{2Re} D^2\big),$$
$$B=\big(\Delta t^{-1} + \frac{1}{2Re} D^2\big).$$
Thus, equations (3) and (4) become:
$$A \vec{u}^\star = -\frac{1}{2} D P^t+B\vec{u}^t-\frac{3}{2}H^t+\frac{1}{2}H^{t-1},$$
$$A(\vec{u}^{t+1}-\vec{u}^\star)=-\frac{1}{2} D P^{t+1}.$$
The last equation for $t \rightarrow t-1$ becomes:
$$A(\vec{u}^{t}-\big(\vec{u}^\star \big)^{-1})=-\frac{1}{2} D P^{t},$$
where $\big(\vec{u}^\star \big)^{-1}$ denotes the intermediate solution at a previous time step.

## Incompressible Navier Stokes Equation and the Kelvin - Helmholtz Instability

The [Kelvin - Helmholtz instability](https://en.wikipedia.org/wiki/Kelvin%E2%80%93Helmholtz_instability) arises due to velocity differences between two fluids or layers of the same fluid. This type of instability is visible in the atmospheres of planets like Jupiter or the Sun. In order to for this instability to arise we have to model a viscous fluid flow model based on the Incompresible Navier - Stokes equation:
$$\frac{\vartheta \vec{u}}{\vartheta t}=-\nabla P + \frac{1}{Re} \Delta \vec{u} - \vec{u} \cdot \nabla  \vec{u},$$
$$\nabla \cdot u = 0,$$
where $P$ is the pressure, $Re$ is the [Reynolds number](https://en.wikipedia.org/wiki/Reynolds_number) and $\vec{u} = (u_1,u_2)$ is the velocity vector with componenets $u_1(x,y,t)$ and $u_2(x,y,t)$ across the $x$ and $y$ directions respectively. The first equation is the conservation of momentum, while the second, acting as a constraint, is the concervation of mass. The constraint is derived from the mass continuity equation:
$$\frac{D\rho}{Dt}+\rho (\nabla \vec{u})=0,$$
by considering that along the flow line the density is constant in an incompressible fluid:
$$\frac{D\rho}{Dt}=0 \rightarrow \nabla \cdot u = 0$$
The operator $\frac{D}{Dt}$ is the [material or sustantive derivative](https://en.wikipedia.org/wiki/Material_derivative) and is defined as follows:
$$\frac{D}{Dt}=\frac{\vartheta}{\vartheta t}+\vec{u}\cdot \nabla$$
defined for any tensor field depending on position and time coordinates only.

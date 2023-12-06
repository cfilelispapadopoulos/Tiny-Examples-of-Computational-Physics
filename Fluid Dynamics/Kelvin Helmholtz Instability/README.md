## Incompressible Navier Stokes Equation and the Kelvin - Helmholtz Instability

The [Kelvin - Helmholtz instability](https://en.wikipedia.org/wiki/Kelvin%E2%80%93Helmholtz_instability) arises due to velocity differences between two fluids or layers of the same fluid. This type of instability is visible in the atmospheres of planets like Jupiter or the Sun. In order to for this instability to arise we have to model a viscous fluid flow model based on the Incompresible Navier - Stokes equation:
$$\frac{\vartheta \vec{u}}{\vartheta t}=-\nabla P + \frac{1}{Re} \Delta \vec{u} - \vec{u} \cdot \nabla  \vec{u}$$,
where $P$ is the pressure, $Re$ is the [Reynolds number](https://en.wikipedia.org/wiki/Reynolds_number) and $\vec{u} = (u_1,u_2)$ is the velocity vector with componenets $u_1(x,y,t)$ and $u_2(x,y,t)$ across the $x$ and $y$ directions respectively.


## Incompressible Navier - Stokes in Vorticity - Velocity formulation
Let us consider the Incompresible Navier - Stokes equation:
$$\frac{\vartheta \vec{u}}{\vartheta t}=-\nabla P + \frac{1}{Re} \Delta \vec{u} - (\vec{u} \cdot \nabla)  \vec{u},  (1)$$
$$\nabla \cdot u = 0,  (2)$$
where $P$ is the pressure, $Re$ is the [Reynolds number](https://en.wikipedia.org/wiki/Reynolds_number) and $\vec{u} = (u_1,u_2)$ is the velocity vector with components $u_1(x,y,t)$ and $u_2(x,y,t)$ across the $x\in[-\ell_x/2,\ell_x/2]$ and $y\in[-\ell_y/2,\ell_y/2]$ directions respectively. In this particular example we will consider the domain to be double periodic.
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

# Light Ray-Tracing Near a Schwarzschild Black Hole
## Schwarzchild metric
In order to simulate the trajectoris of light near a Schwarzchild black hole we need the [Schwarzchild metric](https://en.wikipedia.org/wiki/Schwarzschild_metric):

$$ds^2 = c^2 d\tau^2= - \left( 1-\frac{2GM}{c^2 r} \right) dt^2 + \left( 1-\frac{2GM}{c^2 r} \right)^{-1}dr^2 +r^2d\theta^2+r^2 sin^2 \theta d\phi^2,$$

with $\tau$ denoting the [proper time](https://en.wikipedia.org/wiki/Proper_time) and is an exact solution to the [Einstein Field Equations (EFE)](https://en.wikipedia.org/wiki/Einstein_field_equations) that describes the gravitational field outside a spherically symmetic, non rotating, uncharged mass (Cosmological constant $\Lambda$ is also zero). The chosen metric follows the mostly positive sign convention $(-+++)$. The $G$ is the gravitational constant, $M$ is the mass of the body (spherical mass) and $c$ is the speed of light. It should be noted that the quantity:

$$r_s = \frac{2GM}{c^2},$$

is the Schwarzchild radius that corresponds to the radius of a sphere in flat space that has the same surface area as that of the event horizon of a non rotating, uncharged, spherically symmetric black hole (Schwarzchild black hole) of a given mass.

In order to simplify computations the [geometric units](https://en.wikipedia.org/wiki/Geometrized_unit_system) will be adopted, which translated to $G=c=1$. Thus,

$$ds^2 = - \left( 1-\frac{2M}{r} \right) dt^2 + \left( 1-\frac{2M}{r} \right)^{-1}dr^2 +r^2d\theta^2+r^2 sin^2 \theta d\phi^2.$$

In the chosen system the Schwarzchild radius (event horizon) becomes $r_s = 2M$. Moreover we consider the motion only in the equatorial plane $\left(\theta=\frac{\pi}{2}\right)$, thus the metric becomes:

$$ds^2 = - \left( 1-\frac{2M}{r} \right) dt^2 + \left( 1-\frac{2M}{r} \right)^{-1}dr^2 +r^2 d\phi^2.$$

The metric, without loss of generality, restricts the motion to a plane, resulting in a (2+1)D system.

## Light Rays and Null Geodesics

The light rays follow null geodesics (lightlike or null worldlines), while for massive particles $ds^2<0$ (timelike worldlines), so the spacetime interval satisfies:

$$ds^2 = 0,$$

and the tangent $4$-vector to the geodesics satisfies:

$$g_{\mu\nu} \frac{dx^\mu}{d\lambda} \frac{dx^\nu}{d\lambda}=0,$$

where $\lambda$ is an affine parameter along the geodesic.

## Constants of Motion

We can treat motion as a variational problem using the the following Lagrangian $\mathcal{L}$:

$$\mathcal{L} = \frac{1}{2} g_{\mu\nu} \dot{x}^\mu \dot{x}^\nu = \frac{1}{2} g_{\mu\nu} \frac{dx^\mu}{d\lambda} \frac{dx^\nu}{d\lambda},$$

where $\lambda$ is an affine paramter (like proper time for massive particles or any affine parameter for light). Substituting the Schwarzchild metric given in the previous section we obtain:

$$\mathcal{L}=\frac{1}{2} \left[ -\left(1-\frac{2M}{r} \right)\dot{t}^2 + \left( 1- \frac{2M}{r}\right)^{-1} \dot{r}^2 + r^2\dot{\phi}^2 \right].$$

For null geodesics, $ds^2=0\Rightarrow \mathcal{L}=0$.

[Noether’s theorem](https://en.wikipedia.org/wiki/Noether%27s_theorem) tells us that every symmetry of the metric leads to a conserved quantity along geodesics.

A. Time translation symmetry $\rightarrow$ Energy conservation
The Schwarzschild metric is independent of $t$:

$$\frac{\partial\mathcal{L}}{\partial t} = 0 \Rightarrow \ conserved \ momentun \ p_t.$$

The canonical momentum is computed as follows:

$$p_t = \frac{\partial \mathcal{L}}{\partial \dot{t}} = -\left( 1-\frac{2M}{r} \right) \dot{t},$$

and the conserved energy per unit mass (or per unit affine parameter for light) is as follows:

$$E=-p_t=\left(1-\frac{2M}{r} \right)\dot{t},$$

which is conserved.

B. Rotational Symmetry $\rightarrow$ Angular Momentum Conservation
Similarly, the metric is independent of $\phi$:

$$\frac{\partial \mathcal{L}}{\partial \phi}=0 \Rightarrow \ conserved \ momentum \ p_\phi.$$

$$\frac{\partial \mathcal{L}}{\partial \dot{\phi}} = r^2 \dot{\phi},$$

Which is the conserved angular momentum $L=p_\phi=r^2 \dot{\phi}$.

For light, we use the null condition $\mathcal{L}=0$:

$$ -\left(1-\frac{2M}{r} \right)\dot{t}^2 + \left( 1- \frac{2M}{r}\right)^{-1} \dot{r}^2 + r^2\dot{\phi}^2=0$$

and substituting the conserved quantities:

$$\dot{t}=\frac{E}{1-\frac{2M}{r}},\dot{\phi}=\frac{L}{r^2},$$

we obtain:

$$\dot{r}^2 + \left( 1 - \frac{2M}{r} \right)\frac{L^2}{r^2}=E^2,$$

which is an effective energy equation for radial motion. The quantity $V_{eff}(r)=\left(1-\frac{2M}{r} \right) \frac{L^2}{r^2}$ is called the effective potential. The interpretation of the constants are as follows:

1. $E=\left( 1-\frac{2M}{r} \right) \dot{t}$: Energy at infinity (redshift energy)
2. $L=r^2\dot{\phi}$: Angular momentum per unit affine parameter
3. $b=\frac{L}{E}$: Impact parameter for light (distance of closest approach in flat spacetime)

These constants arise naturally from the geodesic Lagrangian, via Noether's theorem applied to the spacetime symmetries. Using the Impact parameter definition and a substitution of the form $u=\frac{1}{r}$ we obtain:

$$\left( \frac{du}{d\phi} \right)^2 + u^2=\frac{1}{b^2}+2Mu^3.$$

Differentiating with respect to $\phi$ and simplifying:

$$\frac{d^2 u}{d\phi^2}+u=3Mu^2$$

leads to the Orbit equation. This is the key equation we use for simulating light bending — it directly follows from the conservation laws. However, this equation cannot be used to simulate the full light ray path with natural initial contditions. Thus, the [Cartesian Geodesic Equations](https://en.wikipedia.org/wiki/Geodesics_in_general_relativity).

## Cartesian Geodesic Equations
In order to simulate the full light ray path with natural initial conditions, such as the position $(x,y)$, where the light ray starts, we have to solve the first-order geodesic equations in Cartesian-like coordinates or Schwarzschild coordinates with affine parameter $\lambda$:

$$\frac{d^2 x^\mu}{d\lambda^2}+\Gamma^\mu_{\alpha \beta} \frac{dx^\alpha}{d\lambda}\frac{dx^\beta}{d\lambda}=0,$$

where $\Gamma^\mu_{\alpha \beta}$ are the [Christoffel symbols](https://en.wikipedia.org/wiki/Christoffel_symbols) of the Schwarzschild metric:

$$\Gamma^\mu_{\alpha\beta}=\frac{1}{2} g^{\mu \nu}\left(\partial_\alpha g_{\beta\nu} + \partial_\beta g_{\alpha\nu} -\partial_\nu g_{\alpha\beta}\right).$$

Using the geodesic equation we derive the following system of $2$nd order nonlinear Ordinary Differential Equations (ODEs) correspoding to the three coordinates $x^\mu = (t(\lambda),r(\lambda),\phi(\lambda))$:

$$\frac{d^2 t}{d \lambda^2} + \frac{2M}{r^2 \left(1 - \frac{2M}{r}\right)}  \frac{dt}{d\lambda} \frac{d r}{d\lambda} = 0,$$

$$\frac{d^2 r}{d\lambda^2} + \frac{M}{r^2}  \left(1 - \frac{2M}{r}\right) \left(\frac{d t}{d\lambda}\right)^2
            - \frac{M}{r^2 \left(1 - \frac{2M}{r} \right)} \left(\frac{dr}{d\lambda}\right)^2
            - r \left(1 - \frac{2M}{r} \right) \left(\frac{d\phi}{d\lambda} \right)^2 = 0,$$

$$\frac{d^2 \phi}{d\lambda^2}+\frac{2}{r} \frac{dr}{d\lambda} \frac{d\phi}{d\lambda}=0.$$

These equations correspond to the time component, radial component, azimuthial component.

## Numerical Solution
In order to compute the paths of light rays we have to solve the two nonlinear ODEs corresponding to $r$ and $\phi$. However, these equations are second order and have to be converted to first order to solve them. Thus, we set $v=\dot{r}$ and $w=\dot{\phi}$ and we have:

$$\frac{du}{d \lambda} + \frac{2M}{r^2 \left(1 - \frac{2M}{r}\right)}  u v = 0,$$

$$\frac{d v}{d\lambda} + \frac{M}{r^2}  \left(1 - \frac{2M}{r}\right) u^2
            - \frac{M}{r^2 \left(1 - \frac{2M}{r} \right)} v^2
            - r \left(1 - \frac{2M}{r} \right) w^2 = 0,$$

$$\frac{dw}{d\lambda}+\frac{2}{r} v w=0.$$

$$u = \frac{dt}{d\lambda}$$

$$v = \frac{dr}{d\lambda}$$

$$w = \frac{d\phi}{d\lambda}$$

The initial conditions for three out of the six first order nonlinear ODEs in the system are given as input to the code $(t_0,r_0,\phi_0)$. The start time can be set to $0$, while the other two are the initial position where the rays are emitted. The initial conditions for the $(u,v,w)$ variables require a differnt approach. The first and third stem from the corresponding convervation equaitons thus:

$$u_0=\frac{dt}{d\lambda} \Bigg|_{\lambda=0} = \frac{E}{1-\frac{2M}{r_0}},$$

$$w_0=\frac{d\phi}{d\lambda} \Bigg|_{\lambda=0} = \frac{b}{r_0^2}.$$

The initial condition for the second quantity comes from solving $ds^2=0$ and is equal to:

$$v_0=\frac{dr}{d\lambda} \Bigg|_{\lambda=0}= -\sqrt{E^2-\left(1-\frac{2M}{r_0}\right) \frac{L^2}{r^2_0}}.$$

The initial condition for the derivative of time $\dot{t}$ creates an issue in the stability for the solution of the system. The issue arises because the expression for $\frac{dt}{dλ}=\frac{E}{1−2M/r}$ becomes very large near the black hole, where the denominator approaches zero. This leads to numerical instability when integrating the geodesic equations, especially since the second-order system is nonlinear and sensitive to large derivatives. When you set $\frac{dt}{dλ}=0$, the time component effectively decouples from the equations, allowing the spatial part of the system to evolve stably and correctly trace the light paths. Since coordinate time is not needed to determine the spatial trajectory of light rays, it is often excluded from such simulations unless time delays or redshift effects are being studied. Including the time variable requires careful numerical handling or a reformulation using conserved quantities instead of raw second-order derivatives.

The solution is performed with [MATLAB's ode45 routine](https://www.mathworks.com/help/matlab/ref/ode45.html) with the aforementioned initial conditions. The result of the simulation for $20$ rays is given below:

![eh](https://github.com/user-attachments/assets/41c32ce4-1aff-44b1-a3bd-4ef38e7a98c6)

We can observe that for low values of $b$ the light bends significantly and cannot escape the black hole attraction.

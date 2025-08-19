# Light Ray-Tracing Near a Kerr Black Hole
## Kerr metric

The geometry of spacetime around a rotating mass $M$ with angular momentum $J$ is described using the Einstein Field Equations (EFE) and more specifically the [Kerr metric](https://en.wikipedia.org/wiki/Kerr_metric). In order to simulate the trajectoris of light near a Kerr black hole we need the Kerr metric. In [Boyer-Lindquist coordinates](https://en.wikipedia.org/wiki/Boyer%E2%80%93Lindquist_coordinates) it takes the following form:

$$ds^2 = - \left( 1- \frac{r_s r}{\Sigma} \right)c^2 dt - \frac{2 r_s \alpha r sin^2 \theta}{\Sigma} c dt d\phi+\frac{\Sigma}{\Delta} dr^2 + \Sigma d\theta^2 + \frac{A sin^2 \theta}{\Sigma} d\phi^2,$$

where $r_s=\frac{2GM}{c^2}$ is the [Schwarzchild radius](https://en.wikipedia.org/wiki/Schwarzschild_metric). The $r,\theta,\phi$ are expressed in the oblate [spheroidal coordinate system](https://en.wikipedia.org/wiki/Oblate_spheroidal_coordinates). The transforamtion of these coordinates to the cartesian coordinate system is performed using the following equaitons:

$$x=\sqrt{r^2+\alpha^2}sin\theta cos\phi,$$

$$y=\sqrt{r^2+\alpha^2}sin\theta sin\phi,$$

$$z=r cos\theta.$$

The constants $\Sigma$, $\Delta$ and $\alpha$ are defined as follows:

$$\alpha = \frac{J}{Mc},$$

$$\Sigma = r^2 + \alpha^2 cos^2 \theta,$$

$$\Delta = r^2-r_s r + \alpha^2,$$

$$A = \left(r^2+\alpha^2 \right)^2-\alpha^2 \Delta sin^2 \theta,$$

with $J$ denoting the angular momentum and $M$ the rotating mass. It should be noted that $\alpha$ denotes the angular momentum per unit mass. The Kerr metric can be expressed in tensor form as follows $(t,r,\theta,\phi)$:

$$g_{\mu\nu}=\begin{bmatrix}g_{tt} & 0 & 0 & g_{t\phi} \\\ 0 & g_{rr} & 0 & 0 \\\ 0 & 0 & g_{\theta\theta} & 0 \\\ g_{t\phi} & 0 & 0 & g_{\phi\phi}\end{bmatrix},$$

with:

$$g_{tt}=-\left(1 - \frac{r_s r}{\Sigma} c^2 \right),$$

$$g_{t\phi}=-\frac{r_s \alpha r sin^2 \theta}{\Sigma} c,$$

$$g_{rr}=\frac{\Sigma}{\Delta},$$

$$g_{\theta\theta} = \Sigma,$$

$$g_{\phi\phi} = \frac{A sin^2 \theta }{\Sigma}.$$

The inverse Kerr metric $g^{\mu\nu}$ has the following form:

$$g^{\mu\nu}=\begin{bmatrix} -\frac{(r^2+\alpha^2)^2-\alpha^2 \Delta sin^2 \theta}{c^2 \Delta \Sigma} & 0 & 0 & -\frac{\alpha(r^2+\alpha^2-\Delta)}{c\Delta \Sigma} \\\ 0 & \frac{\Delta}{\Sigma} & 0 & 0 \\\ 0 & 0 & \frac{1}{\Sigma} & 0 \\\ -\frac{\alpha(r^2+\alpha^2-\Delta)}{c\Delta \Sigma} & 0 & 0 & \frac{\Delta - \alpha^2 sin^2 \theta}{\Delta \Sigma sin^2 \theta} \end{bmatrix}$$

In order to reduce complexity the geometric units are adopted, thus $G=c=1$.

## Constants of motion

In a Kerr black hole background, a particle (or photon) moving along a geodesic has several conserved quantities due to spacetime symmetries. The first one is energy $E$. Kerr spacetime is stationary, meaning the metric does not depend on the time coordinate $t$, thus:

$$\frac{\partial g_{\mu\nu}}{\partial t}=0.$$

This implies the existence of a timelike [Killing vector](https://en.wikipedia.org/wiki/Killing_vector_field):

$$\xi_{(t)}^\mu = (1,0,0,0).$$

The corresponding conserved quantity along a geodesic is the energy:

$$E=-p_\mu \xi_{(t)}^\mu=-p_t,$$

where $p_\mu = g_{\mu\nu}\frac{dx^\nu}{d\lambda}$ is the covariant $4$-momentum, and $\lambda$ is an affine parameter along the geodesic. The minus sign appears because the metric signature is $(-+++)$. Another conserved quantity is the Angular Momentum $L_z$. Kerr spacetime is axisymmetric, meaning the metric does not depend on the azimuthal angle $\phi$:

$$\frac{\partial g_{\mu\nu}}{\partial \phi}=0.$$

This gives a rotational Killing vector:

$$\xi_{(\phi)}^\mu = (0,0,0,1).$$

The conserved z-component of angular momentum is:

$$L_z = p_\mu \xi_{(\phi)}^\mu = p_{\phi}.$$

However, the motion of a particle (or photon) in Kerr spacetime is not confined to a plane (unlike Schwarzschild). To fully integrate the geodesic equations ($4$ second order differential equations $\rightarrow$ $4$ constants), we need a third constant. This is the Carter constant $Q$. Let us cosnider the [Hamilton-Jacobi equation](https://en.wikipedia.org/wiki/Hamilton%E2%80%93Jacobi_equation) which can be used to derive the geodesic motion:

$$\frac{\partial S}{\partial \lambda} = \frac{1}{2} g^{\mu\nu} \frac{\partial S}{\partial x^\mu}\frac{\partial S}{\partial x^\nu}=0.$$

For Kerr we can assume an additive separation ansatz of the form:

$$S=\frac{1}{2} m^2 \lambda - Et + L_z \phi S_r(r)+S_\theta (\theta),$$

where $m$ is particle mass (0 for photons). Plugging into Hamilton-Jacobi, after separation, gives two equations:

<b>Radial Part</b>

$$\Sigma^2 \left( \frac{dr}{d\lambda} \right)^2 = \left\[(r^2+\alpha^2)E-\alpha L_z \right\]^2-\Delta \left[Q+(L_z-\alpha E)^2 + m^2 r^2 \right],$$

<b>Polar Part</b>

$$\Sigma^2 \left(\frac{d\theta}{d\lambda}\right)^2 = Q-cos^2\theta \left\[ \alpha^2 (m^2-E^2) + \frac{L_z^2}{sin^2\theta} \right\],$$

with $\Sigma=r^2+\alpha^2 cos^2\theta$ and $\Delta = r^2 - 2Mr+\alpha^2$. . From the equation corresponding to coordinate $\theta$:

$$Q=p_\theta^2+cos^2\theta \left\[ \alpha^2 (m^2-E^2) + \frac{L_z^2}{sin^2\theta} \right\],$$

where $p_\theta = \Sigma \frac{d\theta}{d\lambda}$ and $Q\geq 0$ ensures motion remains real in $\theta$ since it guarantees both the radial equation and the polar equation have real solutions, i.e. the square roots are not imaginary. The Carter constant $Q$ measures the part of a particle’s angular momentum that is not aligned with the axis of rotation of the black hole. In Schwarzschild spacetime $(\alpha=0)$, the spacetime is spherically symmetric, so the total angular momentum $L^2$ is conserved. In Kerr, spherical symmetry is broken, leaving only axial symmetry; $L_z$ is conserved but not the full $L^2$. The hidden symmetry encoded by Carter’s constant effectively restores this missing piece: $Q$ generalizes $L^2-L_z^2$, quantifying the motion in the polar direction $(\theta)$ and ensuring the orbit can tilt out of the equatorial plane. Thus, if $Q$ measures the motion "out of the equatorial plane". If $Q=0$ motion is confined to the equatorial plane and larger $Q$ leads to more motion in the $\theta$ direction. 

The three constants are not enough to render the system integrable. Thus, another one is required. This one is derived from the mass-shell condition. 

## Mass-Shell condition
The mass-shell condition (normalization of the $4$-momentum/velocity) is the following:

$$g^{\mu\nu} p_\mu p_\nu = -\mu^2,$$

where $\mu$ is the particle rest mass. For massive particles $\mu=1$ (choosing $c=1$ and $\lambda = \tau$), while for massless (photons) $\mu=0$ (with any affine parameter $\lambda$). Substituting the known conserved quatities $p_t = -E$, $p_\phi = L_z$ and $p_\theta$ throuth the Carter constant condition we derive a condition for $p_r$:

$$\Delta p_r^2+p_\theta^2 + \frac{A}{\Delta}E^2-\frac{4M\alpha r}{\Delta}E L_z+\frac{\Delta-\alpha^2 sin^2\theta}{\Delta sin^2\theta}=-\mu^2 \Sigma.$$

This is the "master equation" relating the radial and polar momenta $p_r$, $p_\theta$ to the constants of motion.

## Equations of motion and geodesics
Let $x^\mu = (t,r,\theta,\phi)$, affine parameter $\lambda$, covariant momenta $p_\mu$, and inverse metric $g^{\mu\nu}(r,\theta)$ we can define the Hamiltonian as follows:

$$H=\frac{1}{2}g^{\mu\nu} p_\mu p_\nu,$$

with the mass-shell constraint $H=-\frac{1}{2}\mu^2$, which for photons becomes $H=0$ since $\mu=0$. From Hamiltons equations we have:

$$\frac{dx^\mu}{d\lambda}=\frac{\partial H}{\partial p_\mu}=g^{\mu\nu} p_\nu$$
$$\frac{dp_\mu}{d\lambda}=-\frac{\partial H}{\partial x^\mu}=-\frac{1}{2} \partial_\mu g^{\alpha \beta} p_\alpha p_\beta,$$

which lead to:

$$\frac{dt}{d\lambda}=g^{tt}p_t+g^{t\phi}p_\phi$$
$$\frac{dr}{d\lambda}=g^{rr}p_r$$
$$\frac{d\theta}{d\lambda}=g^{\theta\theta}p_\theta$$
$$\frac{d\phi}{d\lambda}=g^{t\phi}p_t+g^{\phi\phi}p_\phi$$
$$\frac{dp_r}{d\lambda}=-\frac{1}{2} \partial_r g^{\alpha\beta}p_\alpha p_\beta$$
$$\frac{dp_\theta}{d\lambda}=-\frac{1}{2} \partial_\theta g^{\alpha\beta}p_\alpha p_\beta$$
$$\frac{dp_t}{d\lambda}=0$$
$$\frac{dp_\phi}{d\lambda}=0$$

where the last two show the conserved quantities are constants of motion (stationarity and axial symmetry).

In order to determine the allowed initial momenta the null constraint is required. For a photon we have:

$$g^{\mu\nu}p_\mu p_\nu = 0$$

or

$$g^{tt}p_t^2+2g^{t\phi}p_t p_\phi + g^{\phi\phi} p_\phi^2 + g^{rr}p_r^2+g^{\theta\theta}p_\theta^2=0$$

and because of the conserved quantities:

$$p_r=\pm \sqrt{-\frac{A}{g^{rr}}}$$

with $A=g^{tt}p_t^2+2g^{t\phi}p_t p_\phi + g^{\phi\phi} p_\phi^2$, which helps initialize the radial momentum and ensures the trajectory really is a photon’s trajectory.

It should be noted that derivatives of the inverse metric can be computed using central finite differences of the the form:

$$\frac{\partial g^{\alpha\beta}}{\partial r} = \frac{g^{\alpha\beta}(r+\delta,\theta)-g^{\alpha\beta}(r-\delta,\theta)}{2\delta}+O(\delta^2),$$

which avoids writing extremely long symbolic derivatives, while still being accurate enough for numerical integration. The numerical integrations required can be performed using [ode45 routine in Matlab](https://www.mathworks.com/help/matlab/ref/ode45.html).

## Event horizon and ergosphere



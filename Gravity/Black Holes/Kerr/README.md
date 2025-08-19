# Light Ray-Tracing Near a Kerr Black Hole
## Kerr metric

The geometry of spacetime around a rotating mass $M$ with angular momentum $J$ is described using the Einstein Field Equations (EFE) and more specifically the [Kerr metric](https://en.wikipedia.org/wiki/Kerr_metric). In order to simulate the trajectoris of light near a Kerr black hole we need the Kerr metric. In [Boyer-Lindquist coordinates](https://en.wikipedia.org/wiki/Boyer%E2%80%93Lindquist_coordinates) it takes the following form:

$$ds^2 = -c^2 d\tau^2 = \left( 1- \frac{r_s r}{\Sigma} \right) +\frac{\Sigma}{\Delta} dr^2 + \Sigma d\theta^2 + \left( r^2 + \alpha^2 + \frac{r_s r \alpha^2}{\Sigma} \right) sin^2 \theta d\phi^2-\frac{2r_s r \alpha sin^2 \theta}{\Sigma} c dt d\phi,$$

where $r_s=\frac{2GM}{c^2}$ is the [Schwarzchild radius](https://en.wikipedia.org/wiki/Schwarzschild_metric). The $r,\theta,\phi$ are expressed in the oblate [spheroidal coordinate system](https://en.wikipedia.org/wiki/Oblate_spheroidal_coordinates). The transforamtion of these coordinates to the cartesian coordinate system is performed using the following equaitons:

$$x=\sqrt{r^2+\alpha^2}sin\theta cos\phi,$$

$$y=\sqrt{r^2+\alpha^2}sin\theta sin\phi,$$

$$z=r cos\theta.$$

The constants $\Sigma$, $\Delta$ and $\alpha$ are defined as follows:

$$\alpha = \frac{J}{Mc},$$

$$\Sigma = r^2 + \alpha^2 cos^2 \theta,$$

$$\Delta = r^2-r_s r + \alpha^2$$

with $J$ denoting the angular momentum and $M$ the rotating mass. It should be noted that $\alpha$ denotes the angular momentum per unit mass. In order to reduce complexity the geometric units are adopted, thus $G=c=1$.

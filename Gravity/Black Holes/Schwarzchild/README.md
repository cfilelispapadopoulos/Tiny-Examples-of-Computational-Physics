# Full Ray-Tracing Near a Schwarzschild Black Hole
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

## Light rays and Null Geodesics

The light rays follow null geodesics (lightlike or null worldlines), while for massive particles $ds^2<0$ (timelike worldlines), so the spacetime intercal satisfies:

$$ds^2 = 0,$$

and the tangent $4$-vector to the geodesics satisfies:

$$g_{\mu\nu} \frac{dx^\mu}{d\lambda} \frac{dx^\nu}{d\lambda}=0,$$

where $\lambda$ is an affine parameter along the geodesic.

## Constants of Motion


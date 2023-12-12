## The Cahn - Hilliard equation
The [Cahn - Hilliard equation](https://en.wikipedia.org/wiki/Cahn%E2%80%93Hilliard_equation) describes the phase transition between the two components of a binary fluid. These components separate spontaneously forming pure consetrations of each component. The equation is described as:
$$\frac{\vartheta c}{\vartheta t}=\Delta \left( c^3 - c -\epsilon^2 \Delta c \right),$$
where $\epsilon$ controls the length of the transition regions between domains and $c$ denotes the concentration.
In order to solve the equation numerically a [convexity splitting method](https://onlinelibrary.wiley.com/doi/abs/10.1002/cnm.2597) can be applied in order to ensure energy stability by balancing the terms:
$$\frac{\vartheta c}{\vartheta t}= \left( -\epsilon^2 \Delta^2 + \alpha \Delta \right) c +  \Delta \left( c^3 - (1+\alpha) c  \right)$$
where $\alpha$ is a stabilizing parameter. It should be noted that the term $\Phi(c) = c^3-c-\epsilon^2 \Delta c$ is ofter referred to as the Chemical Potential and different choices have been explored in the literature, such as the [double-well, Lennard Jones or the logarithmic Flory-Huggins potential](https://www.sciencedirect.com/science/article/abs/pii/S0898122123000500). The above Partial Differential Equation is split to a linear and a non-linear part:
$$\frac{\vartheta c}{\vartheta t}=L(c) + N(c),$$
where the operator $L$ denotes the linear part that is going to be handled implicitly with respect to time, while the operator $N$ denotes the non-linear part that is going to be handled explicitly.
For equations with stiff linear part, such as the Cahn -  Hilliard equation, [Exponential Time Differencing (ETD)](https://en.wikipedia.org/wiki/Exponential_integrator) type integrators can be used. For the above equation (noting the solution as $u$ and its Fourier transform as $U$), using ETD formulation, we obtain:
$$U_{n+1}=e^{\mathcal{L}\delta t} U_n + e^{\mathcal{L} \delta t} \int_{0}^{\delta t} e^{-\mathcal{L}\tau} \mathcal{N}(u(t_n+\tau),t_n+\tau)) d\tau,$$
where $t_n=n\delta t.$ The linear part is propagated exactly, while the nonlinear part is propagated by a method such as 4th order Runge - Kutta, through the approximation of the in (integral) in the above equation.
[Cox and Matthews](https://www.sciencedirect.com/science/article/abs/pii/S0021999102969950) presented several methods based on this formulation, however, herewith we utilize the modified approach proposed by [Kassam and Trefethen](https://people.maths.ox.ac.uk/trefethen/publication/PDF/2005_111.pdf) and more specifically the Exponetial Time Differencing Runge - Kutta 4th order (ETDRK4). These integrators benefit significantly when coupled with spectral methods for space discretization (such as FFT), since the matrix $\mathcal{L}$ corresponding to the linear part would be diagonal enabling element wise application. Using the ETDRK4 the integration can be performed as follows:
$$U_{n+1}=e^{\mathcal{L}\delta t}U_{n}+f_u \mathcal{N}(u_n,t_n)+ 2 f_{ab} \left( \mathcal{N}(a_n,t_n+\delta t_n / 2) +\mathcal{N}(b_n,t_n+\delta t /2) \right) + f_c \mathcal{N}(c_n,t_n+\delta t)$$
with:
$$a_n=e^{\mathcal{L}\delta t /2} U_n+\mathcal{L}^{-1} \left( e^{\mathcal{L}\delta t /2} - I \right)\mathcal{N}(u_n,t_n)$$
$$b_n=e^{\mathcal{L}\delta t /2} U_n+\mathcal{L}^{-1} \left( e^{\mathcal{L}\delta t /2} - I \right)\mathcal{N}(a_n,t_n+\delta t /2)$$
$$c_n=e^{\mathcal{L}\delta t /2} a_n+\mathcal{L}^{-1} \left( e^{\mathcal{L}\delta t /2} - I \right)\left(2 \mathcal{N}(b_n,t_n+\delta t /2)-\mathcal{N}(u_n,t_n) \right)$$
and:
$$f_u=\delta t^{-2} \mathcal{L}^{-3} \left( -4 -\mathcal{L}\delta t +e^{\mathcal{L}\delta t} \left( 4-3\mathcal{L}\delta t +(\mathcal{L}\delta t)^2 \right) \right)$$
$$f_{ab}=\delta t^{-2} \mathcal{L}^{-3} \left( 2+\mathcal{L}\delta t + e^{\mathcal{L}\delta t} (-2+\mathcal{L}\delta t) \right)$$
$$f_c=\delta t^{-2} \mathcal{L}^{-3} \left( -4 -3\mathcal{L}\delta t -(\mathcal{L}\delta t)^2 +e^{\mathcal{L}\delta t} (4-\mathcal{L}\delta t) \right),$$
where the term $\mathcal{N}(u_n,t_n)$ is equal to $\mathcal{F}(N(u_n,t_n))$, with $\mathcal{F}$ denoting the Discrete Fourier Transform (DFT).
This particular formulation though is subject [numerical instabilities](https://people.maths.ox.ac.uk/trefethen/publication/PDF/2005_111.pdf). This issue arises for ETD schemes of order greater than two, due to cancellation errors which are more apparent when the eigenvalues of $\mathcal{L}$ are close to zero. The main issue arises during the computation of the terms $f_u$, $f_{ab}$ and $f_c$. This has been tackled in the literature by techniques such as the use of cutoff value and substitution with [Taylor series approximation](https://www.sciencedirect.com/science/article/abs/pii/S0021999102969950). A more effective strategy was proposed form [Kassam and Trefethen](https://people.maths.ox.ac.uk/trefethen/publication/PDF/2005_111.pdf), which utilizes contour integrals based on [Cauchy's integral formula](https://en.wikipedia.org/wiki/Cauchy%27s_integral_formula):
$$f(\mathcal{L})=\frac{1}{2\pi i}\oint_\Gamma f(z)(z I-\mathcal{L})^{-1}dz,$$
where $I$ denotes the identity matrix and the path $\Gamma$ can be chosen arbitrarily as long as it encloses all the eigenvalues of $\mathcal{L}$. These integrals can be computed numerically using techniques such as the trapezoidal rule, [which converges exponentially](https://epubs.siam.org/doi/book/10.1137/1.9780898719598). Due to the choice of Fourier transform for space discretization the the operator $\mathcal{L}$ is diagonal, thus the integrals can be performed individually as:
$$f(\mathcal{L_{j,j}})=\frac{1}{2\pi i}\oint_{\Gamma_j} f(z)(z-\mathcal{L_{j,j}})^{-1}dz.$$
The path $\Gamma_j$ is chosen as the unit circle centered at a diagonal element $\mathcal{L_{j,j}}$ (integrand pole), $\Gamma_j=( \mathcal{L_{j,j}}+e^{2\pi i z},0 < z \leq 1 )$. Thus, the path integrals take the form:
$$f(\mathcal{L_{j,j}})=\int_0^1 f \left( \mathcal{L_{j,j}}+e^{2\pi i z} \right)dz \approx \frac{1}{M} \sum_{\ell=1}^M f\left( \mathcal{L_{j,j}}+e^{2\pi i \ell / M} \right).$$
In practice, values 16, 32 or 64 are sufficient to achieve accuracy close to machine precision. The integration has to be performed for 4 different functions corresponding to the three coefficients $f_u$, $f_{ab}$ and $f_c$ as well as the term $\mathcal{L}^{-1}(e^{\mathcal{L}\delta t/2}-I)$ which is required for the computation of the intermediate solution values $a_n$, $b_n$ and $c_n$.
For an initial condition $u_0$, drawn from a uniform random distribution in the interval $[-1,1]$, residing in the domain $[0,8]^2$ and a time interval $t\in[0,8]$ with $\delta t = 8/1024$, the solution of the Cahn-Hilliard PDE is as follows:

![cahn](https://github.com/cfilelispapadopoulos/Tiny-Examples-of-Computational-Physics/assets/137081674/ba128a49-344a-4fb8-8efa-8c042bdc653b)

The code can be naturally extended to 3D.

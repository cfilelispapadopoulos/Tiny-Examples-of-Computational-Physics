## The Cahn - Hilliard equation
The [Cahn - Hilliard equation](https://en.wikipedia.org/wiki/Cahn%E2%80%93Hilliard_equation) describes the phase transition between the two components of a binary fluid. These components separate spontaneously forming pure consetrations of each component. The equation is described as:
$$\frac{\vartheta c}{\vartheta t}=\Delta \left( c^3 - c -\epsilon^2 \Delta c \right)$$,
where $\epsilon$ controls the length of the transition regions between domains and $c$ denotes the concentration.
In order to solve the equation numerically a [convexity splitting method](https://onlinelibrary.wiley.com/doi/abs/10.1002/cnm.2597) can be applied in order to ensure energy stability by balancing the terms:
$$\frac{\vartheta c}{\vartheta t}= \left( -\epsilon^2 \Delta^2 + \alpha \Delta \right) c +  \Delta \left( c^3 - (1+\alpha) c  \right)$$
where $\alpha$ is a stabilizing parameter. It should be noted that the term $\Phi(c) = c^3-c-\epsilon^2 \Delta c$ is ofter referred to as the Chemical Potential and different choices have been explored in the literature, such as the [double-well, Lennard Jones or the logarithmic Flory-Huggins potential](https://www.sciencedirect.com/science/article/abs/pii/S0898122123000500). The above Partial Differential Equation is split to a linear and a non-linear part:
$$\frac{\vartheta c}{\vartheta t}=L(c) + N(c),$$
where the operator $L$ denotes the linear part that is going to be handled implicitly with respect to time, while the operator $N$ denotes the non-linear part that is going to be handled explicitly.
For equations with stiff linear part, such as the Cahn -  Hilliard equation, [Exponential Time Differencing (ETD)](https://en.wikipedia.org/wiki/Exponential_integrator) type integrators can be used. For the above equation (noting the solution as $u$ and its Fourier transform as $U$), using ETD formulation, we obtain:
$$U_{n+1}=e^{\mathcal{L}\delta t} U_n + e^{\mathcal{L} \delta t} \int_{0}^{\delta t} e^{-\mathcal{L}\delta t} \mathcal{N}(u(t_n+\tau),t_n+\tau)) d\tau.$$
[Cox and Matthews](https://www.sciencedirect.com/science/article/abs/pii/S0021999102969950) presented several methods based on this formulation, however, herewith we utilize the modified approach proposed by [Kassam and Trefethen](https://people.maths.ox.ac.uk/trefethen/publication/PDF/2005_111.pdf) and more specifically the Exponetial Time Differencing Runge - Kutta 4th order (ETDRK4). These integrators benefit significantly when coupled with spectral methods for space discretization (such as FFT), since the matrix $\mathcal{L}$ corresponding to the linear part would be diagonal enabling element wise application. Using the ETDRK4 the integration can be performed as follows:
$$U_{n+1}=e^{\mathcal{L}\delta t}U_{n}+f_u \mathcal{N}(u_n,t_n)+ 2 f_{ab} \left( \mathcal{N}(a_n,t_n+\delta t_n / 2) +\mathcal{N}(b_n,t_n+\delta t /2) \right) + f_c \mathcal{N}(c_n,t_n+\delta t)$$
with:
$$a_n=e^{\mathcal{L}\delta t /2} U_n+\mathcal{L}^{-1} \left( e^{\mathcal{L}\delta t /2} - I \right)\mathcal{N}(u_n,t_n)$$
$$b_n=e^{\mathcal{L}\delta t /2} U_n+\mathcal{L}^{-1} \left( e^{\mathcal{L}\delta t /2} - I \right)\mathcal{N}(a_n,t_n+\delta t /2)$$
$$c_n=e^{\mathcal{L}\delta t /2} a_n+\mathcal{L}^{-1} \left( e^{\mathcal{L}\delta t /2} - I \right)\left(2 \mathcal{N}(b_n,t_n+\delta t /2)-\mathcal(u_n,t_n) \right)$$
and:
$$f_u=\delta t^{-2} \mathcal{L}^{-3} \left( \right)$$
$$f_{ab}$$
$$f_c$$

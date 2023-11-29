## The Cahn - Hilliard equation
The [Cahn - Hilliard equation](https://en.wikipedia.org/wiki/Cahn%E2%80%93Hilliard_equation) describes the phase transition between the two components of a binary fluid. These components separate spontaneously forming pure consetrations of each component. The equation is described as:
$$\frac{\vartheta c}{\vartheta t}=\Delta \left( c^3 - c -\epsilon^2 \Delta c \right)$$,
where $\epsilon$ controls the length of the transition regions between domains and $c$ denotes the concentration.
In order to solve the equation numerically a [convexity splitting method](https://onlinelibrary.wiley.com/doi/abs/10.1002/cnm.2597) can be applied in order to ensure energy stability by balancing the terms:
$$\frac{\vartheta c}{\vartheta t}= \left( -\epsilon^2 \Delta + \alpha \Delta \right) c +  \Delta \left( c^3 - (1+\alpha) c  \right)$$
where $\alpha$ is a stabilizing parameter. It should be noted that the term $\Phi(c) = c^3-c-\epsilon^2 \Delta c$ is ofter referred to as the Chemical Potential and different choices have been explored in the literature, such as the [double-well, Lennard Jones or the logarithmic Flory-Huggins potential](https://www.sciencedirect.com/science/article/abs/pii/S0898122123000500). The above Partial Differential Equation is split to a linear and a non-linear part:
$$\frac{\vartheta c}{\vartheta t}=L(c) + N(c),$$
where the operator $L$ denotes the linear part that is going to be handled implicitly with respect to time, while the operator $N$ denotes the non-linear part that is going to be handled explicitly.

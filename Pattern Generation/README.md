## The Swift-Hohenberg equation
The Swift-Hohenberg equation, initially derived from [thermal convection equations](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.15.319), was proposed in order to describe patterns arising from phenomena in fluids such as the $Rayleigh-B\acute{e}nard$ convection. However, since then it has been utilized many fields spanning from economics to sociology and optics. The equation in its general form is as follows:
$$\frac{\vartheta u}{\vartheta t}=-(\Delta+1)^2 u +\epsilon u + K(u)$$
where $\epsilon$ is a bifurcation parameter and $N(u)$ is a smooth nonlinear function of $u$. In the current case $K(u)=-u^3$. In order to solve the above Partial Differential Equation (PDE) we will use the modified [Exponential Time Differencing / Runge - Kutta 4th order method](https://people.maths.ox.ac.uk/trefethen/publication/PDF/2005_111.pdf) (ETDRK4) similarly to the Cahn - Hilliard case (in folder Fluid Dynamics/Phase_Separation). Similarly space discretization will be performed using the Fourier Transform in a doubly periodic domain $[0,L]^2$. In order to proceed the PDE should be rewritten in the form:
$$\frac{\vartheta u}{\vartheta t}=L(u)+N(u)$$
with
$$L(u)=-(\Delta+1)^2$$ and $$N(u)=\epsilon u - u^3.$$
This splitting better describes the equilibrium dynamics ($\frac{\vartheta u}{\vartheta t}=0,t\rightarrow \infty$) and has been adopted from the [Kyle Novak's book on Scientific Computing](https://www.mathworks.com/academia/books/numerical-methods-for-scientific-computing-novak.html). The application of the modified ETDRK4 method is similar to the Cahn - Hilliard case (in folder Fluid Dynamics/Phase_Separation).

A graphical representation of the solution is as follows:
![sw](https://github.com/cfilelispapadopoulos/Tiny-Examples-of-Computational-Physics/assets/137081674/e023de19-3747-41f8-8f09-66dddcef18d0)

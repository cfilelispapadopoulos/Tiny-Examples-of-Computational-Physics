## The Kuramoto-Sivashinsky equation (or KS or flame equation)

The KS equation in one spatial dimension is described as follows:
$$\frac{\vartheta u}{\vartheta t}+\frac{1}{2} \left( \frac{\vartheta u}{\vartheta x} \right)^2 + \frac{\vartheta^2 u}{\vartheta x^2}+\frac{\vartheta^4 u}{\vartheta x^4}=0$$
with Periodic boundary conditions $u(-L/2,t)=u(L/2,t)$. The compact form of the above equation is:
$$u_t + \frac{1}{2} u_x^2+u_{xx}+u_{xxxx}=0,x \in [-L/2,L/2],t>0$$
and describes phenomena such as thermal instabilities in [flame fronts, trapped ion instabilities, etc](https://en.wikipedia.org/wiki/Kuramoto%E2%80%93Sivashinsky_equation). The solution of the equation leads to complex dynamical patterns and for large choices of $L$ to [bifurcations and chaotic states](https://encyclopediaofmath.org/wiki/Kuramoto-Sivashinsky_equation). Solving this equation numerically poses significant difficulties since it is a nonlinear fourth order differential equation. In order to numerically handle the equation, the differential form is preferable. The differential form is derived by differentiating the KS equation, with respect to $x$, and substituting $v=u_x$:
$$v_t + v_x v + v_{xx}+v_{xxxx}=0.$$
The term $v_{xx}$ creates instabilities is the source of instability in the large scales and the term $v_{xxxx}$ is a damping term for the small scales, while the convective term $v_x v$ is responsible for transferring energy between them in order to [stabilize the equation](https://encyclopediaofmath.org/wiki/Kuramoto-Sivashinsky_equation).
This is a stiff non-linear ODE and Implicit - Explicit schemes based on Crank-Nicolson/[Adams-Bashforth](https://en.wikiversity.org/wiki/Adams-Bashforth_and_Adams-Moulton_methods) ([CNAB2](https://core.ac.uk/download/pdf/12211007.pdf)) for temporal discretization and Fast Fourier Transform for spatial discretization have been considered. Other approaches such as [IMEX BDF](https://www.cs.uoi.gr/~akrivis/AS1.pdf) based schemes and [exponential time differentiation](https://people.maths.ox.ac.uk/trefethen/publication/PDF/2005_111.pdf) have been considered. Herewith, a simpler IMEX approach will be followed. Let us consider the following for the differential form of the KS equation:
$$v_t=Lv+N(v)$$
where $L=-\frac{\vartheta^2}{\vartheta x^2}-\frac{\vartheta^4}{\vartheta x^4}$ is the linear part and $N(v)=-v_x v$ is the non-linear part. In order to avoid linearization of the non-linear equation, we will proceed by discretizing the linear part implicitly and the non-linear part explicitly, using central Finite Difference approximation for all derivatives and taking into account periodic boundary conditions:
$$\frac{\vartheta v}{\vartheta x}=\frac{v_{i+1}-v_{i-1}}{2h}$$
$$\frac{\vartheta^2 v}{\vartheta x^2}=\frac{v_{i+1}-2v_{i}+v_{i-1}}{h^2}$$
$$\frac{\vartheta^4 v}{\vartheta x^4}=\frac{v_{i+2}-4v_{i+1}+6v_{i}-4v_{i-1}+v{i-2}}{h^4}$$
where $h=L/N$ is the mesh size and $N$ is the number of spatial intervals. We can arrange the equations into matrices ($N \times N$) $F$, $D$ and $D2=D^T D$ and write the KS equation as follows:
$$(I+\delta t D + \delta t D2 )v^{t+1}=v^t+\delta t N(v^t)$$,
with initial condition $v^0=cos(x)+0.15 cos(x/8)(1+2*sin(x/8))$ and $N(v^t)=-(Fv^t)\odot v^t$ ($\odot$: elementwise multiplication). The solution of the linear system is performed in a direct fashion at each timestep. It should be noted that the equations corresponding to last point of the grid ared not required to be retained in the linear system since the last point of the grid coincides with the first one. The presented approach is sensitive to the choices of time-step ($\delta t$) and mesh size ($h$).
Using the CNAB2 scheme the above equation would take the form:
$$(I+0.5 \delta t D + 0.5 \delta t D2 )v^{t+1}=(I-0.5 \delta t D - 0.5 \delta t D2 )v^t+1.5 \delta t N(v^t) - 0.5 \delta t N(v^{t-1})$$.
with $v_{-1}=v_{0}$. The 1D solution yields the following result:
![combustion](https://github.com/cfilelispapadopoulos/Tiny-Examples-of-Computational-Physics/assets/137081674/ed40e15f-a61a-4663-bf86-30eae5767111)

The multi-dimensional equation is of the form:
$$v_t+\Delta v + \Delta^2 v+0.5|\nabla v|^2=0.$$
The solution of the 2D equation is performed through FFT in doubly periodic domain $[-L/2,L/2]^2$, utilizing the CNAB2 scheme described above. A plot of the solution in an intermediate timestep is the following:
![2d](https://github.com/cfilelispapadopoulos/Tiny-Examples-of-Computational-Physics/assets/137081674/e599d086-14a6-450e-80f2-9a400f4c2607)

In order to solve the 2D equation we need to transform the derivatives into the frequency space as follows:
$$v_t - (k_x^2+k_y^2) V_{k_x,k_y} +(k_x^2+k_y^2)^2 V_{k_x,k_y} + \mathcal{F}( 0.5 | (\mathcal{F}^{-1}(i k_x V_{k_x,k_y}))^2 + (\mathcal{F}^{-1}(i k_y V_{k_x,k_y}))^2|) = 0$$
where $\mathcal{F}$ is the Fourier transform, $\mathcal{F}^{-1}$ is the inverse Fourier Transform and $V_{k_x,k_y}=v{x,y}$. The equation can be recast to CNAB2 form:
$$(I-0.5 \delta t (k_x^2+k_y^2)+0.5 \delta t (k_x^2+k_y^2)^2) V_{k_x,k_y}^{t+1} = (I+0.5 \delta t (k_x^2+k_y^2)-0.5 \delta t (k_x^2+k_y^2)^2) V_{k_x,k_y}^{t} + 1.5 \delta t \mathcal{N}(V_{k_x,k_y}^{t}) - 0.5 \delta t \mathcal{N}(V_{k_x,k_y}^{t-1})$$

with
$$\mathcal{N}(V_{k_x,k_y}^t)=\mathcal{F}( 0.5 | (\mathcal{F}^{-1}(i k_x V_{k_x,k_y}^t))^2 + (\mathcal{F}^{-1}(i k_y V_{k_x,k_y}^t))^2|).$$

The term $\delta t$ is the time step.

## The Nikolaevskiy equation

The [Nikolaevskiy equation](https://core.ac.uk/download/pdf/56373642.pdf) models seismic waves (longitudal) transmitted in [viscoelastic media](https://academic.oup.com/ptp/article/106/2/315/1878465). The equation is sixth-order and in one space dimension is of the following form:

$$u_t+u\vartheta_x u = -\vartheta_x^2 \left(r-(1+\vartheta_x^2)^2 \right)u= (1-r) \vartheta_x^2 u +2 \vartheta_x^4 u +\vartheta_x^6 u, x\in [-L/2,L/2],t>0$$

where the parameter $r$ controls the instability. For values of $r \leq 0$ the solution of the equation tends to zero, while for $r>0$ gives rise to interesting dynamical patterns. In two space dimensions the equation is of the following form:

$$u_t+0.5|\nabla u|^2 = (1-r) \Delta u +2 \Delta^2 u +\Delta^3 u.$$

In contrast to the Kuramoto-Sivashinsky equation the dissipative terms are the second (larger scales) and sixth order derivatives (smaller scales), while the fourth order derivative term is the source of the instability. This equation is solved numerically following the same approach as the 2D Kuramoto-Sivashinsky using the FFT (in space) and CNAB2 (Crank-Nicolson/Adams-Bashforth) scheme.
A snapshot of the solution is as follows:

![untitled](https://github.com/cfilelispapadopoulos/Tiny-Examples-of-Computational-Physics/assets/137081674/9bb5a975-e25f-4cc6-896f-dbfd17531e3d)

The code naturally extends to 3D.

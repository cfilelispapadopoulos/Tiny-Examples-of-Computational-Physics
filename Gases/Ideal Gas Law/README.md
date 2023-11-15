## 3D parallepiped with internal barier filled with $O_2$ in both compartments

![image](https://github.com/cfilelispapadopoulos/Tiny-Examples-of-Computational-Physics/assets/137081674/55796589-9deb-4fa2-9697-0c54fce49bee)

The parallepiped is of height $h$, width $w$ and length $l$ (across x) and the barrier, of mass $M_b$, is located at position $p$ on $x$ axis. The code is modular for both compartments but the temperature should be prescribed in order to compute the pressure $P$ of each compartment at each step and from that the force on the barrier of known area $h w$.

Sample molecures are added just for visualization purposes. The positions of the molecules following the uniform random distribution. The velocities are computed using the normal distribution with mean $\mu=0$ and standard deviation $\sigma=v_{a,i}=\sqrt{\frac{k_B T}{m_i}}$ (per dimension), where $k_B$ is the Maxwell - Boltzmann distribution, $T$ is the temperature and $m_i$ is the mass of an individual particle:

$$\pi(v_a)dv_a \propto exp(-\frac{v_a^2}{2\sigma^2})dv.$$

The pressure in the two compartments is computed at each step using the Ideal gas law $PV=nRT$, since the volume of compartments changes when the barrier moves. The gas in the two compartments follows isothermal expansion and compression processes, since temperatures are fixed. The force on the barrier is calculated as $F=\Delta P A = (P_1-P_2) h w$, where $P_1$ is the pressure of the left compartment and $P_2$ the pressure of the right one.

The barrier oscilates since the expansion of one gas results to the compression of the other. It should be noted that the walls and the barrier are perfect thermal insulators, which is something that can be improved in the future.

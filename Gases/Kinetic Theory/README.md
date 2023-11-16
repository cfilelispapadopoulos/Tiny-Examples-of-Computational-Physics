## Gas Piston, filled with $O_2$, in 2D in vacuum supported by spring

![image](https://github.com/cfilelispapadopoulos/Tiny-Examples-of-Computational-Physics/assets/137081674/da23d7dc-9024-4845-9c87-4252400b0f3e)

The gas piston, with width $w$, height $h$ and initial position $p$ (across the $x$ axis), is filled with $O_2$, which are assigned random positions inside the piston following the uniform random distribution. The velocities are computed using the normal distribution with mean $\mu=0$ and standard deviation $\sigma=1$ multiplied by $v_{a,i}=\sqrt{\frac{k_B T}{m_i}}$ (per dimension), where $k_B$ is the Maxwell - Boltzmann distribution, $T$ is the temperature and $m_i$ is the mass of an individual particle:

$$\pi(v_a)dv_a \propto exp(-\frac{v_a^2}{2\sigma^2})dv$$

Because of the large number of particles, e.g. 2 mol of $O_2$ are 2 $N_A$ molecules, they cannot be simulated directly. Thus, the molecules were grouped into $N$ groups  with mass equal to $m_i N$ and velocity that of a single molecule. This creates an issue, since the particle (group) has significant (unrealistic) momentum due to the larger than real mass. To solve this issue the force is spread across the whole square retaining the group. The side $d$ of the square retaining the group is estimated through the density of the which can be computed by $\rho= n M / (h p)$ (M is the molar mass of $O_2$, i.e. 0.032 g). Thus, $d = \sqrt(n M / \rho)$ and the force $F=\sum m_i v_i^2 / p$ is scaled as $F d^2$. The acceleration is analogously scaled. Integration of the laws of potion for both the piston and the particles is performed using [Leapfrom Kick-Drift-Kick scheme](https://en.wikipedia.org/wiki/Leapfrog_integration).

It should be noted that the initial position $p_i$ is considered rest length for the spring and the force at which the spring reacts is computed with Hooke's law $F=-k(p-p_i)$.

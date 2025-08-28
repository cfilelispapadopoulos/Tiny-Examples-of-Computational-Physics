# TM modes on a rectilinear waveguide using FDTD (Finite Difference Time Domain)
## Description
Let us consider a 3D rectangular waveguide of the form given in the Figure below.

<img width="862" height="447" alt="image" src="https://github.com/user-attachments/assets/111f1c4f-de4f-4e63-b863-223110823da6" />

A Gaussian wideband pulse is transmitted across the x-axis, reflecting on the internal square obstacle at the center of the waveguide and propagating towards the $-x$ and $+x$ directions. The obstacle is made out Perfect Electric Conductor (PEC). The surrounding faces of the rectilinear waveguide are covered with Perfect Electric Conductor (PEC) which reflects incident waves. The two faces across the x-axis are considered to be absorbing boundaries. These absorbing boundaries are simulated using [Convolutional Perfectly Matched Layer (CPML)](https://onlinelibrary.wiley.com/doi/10.1002/1098-2760(20001205)27:5%3C334::AID-MOP14%3E3.0.CO;2-A) in order to minimize reflections

## Maxwell's equation
In order to simulate wave propagation in a 3D waveguide, the fundamental equations are Maxwell’s equations in free space (no material dispersion except PEC obstacles and PML absorbing layers):

$$\nabla \times E = -\mu_0 \frac{\partial H}{\partial t},$$

$$\nabla \times H = \epsilon_0 \frac{\partial E}{\partial t}+J_{source},$$

$$\nabla \cdot E = 0,$$

$$\nabla \cdot H = 0,$$

where $E=(E_x,E_y,E_z)$ is the electric field $(V/m)$, $H=(H_x,H_y,H_z)$ is the magnetic field $(A/m)$, $\epsilon_0$ is the permitivity of free space, $\mu_0$ is the permeability of free space and $J_{source}$ represents any impressed current source. In the presented example the source is a Gaussian pulse in $E_z$. Moreover, the Electric displacement field can be computed as $D=\epsilon E$ $(C/m^2)$ with $\epsilon = \epsilon_r \epsilon_0$ ($\epsilon_r = 1$ and $\epsilon_0 = 8.854 \times 10^{-12} F/m$), and the Magnetic flux density is computed as $B=\mu H$ $(T)$ with $\mu = \mu_r \mu_0$ ($\mu_r = 1$ and $\mu_0 = 4\pi \times 10^{-7} H/m$). In our case due to free space $\epsilon_r=1$ and $\mu_r = 1$.

## Finite Differences and Yee's cell



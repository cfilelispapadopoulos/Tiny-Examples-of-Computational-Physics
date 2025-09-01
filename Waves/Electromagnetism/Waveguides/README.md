# TM modes on a rectilinear waveguide using FDTD (Finite Difference Time Domain)
## Description
Let us consider a 3D rectangular waveguide of the form given in the Figure below.

<img width="862" height="447" alt="image" src="https://github.com/user-attachments/assets/111f1c4f-de4f-4e63-b863-223110823da6" />

A Gaussian wideband pulse is transmitted across the x-axis, reflecting on the internal square obstacle at the center of the waveguide and propagating towards the $-x$ and $+x$ directions. The obstacle is made out Perfect Electric Conductor (PEC). The surrounding faces of the rectilinear waveguide are covered with Perfect Electric Conductor (PEC) which reflects incident waves. The two faces across the x-axis are considered to be absorbing boundaries. These absorbing boundaries are simulated using [Convolutional Perfectly Matched Layer (CPML)](https://onlinelibrary.wiley.com/doi/10.1002/1098-2760(20001205)27:5%3C334::AID-MOP14%3E3.0.CO;2-A) in order to minimize reflections. In this simulation TM (Transverse Magnetic) propagation is examined thus the $E_z$ component is plotted and $H_z=0$.

## Maxwell's equation
In order to simulate wave propagation in a 3D waveguide, the fundamental equations are Maxwell’s equations in free space (no material dispersion except PEC obstacles and PML absorbing layers):

$$\nabla \times E = -\mu_0 \frac{\partial H}{\partial t},$$

$$\nabla \times H = \epsilon_0 \frac{\partial E}{\partial t}+J_{source},$$

$$\nabla \cdot E = 0,$$

$$\nabla \cdot H = 0,$$

where $E=(E_x,E_y,E_z)$ is the electric field $(V/m)$, $H=(H_x,H_y,H_z)$ is the magnetic field $(A/m)$, $\epsilon_0$ is the permitivity of free space, $\mu_0$ is the permeability of free space and $J_{source}$ represents any impressed current source. In the presented example the source is a Gaussian pulse in $E_z$. Moreover, the Electric displacement field can be computed as $D=\epsilon E$ $(C/m^2)$ with $\epsilon = \epsilon_r \epsilon_0$ ($\epsilon_r = 1$ and $\epsilon_0 = 8.854 \times 10^{-12} F/m$), and the Magnetic flux density is computed as $B=\mu H$ $(T)$ with $\mu = \mu_r \mu_0$ ($\mu_r = 1$ and $\mu_0 = 4\pi \times 10^{-7} H/m$). In our case due to free space $\epsilon_r=1$ and $\mu_r = 1$.

## Finite Differences and Yee's cell
Numerical solution of Maxwell's equations for the described problem can be performed using finite difference schemes and more specifically [Yee's cell](https://en.wikipedia.org/wiki/Finite-difference_time-domain_method). The Yee's cell method relies on a staggered grid where the Electric Field resides in the mid points of the edges of a voxel while the magnetic field resides in the centers of the faces. Moreover, Yee's cell has temporal staggering where the Electric Field $E$ is updated at integer times steps $n\Delta t$ and the Magnetic Field is updated at half-integer time steps $(n+1/2)\Delta t$. This leapfrog scheme guarantees second-order accuracy and stability under the CFL condition. The method is applied to the Maxwell's curl equations:

$$\frac{\partial E}{\partial t} = - \frac{1}{\epsilon_r} \left( \nabla \times H \right)$$

$$\frac{\partial E}{\partial t} = - \frac{1}{\mu_r} \left( \nabla \times E \right)$$

and it results in the following update equations:

$$H_x^{n+1/2}(i,j,k)=H_x^{n-1/2}(i,j,k)-\frac{\Delta t}{\mu_0} \left( \frac{E_z^n (i,j+1,k)-E_z^n (i,j,k)}{\Delta y} - \frac{E_y^n (i,j,k+1)-E_y^n (i,j,k)}{\Delta z}\right)$$

$$H_y^{n+1/2}(i,j,k)=H_y^{n-1/2}(i,j,k)-\frac{\Delta t}{\mu_0} \left( \frac{E_x^n (i,j,k+1)-E_x^n (i,j,k)}{\Delta z} - \frac{E_z^n (i+1,j,k)-E_z^n (i,j,k)}{\Delta x}\right)$$

$$H_x^{n+1/2}(i,j,k)=H_x^{n-1/2}(i,j,k)-\frac{\Delta t}{\mu_0} \left( \frac{E_y^n (i+1,j,k)-E_y^n (i,j,k)}{\Delta x} - \frac{E_x^n (i,j+1,k)-E_x^n (i,j,k)}{\Delta y}\right)$$

$$E_x^{n+1}(i,j,k)=E_x^{n}(i,j,k)+\frac{\Delta t}{\epsilon_0} \left( \frac{H_z^{n+1/2}(i,j,k)-H_z^{n+1/2}(i,j-1,k)}{\Delta y} - \frac{H_y^{n+1/2}(i,j,k)-H_y^{n+1/2}(i,j,k-1)}{\Delta z} \right) + S_x$$

$$E_y^{n+1}(i,j,k)=E_y^{n}(i,j,k)+\frac{\Delta t}{\epsilon_0} \left( \frac{H_x^{n+1/2}(i,j,k)-H_x^{n+1/2}(i,j,k-1)}{\Delta z} - \frac{H_z^{n+1/2}(i,j,k)-H_z^{n+1/2}(i-1,j,k)}{\Delta x} \right) + S_y$$

$$E_z^{n+1}(i,j,k)=E_z^{n}(i,j,k)+\frac{\Delta t}{\epsilon_0} \left( \frac{H_y^{n+1/2}(i,j,k)-H_y^{n+1/2}(i-1,j,k)}{\Delta x} - \frac{H_x^{n+1/2}(i,j,k)-H_x^{n+1/2}(i,j-1,k)}{\Delta y} \right) + S_z$$

where $S_x$, $S_y$ and $S_z$ are the discrete source terms. 

## Boundary conditions

The Perfect Electric Conductor (PEC) Boundaries implementation is straight forward:

$$E_{tangential} = 0$$

on the wall surfaces. Thus:

$$E_x | PEC = E_y | PEC = E_z | PEC = 0,$$

which enforces reflection at the waveguide and obstacle walls. 

CPML is a special absorbing boundary to prevent reflections at open boundaries along $x$ axis. The Split-field approach is adopted for the CPML, thus for the $E_x$ in the $x$-CPML region we have:

$$E_x=E_{xy}+E_{xz}$$

where $E_{xy}$ accounts for the derivative along $y$ axis and $E_{xz}$ accounts for the derivative along $z$ axis. The discrete update equations for the split field are as follows:

$$E_{xy}^{n+1}=E_{xy}^{n}+\frac{\Delta t}{\epsilon_0} \frac{\partial H_z}{\partial y},$$

$$E_{xz}^{n+1}=E_{xz}^{n}-\frac{\Delta t}{\epsilon_0} \frac{\partial H_y}{\partial z}.$$

The CPML requires memory variables $\Psi$ that store the convolution of the field derivative with the damping profile:

$$\Psi_{E_{xy}}(t)=\int_0^t e^{-(\sigma_x/\kappa_x+\alpha_x)(t-\tau)/\epsilon_0} \frac{\partial H}{\partial y} (\tau) d\tau$$

or in discrete form:

$$\Psi_{E_{xy}}^{n+1}=b_x\Psi_{E_{xy}}^{n}+\frac{\Delta t}{\epsilon_0 \kappa_x} \left( \frac{\partial H_z}{\partial y} \right)^n$$

where $b_x=e^{-(\sigma_x / \kappa_x+\alpha_x)\Delta t / \epsilon_0}$, $\sigma_x$ is the CPML conductivity which controls absorption, $\kappa_x$ is the stretching factor to reduce reflection at grazing incidence, $\alpha_x$ is the low-frequency absorption coefficient. Thus, the total $E$ field in across $x$ in the CMPL region is:

$$E_x = E_{xy}+E_{xz}+\Psi_{E_{xy}}+\Psi_{E_{xz}}.$$

To minimize reflection, the CPML damping is graded gradually from $0$ to $\sigma_{max}$ over $N_{pml}$ cells. A quartic profile is often used:

$$\sigma_x (i) = \sigma_{max} \left( \frac{N_{pml}-i+0.5}{N_{pml}} \right)^4, i=1,...,N_{pml}.$$

The Quartic profile gives smooth absorption, reducing numerical reflection compared to linear or cubic. It should be noted that symmetric grading applies at both sides of the boundary. Similarly, the stretching factor $\kappa_x$ is graded:

$$\kappa_x(i)=1+(\kappa_{max}-1) \left( \frac{N_{pml}-i+0.5}{N_{pml}} \right)^4,$$

which increases the wave impedance gradually and helps the absorbtion of waves at grazing angles. It should be noted that polynomials of order $m$ can be used other than the fourth order utilized for this example. The low-frequency term $\alpha$ is often linearly graded:

$$\alpha_x(i)=\alpha_{max} \left( 1- \frac{N_{pml}-i+0.5}{N_{pml}} \right),$$

which enhances absorption of evanescent or low-frequency waves. In their discrete form the equations can be summarized as follows:

$$\Psi_{E_{xy}}^{n+1}=b_x \Psi_{E_{xy}}^n+\frac{\Delta t}{\epsilon_0 \kappa_x} \frac{H_z^n(i,j+1,k)-H_z^n (i,j,k)}{\Delta y},$$

$$\Psi_{E_{xz}}^{n+1}=b_x \Psi_{E_{xz}}^n-\frac{\Delta t}{\epsilon_0 \kappa_x} \frac{H_y^n(i,j,k+1)-H_y^n (i,j,k)}{\Delta y},$$

$$E_x^{n+1}=E_{xy}^{n+1}+E_{xz}^{n+1}+\Psi_{E_{xy}}^{n+1}+\Psi_{E_{xz}}^{n+1}.$$

Similar expressions exist for $E_y$ and $E_z$, splitting derivatives along x in the CPML region.

## Electromangetic Energy

The electromagnetic energy is a good indicator to assess how effective the CPML is and if the simulation is correct in many cases. The EM energy density at each point in space is:

$$u=\frac{1}{2} \left( \epsilon_0 |E|^2 + \mu_0 |H|^2 \right)$$

and the total energy in the computational domain is the sum over all cells:

$$U_{total} = \sum_{i,j,k} \frac{1}{2} \left( \epsilon_0 (E_x^2+E_y^2+E_z^2)+\mu_0 (H_x^2+H_y^2+H_z^2) \right)\Delta_x \Delta_y \Delta_z.$$

Different methods of integration can be implemented to increase the accuracy for computing the total energy, however, in this example we use the simplest one which is direct summation of the energy in all the voxels of the grid.

## Parameters

There are several parameters that affect the simulations such as the grid spacing $(\Delta x,\Delta y, \Delta z)$, the number of intervals $(N_x,N_y,N_z)$, the time step $\Delta t$ and the PML parameters $N_{pml}$, $\sigma_{max}$, $\kappa_{max}$ and $\alpha_{max}$.

The grid spacing and number of intervals can be chosen arbitrarily. However, the time step needs to fullfil the [CFL condition](https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition) of the following form:

$$\Delta t \leq  \frac{0.6}{c_0 \sqrt{\frac{1}{\Delta x^2}+\frac{1}{\Delta y^2}+\frac{1}{\Delta z^2}}},$$

where $c_0$ is the speed of light. The $N_{pml}$ which describes the number of intervals the PML consists of is usually chosen by experimentation. Under the a target one-pass reflection $R$ (power ratio, e.g. $R = 1e-6$ $\rightarrow$ $−60$ dB), the $\sigma_{max}$ can be estimated as:

$$\sigma_{max} = -\frac{(m+1) ln R}{2}\frac{\epsilon_0 c_0}{d_{pml}},$$

where $R=10^{-6}$, $m=4$ is the polynomial order, $d_{pml}=N_{pml}\Delta x$ is the thickness of the PML, $\epsilon_0$ is the permitivity of vacuum and $c_0$ is the speed of light. The typical values for $\alpha_{max}$ and $\kappa_{max}$ are $3$ and $0.05$, respectively.

An important check is to make sure that the dimensionless number $\beta = \frac{\sigma_{max}\Delta t}{\epsilon_0} is between $0.5-2.0$ in the outer edge, since $\beta>>1$ leads to a stiff discrete recusion.


## Results

The field component $E_z$ at different time steps is given below.
![1](https://github.com/user-attachments/assets/bd026189-e5cd-4f73-9059-241cb4099f4d)

![2](https://github.com/user-attachments/assets/72942d4a-76c5-4173-a140-08d1b0b5f62b)

![3](https://github.com/user-attachments/assets/bf598f7c-bb0f-40f1-af00-3a0ab8af57f9)

![4](https://github.com/user-attachments/assets/80f6c988-cea3-443f-a302-2b62b33f9cac)




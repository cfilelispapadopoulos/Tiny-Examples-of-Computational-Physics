# Hourglass simulation with Discrete Element Method

## Definition of the system
We consider a system of $N$ particles moving within an hourglass composed of lines in 2D under the influence of gravity and contact forces. The simulation is performed using the [Discrete Element Method (DEM)](https://en.wikipedia.org/wiki/Discrete_element_method).

Each particle $i$ is characterized by the following quantities:

Position: $x_i \in \mathbb{R}^2$,

Velocity: $v_i = \dot{x}_i$,

Angular Velocity: $\omega_i \in \mathbb{R}$,

Radius: $R_i$,

Mass: $m_i=\rho \pi R_i^2$,

Moment of Inertia for a disk: $I_i=\frac{1}{2} m_i R_i^2$,

Thus its particle is a disk interacting similarly to a rigid body with other particles and the wall.

## The equations of motion

The equations of motion for each particle $i$ should account for two types of motion, i.e. translational and rotational, since each particle is modeled as a disk of density $\rho$. The translational motion can be described using Newton's second law as follows:

$$m_i \frac{d^2 x_i}{dt^2}= F_i^{gravity}+\sum_{j} F_{ij}^{pp} + \sum_{w} F_{iw}^{pw},$$

where $F_i^{gravity}$ is the force of gravity on the particle, $F_{ij}^{pp}$ is the particle-particle contact force between body $i$ and $j$ and $F_{iw}^{pw}$ is the particle-wall contact force between particle $i$ and wall $w$. For the rotational motion using Newton's second law we obtain:

$$I_i \frac{d \omega_i}{d t}=\sum_j \tau_{ij}^{pp}+\sum_{w} \tau_{iw}^{pw},$$

where $\tau_{ij}^{pp}$ and $\tau_{iw}^{pw}$ are the torques due to tangential forces and rolling resistance between the $i$-th and $j$-th particle and the $i$-th particle and $w$-th wall, respectively.

## Contact detection

In order for the various forces and torques to be computed, contact between particles or between particles and walls should be detected. Two particles $i$ and $j$ are in contact if:

$$\delta_{ij}=R_i+R_j-||x_i-x_j||,$$

where $\delta_{ij}>0$ is the overlap. For the particle-wall contact, each wall is a line segment of the form:

$$s(t)=x_1+t(x_2-x_1), t\in[0,1],$$

thus the closest point $x_c$ on this wall and particle $i$ is:

$$t^{'}=clamp ( \frac{(x_i-x_1)(x_2-x_1)}{||x_2-x_1||^2},0,1),$$
$$x_c=x_1+t^{'} (x_2-x_1),$$

where $clamp(x,a,b)$ is a function that restricts $x$ in the interval $[a,b]$. Thus, if $x<a$ then $clamp(x,a,b)=a$, if $x>b$ $clamp(x,a,b)=b$ and for $a\leq x \leq b$ we get $clamp(x,a,b)=x$. The $clamp$ function ensures the closest point lies on the wall segment, not outside it, avoiding contact with non-existent geometry. The overlap with the wall is $\delta_{iw}=R_i-||x_i-x_c||$.

## Relative velocities at contact

The normal relative velocities are required in order to compute compressions or separations. While the tangential relative velocities are required to compute slidingand friction. The normal relative velocity is defined as:

$$\nu_n = (v_j-v_i)\dot n_{ij},$$

where $n_{ij}=\frac{x_j-x_i}{||x_j-x_i||}$ is a unit normal vector pointing from particle $i$ towards particle $j$. The tangential relative velocity, including rotation, is defined as follows:

$$\nu_t=(v_j-v_i)\dot t + R_i \omega_i + R_j \omega_j,$$

with $t=(-n_y,n_x)$. The $n_y$ and $n_x$ are the two components of the $n_{ij}$ vector.

## Normal contact force (soft-sphere model)

The normal contact force consists of two terms: the Hertzian elastic term and a linear viscous damping term:

$$F_n=k_n \delta^{3/2}-\gamma_n \nu_n,$$

where $k_n$ is the Hertzian normal stiffness, $\gamma_n$ is normal damping coefficient, $\delta_{ij}^{3/2}$ is the nonlinear Hertzian elastic response and $\nu_u$ is the normal relative velocity. The force is computed under the constraint $F_n\geq 0$ which clips the normal force in case its negative avoiding attracting forces. The physical interpretation of the first term of the force (elastic term) is that it represents elastic deformation at the contact with nonlinear stiffening with increasing overlap. Moreover, prevents large penetrations when $k_n$ is sufficiently large. The second term (damping term) acts only when particles move relative to each other, dissipated kinetic energy and controls restitution (bounciness). More importantly, damping is directional and it depends on relative motion, not absolute motion.

## Tangential contact force (Coulomb Friction) and rolling resistance

The elastic-damped tangential force can be described as follows:

$$F_t^{trial}=-k_t\nu_t,$$

where $k_t$ is a tangential stiffness/damping parameter. This is not a history-based spring (Cundall–Strack). It is a viscous friction approximation, chosen for simplicity and robustness. Real friction cannot grow arbitrarily large. It is limited by the normal force:

$$|F_t|\leq \mu F_n,$$

where $\mu$ is the friction coefficient and $F_n$ is the normal contact force magnitude. Thus the tangential force law takes the form:

$$F_t = \begin{cases}
F_t^{trial}, & |F_t^{trial}|\leq \mu F_n \\
-\mu F_n sign(\nu_t) & |F_t^{trial}| > \mu F_n
\end{cases},$$

which enforces sticking at low relative tangential velocity and sliding at the Coulomb limit when friction saturates.

Rolling friction introduces a resisting torque of the form:

$$\tau^{roll}_i=-\mu_r F_n R_i sign(\omega_i),$$

where $\mu_r$ is the rolling friction coefficient, $F_n$ is the normal contact force magnitude, $R_i$ is the $i$-th particle radius and $\omega_i$ is the angular velocity. The rolling friction torque stabilizes piles and prevents unrealistically free rolling. This torque always opposes rotation, it is proportional to contact load and ndependent of sliding velocity. The $sign$ term ensures $\tau_i^{roll} \omega_i \leq 0$ guarantees energy dissipation. Power dissipated by rolling resistance $P_{roll}=\tau_i^{roll} \omega_i=-\mu_r F_n R_i |\omega_i|$.

## Force and torque assembly

With the above in mind the total force on particle $i$ is:

$$F_i = \sum_j F_{ij}+\sum_{w} F_{iw}+m_i g,$$

with $F_{ij}=F_{ij}^{(n)}+F_{ij}^{(t)}$, where the superscript $(n)$ denotes the normal forces and $(t)$ the tangential forces. The total torque on particle $i$ is:

$$\tau_i=-\sum_{j} R_i F_{t,ij} - \sum_w R_i F_{t,iw}+ \tau^{roll}.$$

The equations guarantee linear momentum conservation, due to equal and opposite particle–particle forces, and angular momentum consistency, since torques arise only from tangential forces and rolling resistance.

## Time integration (Velocity Verlet)

The method chosen for time integration of the equations of motion is the [Velocity Verlet](https://en.wikipedia.org/wiki/Verlet_integration):

$$x^{n+1}_i=x_i^{n}+v_i^n \Delta t + \frac{1}{2} a_i^n \Delta t^2,$$

$$\omega_i^{n+1/2}=\omega_i^n+\frac{1}{2} \alpha_i \Delta t,$$

$$v_i^{n+1}=v_i^{n}+\frac{1}{2} (a_i^n+a_i^{n+1})\Delta t,$$

$$\omega_i^{n+1}=\omega_i^n+\frac{1}{2}(\alpha_i^n+\alpha_{i+1}^{n+1})\Delta t,$$

with $a_i=F_i/m_i$ and $\alpha_i=\tau_i/I_i$. The time step is chosen such that the [CFL](https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition) condition is satisfied:

$$\Delta t = CFL \; min_i(\sqrt{m_i/k_n}),$$

where $k_n$ is the Hertzian normal stiffness and $CFL$ is the CFL coefficient.

## Initial conditions and relaxation

In order to simulate the granular nature of sand, with different radius per particle, a log normal particle size distribution is utilized, since Many natural granular materials follow approximately log-normal statistics and skewed distributions are common in sands or powders:

$$R~LogNormal(\mu_{log},\sigma_{log},$$

where $\mu_{log}$ mean of the natural logarithm of radius and $\sigma_{log}$ is the standard deviation of the natural logarithm of radius.

Another important step is relaxation. When particles are randomly placed in the domain, they often overlap initially and large initial overlaps create unrealistic forces. Moreover, direct simulation may lead to instabilities or explosions. Thus, relaxation is a preprocessing step that removes initial overlaps, brings particles to a mechanically stable state and repares realistic initial velocities and forces. The relaxation step requires high damping, thus $\gamma_n$ is set to a large value $\gamma_{relax}$, which dissipates kinetic energy quickly and prevents particles from bouncing excessively. For a prescribed number of steps the Velocity Verlet method is used to relax the particles. The forces are computed with high damping and the system slowly relaxes toward equilibrium. The relaxation is typically run for a fixed number of time steps (e.g. 1000) or until the maximum particle velocity falls below a threshold:

$$max_i ||v_i||<\nu_{tol},$$
$$max_i ||\omega_i||<\omega_{tol},$$

and the system is considered mechanically stable. Once relaxation is complete the damping to normal simulation value $(\gamma_n)$.

Some snapshots of the simulation with the following parameters are given below:

$$R_{mean}=0.003 m$$

$$R_{std}=0.001 m$$

$$\rho=2500 kg/m^3$$

$$k_n=10^5 N/m$$

$$k_t=2\times 10^4 N/m$$

$$\gamma_n = 50 Ns/m$$

$$\gamma_{relax}=200 Ns/m$$

$$\mu=0.4$$

$$\mu_{roll}=0.02$$

It should be noted that neighboring points are discovered using the [KDtree algorithm](https://en.wikipedia.org/wiki/K-d_tree).

![2](https://github.com/user-attachments/assets/5288656e-e53d-464b-a613-4e722c15c38b)
![1](https://github.com/user-attachments/assets/75ae9ae7-5f7a-4b0d-a314-796a765027c2)


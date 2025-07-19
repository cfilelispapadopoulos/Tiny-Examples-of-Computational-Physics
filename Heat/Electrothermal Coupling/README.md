# Electrothermal Coupling

## 2D Electrothermal Simulation of an Integrated Circuit Interconnect

The temperature distribution of an IC interconnect network due to Joule heating is a complex phenomenon that couples the electrical potential and temperature fields. Thus, the couple solution of two Partial Differential Equations simultaneously is required. In order to showcase this complex phenomenon, initially, an model integrated circuit will be defined.

### Model Integrated Circuit

The considered integrated circuit is given in the Figure below.

<img width="515" height="301" alt="icm" src="https://github.com/user-attachments/assets/63dfadea-41b1-4452-8873-0ca92e7a6590" />

The IC is composed of three resistors, where the first one is connected, in series, to the other two which are connected in parallel. All three stainless steel resistors with [electrical conductivity](https://en.wikipedia.org/wiki/Electrical_resistivity_and_conductivity) $\sigma_R^0=1.45\times 10^{-6} (S/m)$ lie on a substrate ($0.2 m \times 0.2 m$) with electrical conductivity $\sigma_S = 10^{-6} (S/m)$. All components are connected with copper electrical wire (red) with electrical conductivity $\sigma_W=5.8\times 10^7 (S/m)$.

The conductivity of the wire as well as the substrate are considered to be independent of the temperature, while the conductivity of the resistors is considered to be varying with temperature:

$$\sigma_R(T) = \frac{\sigma_R^0}{1+\alpha (T-T_0)},$$

where $\alpha$ is the temperature coefficient ($K^{-1}$) (for steel is approximately $3 \times 10^{-3}$), $T$ is the current temperature of the resistor and $T_0$ is the reference temperature. It should be noted that the width of the wires and the length and width of resistors are significant and are considered parameters in the numerical modelling process that follows. In order to make the modelling process easier, the distance of the resistors and length of the cables is considered arbistrary, since it does not affect significantly the results of this example. The governing equations will be presented in the next part.

### Electrical Potential

Initially the electrical potential of the IC is required. The Electrical Potential can be computed using the following governing equation:

$$-\nabla \cdot (\sigma(T) \nabla T) = 0.$$ 

This is a linear PDE which can be solved by either the Finite Differences or the Finite Elements method. Due to the normal geometry of the substrate (square) the Finite Differences method will be considered. In order to apply the FD method the substrate is discretized using a square grid composed of $n_x+1$ points in the $x$-direction and $n_y+1$ points in the $y$-direction. 

![grid](https://github.com/user-attachments/assets/6f3463de-5021-4a31-9f1e-80e5e06259a0)

It should be noted that the axis are switched because MATLAB utilizes column-wise storage. Thus, the column wise indexing corresponds to the $x$ direction. Moreover, the number of intervals per dimensions is denoted as $n_x$ for the $x$-axis and $n_y$ for the $y$-axis.

In order to solve the PDE the regular five point stencil is not going to yield a physically accurate solution. Instead, a conservative Flux-based discretization will be chosen, since the aforementioned PDE is derived from the charge conservation:

$$-\nabla \cdot (\sigma(T) \nabla T) = 0 \Rightarrow \nabla \cdot J = 0,$$

Thus, the net current into each cell is exactly zero. Another important advantage of Flux conservative FD is that this type of discretization preserves continuity $\sigma_1 \frac{\partial V}{\partial n} = \sigma_2 \frac{\partial V}{\partial n}$ in case of material discontinuities, such as conductivity discontinuity, especially when combined with harmonic averaging:

$$\sigma_{harmonic} = \frac{2\sigma_1 \sigma_2}{\sigma_1 + \sigma_2},$$

in neighboing cells with different conductivity $\sigma$.

Let us consider the PDE on the above domain:

$$-\nabla \cdot (\sigma(x,y) \nabla T) = 0.$$

Let us also define a control surface centered at a vertex $(i,j)$. Integrating over the control volume we have:

$$\int_{\Omega_{i,j}} \nabla \cdot (\sigma \nabla V) dA = 0 \Rightarrow \oint_{\partial \Omega_{i,j}} \sigma \nabla V \cdot n ds,$$

where $\partial \Omega_{i,j}$ is the boundary aroung the volume (in our case a square) consisting of four faces (sides), namely East, West, North and South. For the east face (side) between vertices $(i,j)$ and $(i+1,j)$ the flux (current density) is defined as follows, using central differences:

$$J_x^{i+1/2,j}=-\sigma_{i+1/2,j} \frac{V_{i+1,j}-V_{i,j}}{\Delta x}$$

with

$$\sigma_{i+1/2,j}=\frac{2\sigma_{i,j}\sigma_{i+1,j}}{\sigma_{i,j}+\sigma_{i+1,j}}$$

computed throught the harmonic average. For the West face (side), between vertices $(i-1,j)$ and $(i,j)$, the flux is as follows:

$$J_x^{i-1/2,j}=-\sigma_{i-1/2,j}\frac{V_{i,j}-V_{i-1,j}}{\Delta x}$$

with 

$$\sigma_{i-1/2,j} = \frac{2\sigma_{i,j} \sigma_{i-1,j}}{\sigma_{i,j}+\sigma{i-1,j}}.$$

For the North face (side), between vertices $(i,j+1)$ and $(i,j)$, the flux is as follows:

$$J_y^{i,j+1/2}=-\sigma_{i,j+1/2} \frac{V_{i,j+1}-V_{i,j}}{\Delta y}$$

with 

$$\sigma_{i,j+1/2}=\frac{2\sigma_{i,j}\sigma_{i,j+1}}{\sigma_{i,j}+\sigma_{i,j+1}}.$$

For the South face (side), between vertices $(i,j)$ and $(i,j-1)$ teh flux is as follows:

$$J_y^{i,j-1/2}=-\sigma_{i,j-1/2} \frac{V_{i,j}-V_{i,j-1}}{\Delta y}$$

with

$$\sigma_{i,j-1/2}=\frac{2\sigma_{i,j}\sigma_{i,j-1}}{\sigma_{i,j}+\sigma_{i,j-1}}.$$

Using the fluxes over the faces we can enforce conservation at cell $(i,j)$ as follows:

$$\frac{J_x^{i+1/2,j}-J_x^{i-1/2,j}}{\Delta x}+\frac{J_y^{i,j+1/2}-J_y^{i,j-1/2}}{\Delta y}=0.$$

Substituting the fluxes corresponding to the faces to the conservation of flux above we obtain the five point conservative stencil:

$$-\frac{\sigma_{i+1/2,j}}{\Delta x^2} V_{i+1,j}-\frac{\sigma_{i,j+1/2}}{\Delta y^2} V_{i,j+1} + \left( \frac{\sigma_{i+1/2,j}+\sigma_{i-1/2,j}}{\Delta x^2} + \frac{\sigma_{i,j+1/2}+\sigma_{i,j-1/2}}{\Delta y^2} \right) V_{i,j} -\frac{\sigma_{i-1/2,j}}{\Delta x^2} V_{i-1,j}-\frac{\sigma_{i,j-1/2}}{\Delta y^2} V_{i,j-1} = 0.$$

All conductivities appearing in the stencil are given by the equations used in the definition of flux for each face (side) above. The derived scheme has some distinct advantages over standard FD:

1. Enforces conservation,
2. Handles discontinuous $\sigma$,
3. It retains physical fidelity,
4. It is more suitable for PDEs that are expressed in divergence form $(\nabla \cdot F = S)$, such as those arising in Electrostatics, Heat conduction, Fluid dynamics and Groundwater flow.

It should be noted that the quantities $\Delta x = \frac{L_x}{n_x}$ and $\Delta y=\frac{L_y}{n_y}$ are the mesh sizes corresponding to each axis.

After computing the potential distribution the temperature field should be computed. This is going to be achieved by considering the steady-state thermal equation.

One of the key differences of Flux conservative FD, as opposed to traditional FD, is that without boundary treatment the numerical scheme considers Neumann type boundary conditions, while traditional FD considers Dirichlet boundary conditions.

The physical interpretation of homogeneous Neumann boundary conditions, for the electrical potential PDE, correspond to current density $j_n$ through the boundary:

$$n\cdot \sigma (T) \nabla V = - j_n,$$

while in the case of electrically insulated surface, such as the one considered here due to the substrate, the boundary conditions are as follows:

$$n \cdot \nabla V = 0 \Rightarrow \frac{\partial V}{\partial n} = 0.$$

In our setup there are also the voltage pads, where the voltage is known a priori. In the left part (positive side) of the domain these points that correspond to the pad in the boundary are considered to have voltage equal to $1$ (Dirichlet). On the right side (negative side), the known voltage of the pad is considered to be $0$ (Dirichlet). In order to enforce these boundaries to the coefficient matrix the rows corresponding to these points should be substituted with the corresponding rows of the identity matrix $I$ of the same order as the coefficient matrix. This ensures homogeneous Neumann boundaries everywhere on the sides of the domain except the positive and negative pads.

### Steady-state thermal equation

The steady-state thermal equation is considered instead of the transient one:

$$\nabla \cdot (k(x,y) \cdot \nabla T) = - Q = - \sigma (x,y) |\nabla V|^2,$$

where $k(x,y)$ in $(W/(m \cdot K))$ is the [thermal conductivity](https://en.wikipedia.org/wiki/Thermal_conductivity_and_resistivity), which is different for every element inside the domain. The steady-state thermal equation is considered since:

* Joule heating has been applied long enough that temperature has stabilized,
* A DC or low-frequency problem, is analyzed, where thermal inertia dominates,
* If one cares only about the final temperature distribution and not how fast it heats up.

Because of the spatial variance of thermal conductivity, and more specifically the discontinuities due to the different materials. Thus, a similar approach for discretization  as the one followed for voltage (electrical potential) is going to be considered here also. Using the same grid and conservative FD discretization the thermal stencil has the following form:

$$\frac{k_{i+1/2,j}}{\Delta x^2} T_{i+1,j}+\frac{k_{i-1/2,j}}{\Delta x^2} T_{i-1,j}+\frac{k_{i,j+1/2}}{\Delta y^2} T_{i,j+1}+\frac{k_{i,j-1/2}}{\Delta y^2} T_{i,j-1}-\left( \frac{k_{i+1/2,j}+k_{i-1/2,j}}{\Delta x^2}+\frac{k_{i,j+1/2}+k_{i,j-1/2}}{\Delta y^2} \right) T_{i,j} = -Q_{i,j},$$

with:

$$Q_{i,j}=\sigma_{i,j} \left( \left(\frac{V_{i+1,j}-V_{i-1,j}}{2\Delta x} \right)^2 + \left( \frac{V_{i,j+1}-V_{i,j-1}}{2\Delta y} \right)^2 \right),$$

and:

$$k_{i+1/2,j} = \frac{2 k_{i,j} k_{i+1,j}}{k_{i,j} + k_{i+1,j}},$$

$$k_{i-1/2,j} = \frac{2 k_{i,j} k_{i-1,j}}{k_{i,j} + k_{i-1,j}},$$

$$k_{i,j+1/2} = \frac{2 k_{i,j} k_{i,j+1}}{k_{i,j} + k_{i,j+1}},$$

$$k_{i,j-1/2} = \frac{2 k_{i,j} k_{i,j-1}}{k_{i,j} + k_{i,j-1}}.$$

Again, harmonic mean is used to ensure conservation and continuity accross neighboring element with different thermal conductivity coefficients.

Heat dissipation is perfromed through the boundaries. This can be modelled by considering [Robin boundary conditions](https://en.wikipedia.org/wiki/Robin_boundary_condition) of the following form:

$$-k n\cdot \nabla T = h (T-T_{\infty})$$

where $T_{\infty}$ is the ambient temperature, which is set equal to the reference temperature $T_0=293 K$, and $h$ is the convective heat transfer coefficient $(W/(m^2 \cdot K))$. The value of the convective heat transfer coefficient is set to $5$ $W/(m^2 \cdot K))$ approximating natural convection of a vertical plate facing up. 

However, heat dissipation through the boundary is not enough to cool the IC. In practice heat dissipation also happens from the surface of the IC. Thus, a modification to the PDE is applied, namely:

$$\nabla \cdot (k(x,y) \cdot \nabla T) + h_{vol} (T-T_{\infty})= - Q = - \sigma (x,y) |\nabla V|^2,$$

where $h_{vol} = \frac{2h}{d}$ is the volumetric heat loss rate $(W/(m^3 \cdot K))$. The parameter $d$ is the thickness of the substrate, which was set to $3 mm$. This modification to the PDE leads to the addition of constant terms to the diagonal and the RHS of the sparse linear system.

The incorporation of Robin boundary conditions happens by substituting their discrete version into the corresponding coefficient matrix equations. The coefficient matrix however is already incorporating Neumann Boundary conditions, thus only the diagonal element and the RHS are affected. More specifically, the diagonal elements are increased by the quantity:

$$\frac{h \cdot h_x}{k_{i,j}}$$

or 

$$\frac{h \cdot h_y}{k_{i,j}}$$

depending on the boundary. While the RHS is increased by the quantities:

$$\frac{h \cdot h_x \cdot T_{\infty}}{k_{i,j}}$$

or 

$$\frac{h \cdot h_y \cdot T_{\infty}}{k_{i,j}}$$

depending on the boundary. 

An important note is that the Joule heating term $Q$ needs to be multiplied (scaled) by $h_x \cdot h_y$ to account for heat inside a cell. This returns the correct units required to acquire a physically meaningful result.

### Solution of the coupled system

The process of solving the coupled system involves the following steps:

1. Compute V
2. Compute T
3. Update $\sigma$ where required

This process is repeated until a steady solution is reached. The chosen stopping criterion for the iterative process was set to $|max(T^i)-max(T^{i+1}|<10^{-3}$. Some indicative results for the above setup are given below.

![4](https://github.com/user-attachments/assets/5ad1dca3-29c7-4df5-8663-94f07948b84e)
![3](https://github.com/user-attachments/assets/0c46c732-4b75-4a1d-8807-d0e937fa5cee)
![2](https://github.com/user-attachments/assets/7fb8c1d1-827e-4b9e-8874-6c1dd0902d81)
![1](https://github.com/user-attachments/assets/04d7494e-28cd-42b5-95f4-c7a785242a9f)


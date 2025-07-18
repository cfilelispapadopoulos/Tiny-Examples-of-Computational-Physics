# Electrothermal Coupling

## 2D Electrothermal Simulation of an Integrated Circuit Interconnect

The temperature distribution of an IC interconnect network due to Joule heating is a complex phenomenon that couples the electrical potential and temperature fields. Thus, the couple solution of two Partial Differential Equations simultaneously is required. In order to showcase this complex phenomenon, initially, an model integrated circuit will be defined.

### Model Integrated Circuit

The considered integrated circuit is given in the Figure below.

<img width="515" height="301" alt="icm" src="https://github.com/user-attachments/assets/63dfadea-41b1-4452-8873-0ca92e7a6590" />

The IC is composed of three resistors, where the first one is connected, in series, to the other two which are connected in parallel. All three stainless steel resistors with electrical conductivity $\sigma_R^0=1.45\times 10^{-6} (S/m)$ lie on a substrate ($0.2 m \times 0.2 m$) with electrical conductivity $\sigma_S = 10^{-6} (S/m)$. All components are connected with copper electrical wire (red) with electrical conductivity $\sigma_W=5.8\times 10^7 (S/m)$.

The conductivity of the wire as well as the substrate are considered to be independent of the temperature, while the conductivity of the resistors is considered to be varying with temperature:

$$\sigma_R(T) = \frac{\sigma_R^0}{1+\alpha (T-T_0)},$$

where $\alpha$ is the temperature coefficient (per K), $T$ is the current temperature of the resistor and $T_0$ is the reference temperature. It should be noted that the width of the wires and the length and width of resistors are significant and are considered parameters in the numerical modelling process that follows. In order to make the modelling process easier, the distance of the resistors and length of the cables is considered arbistrary, since it does not affect significantly the results of this example. The governing equations will be presented in the next part.

### Electrical Potential

Initially the electrical potential of the IC is required. The Electrical Potential can be computed using the following governing equation:

$$\nabla \cdot (\sigma(T) \nabla T) = 0.$$ 

This is a linear PDE which can be solved by either the Finite Differences or the Finite Elements method. Due to the normal geometry of the substrate (square) the Finite Differences method will be considered. In order to apply the FD method the substrate is discretized using a square grid composed of $n_x+1$ points in the $x$-direction and $n_y+1$ points in the $y$-direction. 

![grid](https://github.com/user-attachments/assets/6f3463de-5021-4a31-9f1e-80e5e06259a0)

It should be noted that the axis are switched because MATLAB utilizes column-wise storage. Thus, the column wise indexing corresponds to the $x$ direction. Moreover, the number of intervals per dimensions is denoted as $n_x$ for the $x$-axis and $n_y$ for the $y$-axis.

In order to solve the PDE the regular five point stencil is not going to yield a physically accurate solution. Instead, a conservative Flux-based discretization will be chosen, since the aforementioned PDE is derived from the charge conservation:

$$\nabla \cdot (\sigma(T) \nabla T) = 0 \Rightarrow \nabla \cdot J = 0,$$

Thus, the net current into each cell is exactly zero. Another important advantage of Flux conservative FD is that this type of discretization preserves continuity $\sigma_1 \frac{\partial V}{\partial n} = \sigma_2 \frac{\partial V}{\partial n}$ in case of material discontinuities, such as conductivity discontinuity, especially when combined with harmonic averaging:

$$\sigma_{harmonic} = \frac{2\sigma_1 \sigma_2}{\sigma_1 + \sigma_2},$$

in neighboing cells with different conductivity $\sigma$.

Let us consider the PDE on the above domain:

$$\nabla \cdot (\sigma(x,y) \nabla T) = 0.$$

Let us also define a control surface centered at a vertex $(i,j)$. Integrating over the control volume we have:

$$\int_{\Omega_{i,j}} \nabla \cdot (\sigma \nabla V) dA = 0 \Rightarrow \oint_{\partial \Omega_{i,j}} \sigma \nabla V \cdot n ds,$$

where $\partial \Omega_{i,j}$ is the boundary aroung the volume (in our case a square) consisting of four faces (sides), namely East, West, North and South. For the east face (side) between vertices $(i,j)$ and $(i+1,j)$ the flux is defined as follows, using central differences:

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

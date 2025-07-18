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

$$\nabla \cdot (\sigma(T) \nabla T) = 0$$ 

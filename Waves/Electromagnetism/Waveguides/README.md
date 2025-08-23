# TM modes on a rectilinear waveguide using FDTD (Finite Difference Time Domain)
## Description
Let us consider a 3D rectangular waveguide of the form given in the Figure below.

<img width="862" height="447" alt="image" src="https://github.com/user-attachments/assets/111f1c4f-de4f-4e63-b863-223110823da6" />

A Gaussian wideband pulse is transmitted across the x-axis, reflecting on the internal square obstacle at the center of the waveguide and propagating towards the $-x$ and $+x$ directions. The obstacle is made out Perfect Electric Conductor (PEC). The surrounding faces of the rectilinear waveguide are covered with Perfect Electric Conductor (PEC) which reflects incident waves. The two faces across the x-axis are considered to be absorbing boundaries. These absorbing boundaries are simulated using [Convolutional Perfectly Matched Layer (CPML)](https://onlinelibrary.wiley.com/doi/10.1002/1098-2760(20001205)27:5%3C334::AID-MOP14%3E3.0.CO;2-A) in order to minimize reflections

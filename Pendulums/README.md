# Different (not so simple) pendulum cases

The three example files correspond to the following:
1) Not so simple pendulum (simple_pendulum.m): Solving the ODE with Kick-Drift-Kick method without small angle approximation.
2) Not so simple not so aerodynamic pendulum (simple_pendulum_drag.m): Solving the ODE, but this time aerodynamic drag is considered, with Kick-Drift-Kick method without small angle approximation.
3) The double pendulum (double_pendulum.m): The system of 4 ODEs is solved using explicit RK4 (ode45) for non-stiff differential equations.

 All provided codes are modular and can be executed with arbitrary initial conditions. Also, they output graphs related to the simulation and variables of interest for the simulated phenomena.

## Not so simple pendulum
Let us consider a bob of mass $M$ suspended from a wall using a massless thread with length equal to $L$, as shown below:
![image](https://github.com/cfilelispapadopoulos/Tiny-Examples-of-Computational-Physics/assets/137081674/c6372223-a462-4823-bcca-ffe300393c40)

The forces acting on the bob is gravity $(mg)$ and tension from the thread. From Newton's second law we have the following:
$$F=ma=-mg sin\theta$$
or equivalently:
$$a=-gsin\theta \iff \frac{d^2 s}{dt^2}=-gsin\theta$$,
where $s$ denotes the arc length. But $s=L\theta$, so we can rewrite the above equation as follows:
$$\frac{d^2\theta}{dt^2}=-\frac{g}{L}sin\theta \iff \ddot{\theta}=-\frac{g}{L}sin\theta=f(t,\theta)$$.
If we avoid linearization using the small angle approximation $\theta\approx sin\theta$, then the above ODE can be solved using numerical methods, such as the [Leapfrom Kick-Drift-Kick](https://en.wikipedia.org/wiki/Leapfrog_integration) scheme, which is advantageous (in terms of stability) for oscillatory motion. Thus, we have:
$$\omega_{i+1/2}=\omega_i+f(t_i,\theta_i)\frac{\Delta t}{2},$$
$$\theta_{i+1}=\theta_i+\omega_{i+1/2}\Delta t,$$
$$\omega_{i+1}=\omega_{i+1/2}+f(t_{i+1},\theta_{i+1})\frac{\Delta t}{2},$$
where $a_i=\ddot{\theta}_i=f(t_i,\theta_i)$ is the angular acceleration during the $i$-th time step, $\omega_i=\dot{\theta}_i$ is the angular velocity during the $i$-th time step and $\Delta t$ is the time step. The solution at each time step is produced by evaluating the three last equations at each time step until the maximum prescribed time. The pendulum performs simple harmonic motion.

## Not so simple not so aerodynamic pendulum
This case is similar to the above with the only difference being the addition of aerodynamic [drag force](https://en.wikipedia.org/wiki/Drag_(physics)) $F_d=\frac{1}{2} \rho C_d A v^2$, where $\rho$ is the mass density of the fluid, $C_d$ is the drag coefficient (in the case of the ball $C_d=0.47$), $A$ is the cross sectional area (in the case of the ball $\pi r^2$) and $v$ is the flow speed (velocity) of the object relative to the fluid.

![image](https://github.com/cfilelispapadopoulos/Tiny-Examples-of-Computational-Physics/assets/137081674/f0c4a5f4-91ae-4d2d-b819-e87ed60f79b5)

Thus, analogously the differential equation has the following form:
$$\frac{d^2\theta}{dt^2}=-\frac{g}{L}sin\theta - sgn(\omega) \frac{L}{2M} \rho C_d A \omega^2 \iff \ddot{\theta}=f(t,\theta,\omega)$$,
where $sgn(x)$ denotes the [sign function](https://en.wikipedia.org/wiki/Sign_function) and $sgn(\omega)$ ensures that the drag force points always in the opposite direction of the tangential velocity $v=\omega L$. Using again the Leapfrom KDK scheme, as above, we can solve the ODE:
$$\omega_{i+1/2}=\omega_i+f(t_i,\theta_i,\omega_i)\frac{\Delta t}{2},$$
$$\theta_{i+1}=\theta_i+\omega_{i+1/2}\Delta t,$$
$$\omega_{i+1}=\omega_{i+1/2}+f(t_{i+1},\theta_{i+1},\omega_{i+1/2})\frac{\Delta t}{2}.$$
It shoule be noted that in this case the pendulum performes damped harmonic motion.

## Double Pendulum
The double pendulum is composed of a pendulum mounted on a wall and the second mounted onto the bob of the first. The two threads might be of different length ($L_1$ and $L_2$) and the bobs might also be of different mass ($M_1$ and $M_2$). An example of a double pendulum is given in the figure below.

![image](https://github.com/cfilelispapadopoulos/Tiny-Examples-of-Computational-Physics/assets/137081674/873fbe9d-3dd2-4f06-b42e-3cd4156a4984)

There are many ways to derive the ODEs governing the motion of a double pendulum including direct derivation through [kinematics](https://web.mit.edu/jorloff/www/chaosTalk/double-pendulum/double-pendulum-en.html) or through Calculus of Variations using either the [Lagrangian or Hamiltonian formulations](https://scienceworld.wolfram.com/physics/DoublePendulum.html). Because of the nature of this example (supposed to be small and easy to follow) we will follow the direct kinematic approach and leave the other two for the interested reader or more appropriate future examples.

TBC soon

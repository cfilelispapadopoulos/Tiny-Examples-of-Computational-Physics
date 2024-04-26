# Orbits

## Three - body problem and the Sun-Earth-Moon system
The [three body problem](https://en.wikipedia.org/wiki/Three-body_problem) is a famous problem in classical mechanics concerning the prediction of trajectories of three point masses knowing the initial positions and velocities. The bodies interact with each other only through mutual gravitational forces. This system exhibits chaotic behavior for most of the choices of initial conditions. Several researchers over the last 50 years, apart from the works of Euler and Lagrange, have proven different families of stable solutions such as the work of [Ivan Hristov et al](https://arxiv.org/abs/2308.16159) in 2023 or [Xiaoming Li and Shijun Liao](https://arxiv.org/abs/1705.00527) in 2017.

Here we will focus on a the stable system of the Sun - Earth - Moon which will be described using the [Lagrangian mechanics](https://en.wikipedia.org/wiki/Lagrangian_mechanics) approach in order to estimate the maximum distance between Sun and Earth and the trajectories of Earth and Moon. Initially, we have to define the Lagrangian of the system in Cartesian coordinates:
$$L(t,\vec{x},\vec{y},\vec{\dot{x}},\vec{\dot{y}}) = T(\vec{\dot{x}},\vec{\dot{y}})-V(\vec{x},\vec{y})$$
where $\vec{x},\vec{y}$ are vectors retaining the (time-dependent) coordinates, for all bodies. The Kinetic Energy of the system is of the form:
$$T = \frac{1}{2} m_1 (\dot{x}^2_1(t)+\dot{y}^2_1(t))+\frac{1}{2} m_2 (\dot{x}^2_2(t)+\dot{y}^2_2(t))+\frac{1}{2} m_3 (\dot{x}^2_3(t)+\dot{y}^2_3(t))=\frac{1}{2} \sum_{i=1}^3 \left[m_i (\dot{x}^2+\dot{y}^2) \right].$$
The potential energy of the three body system is of the form:
$$V(\vec{r})=-G\frac{m_1 m_2}{|r_1-r_2|}-G\frac{m_1 m_3}{|r_1-r_3|}-G\frac{m_2 m_3}{|r_2-r_3|}$$
or equivalently,
$$V(\vec{x},\vec{y})=-G\frac{m_1 m_2}{\sqrt{(x_1-x_2)^2+(y_1-y_2)^2}}-G\frac{m_1 m_3}{\sqrt{(x_1-x_3)^2+(y_1-y_3)^2}}-G\frac{m_2 m_3}{\sqrt{(x_2-x_3)^2+(y_2-y_3)^2}}.$$
The Euler - Lagrange equations for each dependent variable, i.e. $x_i (t),y_i (t)$, yield:
$$\frac{\partial L }{\partial x_i}-\frac{d}{dt} \frac{\partial L}{\partial \dot{x}_i}=0,i=1,2,3$$
$$\frac{\partial L }{\partial y_i}-\frac{d}{dt} \frac{\partial L}{\partial \dot{y}_i}=0,i=1,2,3$$
or equivalently:
$$\ddot{x}_1+G\frac{m_2}{|r_1-r_2|^3}(x_1-x_2)+G\frac{m_3}{|r_1-r_3|^3}(x_1-x_3)=0$$
$$\ddot{y}_1+G\frac{m_2}{|r_1-r_2|^3}(y_1-y_2)+G\frac{m_3}{|r_1-r_3|^3}(y_1-y_3)=0$$
$$\ddot{x}_2-G\frac{m_1}{|r_1-r_2|^3}(x_1-x_2)+G\frac{m_3}{|r_2-r_3|^3}(x_2-x_3)=0$$
$$\ddot{y}_2-G\frac{m_1}{|r_1-r_2|^3}(y_1-y_2)+G\frac{m_3}{|r_2-r_3|^3}(y_2-y_3)=0$$
$$\ddot{x}_3-G\frac{m_1}{|r_1-r_3|^3}(x_1-x_3)-G\frac{m_2}{|r_2-r_3|^3}(x_2-x_3)=0$$
$$\ddot{y}_3-G\frac{m_1}{|r_1-r_3|^3}(y_1-y_3)-G\frac{m_2}{|r_2-r_3|^3}(y_2-y_3)=0$$
In order to solve this system of six second order nonlinear ODEs we will convert them to a system of twelve first order nonlinear ODEs by substituting $v_i=\dot{x}_i,u_i=\dot{y}_i$:
$$\dot{x}_1=v_1$$
$$\dot{x}_2=v_2$$
$$\dot{x}_3=v_3$$
$$\dot{y}_1=u_1$$
$$\dot{y}_2=u_2$$
$$\dot{y}_3=u_3$$
$$\dot{v}_1=-G\frac{m_2}{|r_1-r_2|^3}(x_1-x_2)-G\frac{m_3}{|r_1-r_3|^3}(x_1-x_3)$$
$$\dot{v}_2= G\frac{m_1}{|r_1-r_2|^3}(x_1-x_2)-G\frac{m_3}{|r_2-r_3|^3}(x_2-x_3)$$
$$\dot{v}_3= G\frac{m_1}{|r_1-r_3|^3}(x_1-x_3)+G\frac{m_2}{|r_2-r_3|^3}(x_2-x_3)$$
$$\dot{u}_1=-G\frac{m_2}{|r_1-r_2|^3}(y_1-y_2)-G\frac{m_3}{|r_1-r_3|^3}(y_1-y_3)$$
$$\dot{u}_2= G\frac{m_1}{|r_1-r_2|^3}(y_1-y_2)-G\frac{m_3}{|r_2-r_3|^3}(y_2-y_3)$$
$$\dot{u}_3= G\frac{m_1}{|r_1-r_3|^3}(y_1-y_3)+G\frac{m_2}{|r_2-r_3|^3}(y_2-y_3)$$

We can solve this system using an integrator of the MATLAB environment such as the [ode45](https://uk.mathworks.com/help/matlab/ref/ode45.html). The masses of the Sun, Earth and Moon are, respectively, as follows:
$$m_1 = 1.989\times 10^{30}\ kg, m_2 = 5.972\times 10^{24}\ kg, m_3 = 7.34767309\times 10^{22}\ kg$$
The position of the Earth at $t=0$, relative to the Sun positioned at the origin $x_1(0)=0\ km,y_1(0)=0\ m$, is $x_2(0)=147098074\times 10^3\ m,y_2(0)=0\ m$. The position of the Moon relative to the Earth is at $t=0$ is $x_3(0) = 384400\times 10^3\ m, y_3 (0) = 0\ m$. These positions correspond to the case where the Earth is closest to the Sun. The corresponding velocities of the bodies at $t=0$ are $v_1 (0) = 0\ m/s, v_2 (0) = 0\ m/s, v_3(0) = 0\ m/s$ and $u_1(0)=0\ m/s, u_2(0)=30.3\times 10^3\ m/s, u_3(0) = u_2(0) + 1022\ m/s$. It should be noted that the relative velocity of the Moon with respect to earth is approximatelly $1022\ m/s$.

The results of the simulation are:

![orb](https://github.com/cfilelispapadopoulos/Tiny-Examples-of-Computational-Physics/assets/137081674/b35a65ce-c90f-4d01-b0a3-96b9003545db)
![dist1](https://github.com/cfilelispapadopoulos/Tiny-Examples-of-Computational-Physics/assets/137081674/a1d7ed66-d444-461d-8644-20a1e3cf089d)

The simulation leads to a maximum distance between Earth and Moon equal to $1.52551 \times 10^8\ km$, while the actual distance is $1.521\times 10^8\ km$, which leads to an estimation error of $0.3\\%$.

## Hofmann transfer orbit and the shooting method
The [Hofmann transfer obrit](https://en.wikipedia.org/wiki/Hohmann_transfer_orbit) is an orbital maneuver used to transfer a spacecraft or a satellite from one elliptic orbit of lower altitude, around a massive body, to another of higher altitude. The Hofmann transfer orbit is optimal in the sense of propelant required to accelerate the craft through two impulses, i.e. $\Delta v$, in order for it to reach the required velocity. In order to compute the two impulses $\Delta v_1$ and $\Delta v_2$ required for exiting the first orbit and entering the second, analytic formulas exist, such as the [vis-viva equation](https://en.wikipedia.org/wiki/Vis-viva_equation):
$$v^2 = GM \left(\frac{2}{r} - \frac{1}{a}\right)$$
where $G$ is the gravitational constant, $M$ is the mass of the body orbited by the craft, $r$ is the the distance between the centers of mass of the two bodies and $a$ is the length of the semi-major axis of the elliptic transfer trajectory.

Here we will follow a different approach, and a little more general, leveraging the Shooting method to solve the system of ODEs derived from the Lagrangian:
$$L = \frac{1}{2} m \left(\dot{r}^2+r^2\dot{\theta}^2 \right)+G \frac{mM}{r}$$
with various different initial velocities, in order to pinpoint the trajectory between to prescribed points, that do not necessarily lie on the semi-major axis. After, application of the two Euler-Lagrange equations we have:
$$\ddot{r}=r\dot{\theta}^2-G\frac{M}{r^2}$$
$m r^2 \dot{\theta} = c$
The second equation is the conservation of Angular momentum and by substitution to the first one we have:
$$\ddot{r}=\frac{L^2}{m^2 r^3}-G\frac{M}{r^2}$$
The solution of this equation can be performed as before with an integrator of the MATLAB environment such as the [ode45](https://uk.mathworks.com/help/matlab/ref/ode45.html). The initial and final point of the transfer are well known and the tangential velocities across these two trajectories can be computed at any point using Finite-Differences, e.g. at a point $r_i=r(t_i),\theta_i=\theta(t_i)$ of a trajectory:
$$v_i \approx \frac{\theta_{i+1}-\theta_{i-1}}{t_{i+1}-t_{i-1}} r_i.$$
The main idea is computing the trajectories for several (incrementally larger) initial velocities $v'_1 = v_1 + \Delta v_1$ and retain as an answer the one that is closer to the endpoint. The second impulse can be computed as $\Delta v_2 = v_2 - v'_2$, where $v_2$ is the velocity required at the prescribed point to retain the orbit, while $v'_2$ is the velocity that the craft attains after the transfer. Thus, technically two impulse burns are required in order to change orbit. For simplicity we do not consider any changes in mass caused by the impulses (which in reality require propellant). Moreover, we consider that the impulses (changes in velocity) are instant. It should be noted that the smaller the total $\Delta v = \Delta v_1 +\Delta v_2$ the more efficient the tranfer is. Thus, in the case of our approach if the beginning and ending of a transfer orbit lie on the semi - major axis then this is a Hohmann transfer orbit. We will denote the start and end points as percentages of the length of the total orbit. A Hofman orbit with total impulse $\Delta v=3572.61 m/s$ (marked with the red dashed line is as follows):

![hofmann](https://github.com/cfilelispapadopoulos/Tiny-Examples-of-Computational-Physics/assets/137081674/d6985821-670c-4410-a102-89a4895b91d7)

The circle in the center represents the Sun, while the mass of the craft is chosen to be $500 kg$. The other constants are described in the MATLAB script. The script is flexible allowing for changes in the orbits, in the velocities, etc. A more complicated transfer orbit, with $\Delta v = 3871.83 m/s$, is shown below:

![traj](https://github.com/cfilelispapadopoulos/Tiny-Examples-of-Computational-Physics/assets/137081674/1cf6de59-a62a-4ccb-ba79-fc263542369e)

It should be mentioned that the approach followed here for a shooting method is a trivial one and is based on brute force; more elaborate and powerful techniques exist in the literature.

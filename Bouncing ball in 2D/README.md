# Bouncing ball

In this example a bouncing ball of mass $M$ subject to the force of gravity is simulated. The ball is composed of a collection of $N$ points lying on a circle of user prescribed radius. Since the ball is not simulated as a point particle the code keeps track of the movement of the center of mass $cm$ as wel as every point of the ball. Upon impact the ball is deformed and the process of returning to its original shape is simulated as a series of massless damped springs stemming from the center of mass to the radius of the circle, with respect to Newton's 2nd law:

$$ m_p\frac{d^2 x}{dx^2 }= -bu - kx=-b \frac{dx}{dt}-kx$$

where $m_p$ is the mass of each point $m_p = \frac{m}{N}$, $k$ is the spring coefficient and $b$ is the damping parameter set to $\sqrt{4m_p k}$. The $x$ denotes the position of the spring attached mass with respect to the center of mass and $v$ the velocity relative to the center of mass. The springs are allowed to move only on the radial direction for simplicity. 

The equations of motion for the center of mass $cm$ and every point on the circle are integrated using the Leapfrom integrator in the Kick-Drift-Kick (KDK) form (second order accurate):

$$ v_{i+1/2} = v_{i}+a_i \frac{\Delta t}{2},$$

$$ x_{i+1} = x_{i}+v_{i+1/2} \Delta t,$$

$$ v_{i+1} = v_{i+1/2}+a_{i+1} \frac{\Delta t}{2}$$

The KDK form has the additive advantage that it is stable for variable time stepping apart from the stability for oscillatory motion. Variable timestepping is required to avoid breakdowns in cases of large accelerations and velocities.

## Initial Conditions and Constants
Initially the system starts at rest, but the code enables starting with initial velocity or acceleration. The maximum time (T_max) for the simulation and the time step (t_step) is also selected by the user and is measured in seconds. The value of the spring coefficient (spr_coeff) and the mass mass of the ball can be also selected by the user, however the damping coefficient (spr_damp) is computed based on the spring coefficient to ensure critical damping. The position of the ball (initial center of mass) (xcm = (xc,yc)) and the position of the center of the bounding box (xb,yb), as well as the radius of the ball (r) are given as inputs. All units follow the SI system.

## Updating velocity and position based on the acceleration
At each timestep the acceleration vector needs to be defined in order to compute the new position and velocities. All points of the ball and its center of mass are subject to gravity acceleration:

$$\vec{a}_p=\vec{a}_g=(0,-g)$$

$$\vec{a}_{cm}=\vec{a}_g=(0,-g)$$

where $g=9.81\frac{m}{s^2}$ is the acceleration of gravity. The points on the ball have another type of acceleration if for some reason they are displaced, either closer or away, more or less than the radius $r$, from the position of the center of mass (e.g. from a collision). This acceleration stems from Hooke's Law and the fact that we are treating them as points attached to massless damped springs. Thus, for a point $p$ that, at some point in time, is further away than rest position $r$, its acceleration is:

$$\vec{a}_p=\vec{a}_g-\frac{k}{M}(||\Delta\vec{x}_p||_2-r) \frac{\Delta\vec{x}_p}{||\Delta\vec{x}_p||_2}-\frac{b}{M}||\Delta\vec{v}_p||_2 \frac{\Delta\vec{x}_p}{||\Delta\vec{x}_p||_2}$$,

where

$$\Delta\vec{x_p}=\vec{x_p}-\vec{x}_{cm}$$

and

$$\Delta\vec{v_p}=\vec{v_p}-\vec{v}_{cm}$$.

The vectors and $\vec{x_p}$, $\vec{v_p}$ are the absolute position and velocity of a point (with respect to the box) while $\Delta\vec{x_p}$ and $\Delta\vec{v_p}$ are the relative position of a point with respect to the center of mass (where the spring is attached). It should be noted that the reaction force to the center of mass is not considered for simplicity. 

## Disclaimer
This piece of software is intended as an educational example showcasing the connections between theory and practice, avoiding intricate phenomena, and not as a real piece of research.

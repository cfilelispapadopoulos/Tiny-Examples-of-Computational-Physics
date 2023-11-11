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

$$\vec{a}_p=\vec{a}_g-\frac{k}{M}(||\Delta\vec{x}_p||_2-r) \frac{\Delta\vec{x}_p}{||\Delta\vec{x}_p||_2}-\frac{b}{M}||\Delta\vec{v}_p||_2 \frac{\Delta\vec{v}_p}{||\Delta\vec{v}_p||_2}$$,

where

$$\Delta\vec{x_p}=\vec{x_p}-\vec{x}_{cm}$$

and

$$\Delta\vec{v_p}=\vec{v_p}-\vec{v}_{cm}$$.

The vectors and $\vec{x_p}$, $\vec{v_p}$ are the absolute position and velocity of a point (with respect to the box) while $\Delta\vec{x_p}$ and $\Delta\vec{v_p}$ are the relative position of a point with respect to the center of mass (where the spring is attached). It should be noted that the reaction force to the center of mass is not considered for simplicity. The above formula is impractical and requires tracking positions, velocities and accelerations with respect to the box. The process of updating these can be simplified by separating the frames of reference.

## Two separate frames of reference
We can track the movement of the system by its center of mass (center of mass of the sphere) and then the particles on the reference frame of the center of mass. Thus, we will advance in time the center of mass, separately from the sourounding particles, and we will account for forces between them. That way we can simplify all equations and avoid errors stemming from either rounding or choice of timestep. The center of mass is under motion with variable acceleration, this is because of the gravity acceleration and the acceleration induced by the reaction force of the deformed points attached to the springs. For the center of mass acceleration can be computed through:

$$M\vec{a_{cm}}=-M\vec{g}-\vec{F_r}$$

where $\vec{F_r}$ is the net reaction force computed as the sum of the reaction forces of the deformed springs. The points around the center of mass are moving in an accelerated frame of reference. Their acceleration can be computed as:

$$m_p\vec{a_p}=-k(\vec{x}-\vec{x_0})-b\vec{v_p}$$

where $\vec{x_0}$ is the vector corresponding to the rest length of the springs with respect to the initial position of the center of mass which in its frame of reference lies always at (0,0). All quantities in the above equation are defined with respect to the reference frame corresponding to the initial position of the center of mass. In case some of the springs are compressed or elongated we should account for the net reaction force at the point where they are mounted (0,0). This can be computed as:

$$\vec{F_r}=\sum_p m_p \vec{a_p}$$.

Thus, thus upon collision the springs push back at the mounting point.

## Collisions with the Walls
All collisions, since the ball is deformable, are inellastic. The process of modelling collisions is straightforward, following the above simplifications. When the ball meets with the wall the points that touch the wall have velocity equal to zero and position equal to:

$$2 y_{max} - x_2 - 2 (x_{cm})_2$$,

$$2 y_{min} - x_2 - 2 (x_{cm})_2$$,

for the upper and lower walls while for the side walls, analogously:

$$2 x_{max} - x_1 - 2 (x_{cm})_1$$,

$$2 x_{min} - x_1 - 2 (x_{cm})_1$$.

The bodies are stuck to the wall until the reaction force of all the corresponding springs becomes greated than gravity, in order to accelerate the body to the opposite direction. Thus, collisions are affected significantly by the constants of the springs.

## Improvements and changes
The code can be further impoved by:
1) Taking into account the tension forces created by the movement of the pieces of material
2) Adding air resistance
3) Implement variable time stepping for large velocities and accelerations
4) Implement maximum compression - elongation limits
5) Simulate gas particle dynamics inside ball
and many more.

## Disclaimer
This piece of software is intended as an educational example showcasing the connections between theory and practice, avoiding intricate phenomena, and not as a real piece of research.

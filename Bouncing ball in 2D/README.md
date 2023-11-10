In this example a bouncing ball subject to the force of gravity is simulated. The ball is composed of a collection of $N$ points lying on a circle of user prescribed radius.
Since the ball is not simulated as a point particle the code keeps track of the movement of the center of mass $cm$ as wel as every point of the ball. Upon impact the ball is deformed
and the process of returning to its original shape is simulated as a series of massless damped springs stemming from the center of mass to the radius of the circle, with respect to Newton's 2nd law:

$$ m_p\frac{d^2 x}{dx^2 }= -bu - kx$$

where $m_p$ is the mass of each point $m_p = \frac{m}{N}, $k$ is the spring coefficient and $b$ is the damping parameter set to $\sqrt(4m_p k)$


This piece of software is intended as an educational example showcasing the connections between theory and practice, avoiding intricate phenomena, and not as a real piece of research.

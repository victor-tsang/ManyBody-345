# ManyBody-345
Classical gravitionational 3-body problem, also known as 
__*Pythagorean Three-Body Problem*__. 

In 1893, Meissel stated what is now called the Pythagorean three-body problem: 
Three masses in the ratio 3:4:5 are placed at rest at the vertices 
of a 3:4:5 right triangle. 
Burrau further investigated this problem in 1913. 
In 1967 Victor Szebehely and C. Frederick Peters established eventual 
escape for this problem using numerical integration, while at the same 
time finding a nearby periodic solution. 
https://en.wikipedia.org/wiki/Three-body_problem

Three masses in the ratio 3:4:5 are placed at rest at the vertices of
a 3:4:5 right triangle.

# Verlet integration

Equation of motion for classical conservative system:

`V(r(t))` = potential energy function

Force `F = - grad (V(r(t)) )`

Velocity Verlet:

Given initial conditions `r(0)`, `v(0)`, compute acceleration `a(0)`.

```text
for time step t

v(t + 0.5 dt) := v(t) + 0.5*dt*a(t)
r(t + dt) := r(t) + v(t + 0.5 dt)*dt

compute a(t+dt) based on computed r(t+dt).

v(t+dt) := v(t+0.5dt) + 0.5*dt*a(t+dt)
```

N.B. acceleration `a(t+dt)` can’t depend on `v` because `v(t+0.5dt)` is 
at different time point.


Ref: https://en.wikipedia.org/wiki/Verlet_integration

# Scales

The masses are in ratio 3:4:5 in mass and separation,
and located in relative positions 

| particle | position | mass | 
|----------|----------|------|
| M1       | (1,3)    | 3    |
| M2       | (-2,-1)  | 4    |
| M3       | ( 1,-1)  | 5    |

such that the centre of mass is (0,0).


# Task 1: Study of convergence 

Use several threads to run the integrations with different 
step size to see how they will converge when the step is 
"sufficiently small".

In Particle.cpp, test002(), which uses different dt to
approach a given time, then compare the total energy of the 
system, found that, the energy varies very much.

## Side-track to Halley Comet
Try to study a simple system, Halley comet, with the same Vector and
Particle classes as the 3-4-5 system.

When monitoring the delta of r, v, and a during a dt time step, if the
dt is too "coarse", the delta will be too large which may make the energy
changes "too much", i.e. not conserve. The total energies of the system, for 
the comet moves from the aphelion to the perihelion, are hugh different. 
However, when the comet comes back from the perihelion to aphelion, the
total energy restored to the value very similar to the starting value.

## The velocity Verlet method and variable time steps
It seems people studied the possibility of using variable time steps
such that, when the delta r, v are too large, it can decrease the time step to minimize
the error... 
[The velocity Verlet method and variable time steps](https://scicomp.stackexchange.com/questions/34491/the-velocity-verlet-method-and-variable-time-steps)
But the result seems not quite positive...

## Compare the total energy during the steps
Is the severe energy variance due to close "head-on" collision?

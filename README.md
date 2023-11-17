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


## Relevant works

### AY 212: Dynamical Astronomy Project 2
In [project two -- Burrau's problem of three bodies](https://www.ucolick.org/~laugh/oxide/projects/burrau.html)
some discussion is presented. From its links, we can eventually find a page of a
dynamical astronomy class [codes](https://www.ucolick.org/~laugh/oxide/codes/index.html) which
calculates the Pythagorean three-body problem. The code is
[fewbody.f](https://www.ucolick.org/~laugh/oxide/codes/fewbody.f.txt)
instead of [integrator.f](https://www.ucolick.org/~laugh/oxide/codes/integrator.f.txt) as originally stated in the web page.

The project employed [Bulirsch–Stoer algorithm](https://en.wikipedia.org/wiki/Bulirsch%E2%80%93Stoer_algorithm) to do the
number crunching. This algorithm is also available in boost C++ library ([Boost.Numeric.Odeint](https://www.boost.org/doc/libs/1_83_0/libs/numeric/odeint/doc/html/index.html)):
[boost/numeric/odeint/stepper/bulirsch_stoer.hpp](https://www.boost.org/doc/libs/1_83_0/boost/numeric/odeint/stepper/bulirsch_stoer.hpp)

See also [What the difference is between Størmer Verlet and regular Verlet method?](https://physics.stackexchange.com/questions/447056/what-the-difference-is-between-st%C3%B8rmer-verlet-and-regular-verlet-method)


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

# Runge-Kutta with Adaptive step size
Boost [odeint library](https://www.boost.org/doc/libs/1_83_0/libs/numeric/odeint/doc/html/index.html)
provides handy ODE initial value problem solvers. We can
use its adaptive step size algorithms to better control the error.

See the test003(), test004() and test005() in HalleyComet.cpp. 
Also, the odeint supports Boost multiprecision, which can be a very
powerful tools.

During the hands on tests, odeint tries to identify the user-defined types are
containers or not. To tell odeint the user-defined type, such as boost::multiprecision::number<> is not 
a container type, defines a class like the following:

```c++

// See also https://github.com/victor-tsang/boost_odeint_hands_on/blob/main/README.md

#include<boost/numeric/odeint.hpp>
#include<boost/multiprecision/cpp_dec_float.hpp>
#include<boost/multiprecision/mpfr.hpp>
//#include<boost/multiprecision/gmp.hpp>

// https://stackoverflow.com/questions/58324974/why-does-odeint-fail-with-the-newer-versions-of-odeint
//
// The problem seems to be that boost::multiprecision::number,
// what makes mpf_float_100 (and every other Boost.Multiprecision type) work,
// has an associated value_type since Boost 1.68 and because of that Boost.Numeric.Odeint
// treats it as a container when it is not.
// The way that Odeint checks whether a type is a container is by
// using a trait: has_value_type, specializing that trait for number should work:
template< typename Backend, boost::multiprecision::expression_template_option Option >
struct has_value_type<boost::multiprecision::number<Backend,Option> >:boost::mpl::false_{};
```
Make the class `has_value_type< YourType >::value` gives `false`.

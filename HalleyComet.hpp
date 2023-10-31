//
// Created by Victor Tsang on 2023/10/31.
//

#ifndef MANYBODY_345_HALLEYCOMET_HPP
#define MANYBODY_345_HALLEYCOMET_HPP

#include"Vector.hpp"
#include"Particle.hpp"

#include<numbers>

namespace halley_comet
{
// Taken from victor@pingu.local:/opt/Programming/NumberCrunching/astronomy/celestial001.cpp

// In a simplified model of our solar system we simulate the motion
// of the Sun, the Earth, the planet Jupiter, and a comet similar
// to Halley’s Comet. We restrict the orbits to a two-dimensional
// plane and put the sun in the plane’s origin. Each astronomical
// body is represented by one particle. Between any two particles
// acts the gravitational potential (2.42). A set of initial
// values is given below:
//
// item     mass     r0         v0
// -------  -------  ---------  ----------
// Sun      1        (0,0)      (0,0)
// Earth    3.0e-6   (0,1)      (-1,0)
// Jupiter  9.55e-4  (0,5.36)   (-0.425,0)
// Halley	  1e-14    (34.75,0)  (0,0.0296)
//
//  delta t = 0.015, total t=468.5
//
//
// To build:
//   clang++ -std=c++14 -Wall -O2 -o bin/celestial001_clang celestial001.cpp
//   g++-10 -std=c++14 -Wall -O2 -o bin/celestial001_gcc celestial001.cpp
//

//  plot "output/HalleyComet3.data" using 3:4 w points lc 1 title 'delta T = 1 hour',  "output/HalleyComet.data" using 3:4 with points lc 2 title 'delta T = 1 day' , "output/HalleyComet2.data" u 3:4 w points lc 3 title 'delta T = 1 week'
//  plot from dense to less dense...
//


  struct system
  {
    Particle sun, earth,jupiter,halley;
  };


  void test();

}; // namespace halley_comet


#endif //MANYBODY_345_HALLEYCOMET_HPP

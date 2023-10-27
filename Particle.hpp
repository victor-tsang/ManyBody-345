//
// Created by Victor Tsang on 2023/10/27.
//

#ifndef MANYBODY_345_PARTICLE_HPP
#define MANYBODY_345_PARTICLE_HPP

#include"Vector.hpp"

template<typename T>
struct basic_particle
{
  using value_type=T;
  using vector_type=basic_vector2<T>;

  vector_type r;
  vector_type v;
  vector_type a;
  value_type mass;
};

using Particle=basic_particle<double>;

#endif //MANYBODY_345_PARTICLE_HPP

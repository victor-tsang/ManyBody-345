//
// Created by Victor Tsang on 2023/11/5.
//

#ifndef MANYBODY_345_COLLISIONDETECTION_HPP
#define MANYBODY_345_COLLISIONDETECTION_HPP

#include"Particle.hpp"

namespace CollisionDetection
{

  class Collision
  {
    double dt_;

    bool possible_collide_1st_degree_(const Particle &a, const Particle &b) const noexcept;

  public:
    Collision(double dt);
    ~Collision() = default;

    int collide(const Particle &a, const Particle &b);
  };

}; // namespace CollisionDetection


#endif //MANYBODY_345_COLLISIONDETECTION_HPP

//
// Created by Victor Tsang on 2023/11/5.
//

#include"CollisionDetection.hpp"

CollisionDetection::Collision::Collision(double dt)
:dt_(dt)
{
}

int CollisionDetection::Collision::collide(const Particle &a, const Particle &b)
{
  if(possible_collide_1st_degree_(a, b))
  {
    return 1;
  }
  return 0;
}

bool CollisionDetection::Collision::possible_collide_1st_degree_(const Particle &a, const Particle &b) const noexcept
{
  //                                                   ->                ->
  // Take each particle as a circle/disc with radius | v |*dt, centre at r.
  // If any 2 discs overlap, they may collide.
  // two discs overlap when the distance between 2 centres are less than the sum
  // of their radii: |r_ab| < |v_a|*dt + |v_b|*dt
  return (a.r-b.r).norm() < dt_ * (a.v.norm()+b.v.norm());
}

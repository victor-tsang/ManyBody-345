//
// Created by Victor Tsang on 2023/10/27.
//

#include"Particle.hpp"
#include<iostream>
#include<stdexcept>
#include<cassert>

namespace {

  const Particle M1{{ 1, 3},{0,0},{0,0},3};
  const Particle M2{{-2,-1},{0,0},{0,0},4};
  const Particle M3{{ 1,-1},{0,0},{0,0},5};

  class Runner
  {
  protected:
    Particle m1_, m2_, m3_;
    Particle::value_type G_;
    Particle::value_type dt_;
    Particle::value_type half_dt_;

    void update_acceleration_()
    {

    }

  public:
    Runner(Particle::value_type dt, Particle::value_type G=1):
    m1_(M1), m2_(M2), m3_(M3), G_(G), dt_(dt), half_dt_(0.5*dt)
    {
      if(dt<0)
      {
        throw std::invalid_argument("dt is negative.");
      }
      assert(m1_.mass==3);
      assert(m2_.mass==4);
      assert(m3_.mass==5);

      assert(m1_.r.x==1);
      assert(m1_.r.y==3);
      assert(m2_.r.x==-2);
      assert(m2_.r.y==-1);
      assert(m3_.r.x==1);
      assert(m3_.r.y==-1);

      assert(m1_.v.x==0);
      assert(m1_.v.y==0);
      assert(m2_.v.x==0);
      assert(m2_.v.y==0);
      assert(m3_.v.x==0);
      assert(m3_.v.y==0);

      update_acceleration_();
    }

    virtual ~Runner()=default;

  };

  void test001()
  {
    Particle a{{0,0},{0,0},{0,0},0};
  }
};

void test_particle()
{
  test001();
}

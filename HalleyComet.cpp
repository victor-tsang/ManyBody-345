//
// Created by Victor Tsang on 2023/10/31.
//

#include"HalleyComet.hpp"
#include"Vector.hpp"
#include"Particle.hpp"

#include<iostream>
#include<functional>

namespace halley_comet
{
// item     mass     r0         v0
// -------  -------  ---------  ----------
// Sun      1        (0,0)      (0,0)
// Earth    3.0e-6   (0,1)      (-1,0)
// Jupiter  9.55e-4  (0,5.36)   (-0.425,0)
// Halley	  1e-14    (34.75,0)  (0,0.0296)

  const Particle Sun{{0,0},{0,0},{0,0},1.};
  const Particle Earth{{0,1},{-1,0},{0,0},3e-6};
  const Particle Jupiter{{0,5.36},{-0.425,0},{0,0},9.55e-4};
  const Particle Halley{{34.75,0},{0,0.0296},{0,0},1e-14};


  namespace two_body
  {
    // values are in metric units.
    const Particle Sun{{0,0},{0,0},{0,0},2e30}; // assume fixed point mass
    //const Particle Halley{{5.28e12,0},{0,9.13e2},{0,0},1e-14/2e30}; // at aphelion 遠日點 (N.B: 近日點： perihelion)
    const Particle Halley{{5.28e12,0},{0,9.13e2},{0,0},2.2e14}; // at aphelion 遠日點 (N.B: 近日點： perihelion)
    const double G{6.67428e-11};
    const double GMsun{G*Sun.mass};

  }; // namespace halley_comet::two_body

}; // namespace halley_comet

namespace {


  class two_body_halley_runner
  {
    std::function<void(double,const Vector2 &,const Particle &)> change_direction_callback_{};
    Particle halley_;
    double t_;
    double dt_;
    double half_dt_;
    int x_movement_direction_; // +ve = RHS, increasing x, -ve = LHS, decreasing x.

    void update_acceleration_()
    {
      halley_.a = halley_.r*(-halley_comet::two_body::GMsun*halley_.r.inverse_norm_cubic());
    }

  public:
    // void(double time,const Vector2 &previous_position,const Particle &halley)
    using change_direction_observer=std::function<void(double,const Vector2 &,const Particle &)>;

    two_body_halley_runner(double dt)
    : halley_(halley_comet::two_body::Halley),t_(0),dt_(dt),half_dt_(dt/2.),x_movement_direction_(-1)
    {
      update_acceleration_();
    }

    ~two_body_halley_runner() = default;

    void set(const change_direction_observer &r)
    {
      change_direction_callback_=r;
    }

    void set(change_direction_observer &&r)
    {
      change_direction_callback_=std::move(r);
    }

    void step()
    {
      halley_.v += halley_.a*half_dt_;
      auto previous_r = halley_.r;
      halley_.r += halley_.v * dt_;
      update_acceleration_();
      halley_.v += halley_.a*half_dt_;
      t_ += dt_;

      if(halley_.r.x > previous_r.x)
      {
        if(x_movement_direction_<=0)
        {
          if(change_direction_callback_)
          {
            change_direction_callback_(t_,previous_r,halley_);
          }
        }
        x_movement_direction_=1;
      }
      else if(halley_.r.x < previous_r.x)
      {
        if(x_movement_direction_>=0)
        {
          if(change_direction_callback_)
          {
            change_direction_callback_(t_,previous_r,halley_);
          }
        }
        x_movement_direction_ = -1;
      }
      else
      {
        x_movement_direction_=0;
      }
    }

    double time() const noexcept
    {
      return t_;
    }

    double potential_energy() const noexcept
    {
      return -halley_comet::two_body::GMsun*halley_.mass/halley_.r.norm();
    }

    double kinetic_energy() const noexcept
    {
      return halley_.ke();
    }

    double total_energy() const noexcept
    {
      return potential_energy() + kinetic_energy();
    }
  };

  void test001()
  {
    std::cout<<"halley_comet::(anonymous)::test001: "<<std::endl;

    std::cout<<"Test two-body Halley comet celestial motion:"<<std::endl;

    two_body_halley_runner::change_direction_observer observer=[](double time, const Vector2 previous_position, Particle const &halley){
      std::cout<<"*** at time "<<time<<" ("<<(time/(365.25*24*3600))<<" years)"<<": "<<previous_position<<" -> "<<halley.r<<", v = "<<halley.v<<" ("<<halley.v.norm()<<")."<<std::endl;

      // PE = - GMm/r
      double pe = -halley_comet::two_body::GMsun*halley.mass/halley.r.norm();
      double ke = halley.ke();
      std::cout<<"    PE = "<<pe<<", KE = "<<ke<<", total energy = "<<(pe+ke)<<std::endl;
    };

    double target_time=90*365.25*24*3600;
    for(double dt:{3600, 1800, 60})
    {
      two_body_halley_runner runner(dt);
      runner.set(observer);
      //runner.set([](double t,auto previous_position, auto halley){});
      std::cout<<"looping dt = "<<dt<<"s for "<<target_time<<"s..."<<std::endl;
      std::cout<<"    PE = "<<runner.potential_energy()<<", KE = "<<runner.kinetic_energy()<<", total energy = "<<runner.total_energy()<<std::endl;
      do
      {
        runner.step();
      }
      while(runner.time()<target_time);
      std::cout<<"    at time "<<runner.time()<<", PE = "<<runner.potential_energy()<<", KE = "<<runner.kinetic_energy()<<", total energy = "<<runner.total_energy()<<std::endl;
    }

    std::cout<<"halley_comet::(anonymous)::test001: done."<<std::endl;
  }

  // molecular_dynamics_2013.pdf
  // Chapter 10, Section 10.1
  void test002()
  {
    std::cout<<"halley_comet::(anonymous)::test002: "<<std::endl;

    // x(t + dt) = x(t) + v(t) dt + 1/2 a(t) dt^2
    // v(t + dt) = v(t) + [ a(t) + a(t + dt)] dt/2
    // a = GMsun / r^2

    double target_time=90*365.25*24*3600;

    for(double dt:{3600, 1800, 60})
    {
      std::cout<<"looping dt = "<<dt<<"s for "<<target_time<<"s..."<<std::endl;

      Vector2 r=halley_comet::two_body::Halley.r;
      Vector2 v=halley_comet::two_body::Halley.v;
      Vector2 a=halley_comet::two_body::Halley.a;
      Vector2 previous_a=a;
      Vector2 previous_r=r;
      double t=0;
      double pe=-halley_comet::two_body::GMsun*halley_comet::two_body::Halley.mass/r.norm();
      double ke=0.5*halley_comet::two_body::Halley.mass*v.norm_square();
      int x_moving_direction{-1};
      int y_moving_direction{1};
      double half_dt_square=0.5*dt*dt;
      double half_dt=dt/2;

      auto update_energy=[&r, &v, &a,&pe,&ke]{
        pe=-halley_comet::two_body::GMsun*halley_comet::two_body::Halley.mass/r.norm();
        ke=0.5*halley_comet::two_body::Halley.mass*v.norm_square();
      };

      auto dump_energy=[&r, &v, &a,&t, &previous_r,&pe,&ke,&update_energy]{
        update_energy();

        std::cout<<"*** at time "<<t<<" ("<<(t/(365.25*24*3600))<<" years)"<<": "<<previous_r<<" -> "<<r<<", v = "<<v<<" ("<<v.norm()<<")."<<std::endl;
        std::cout<<"    PE = "<<pe<<", KE = "<<ke<<", total energy = "<<(pe+ke)<<std::endl;
      };

      std::cout<<"    PE = "<<pe<<", KE = "<<ke<<", total energy = "<<(ke+pe)<<std::endl;
      while(t<target_time)
      {
        // x(t + dt) = x(t) + v(t) dt + 1/2 a(t) dt^2
        previous_r=r;
        r += v * dt + a * half_dt_square;

        // a = GMsun / r^2
        previous_a=a;
        a = r*(-halley_comet::two_body::GMsun*r.inverse_norm_cubic());

        // v(t + dt) = v(t) + [ a(t) + a(t + dt)] dt/2
        v += (a+previous_a)*half_dt;
        t+=dt;

        if(r.x>=previous_r.x) // move to right
        {
          if(x_moving_direction<=0)
          {
            dump_energy();
          }
          x_moving_direction=1;
        }
        else if(r.x<previous_r.x) // move to left
        {
          if(x_moving_direction>=0)
          {
            dump_energy();
          }
          x_moving_direction=-1;
        }
      }
      pe=-halley_comet::two_body::GMsun*halley_comet::two_body::Halley.mass/r.norm();
      ke=0.5*halley_comet::two_body::Halley.mass*v.norm_square();
      std::cout<<"    at time "<<t<<", PE = "<<pe<<", KE = "<<ke<<", total energy = "<<(ke+pe)<<std::endl;
    }

    //std::cout<<"Halley mass = "<<halley_comet::two_body::Halley.mass<<" (Wikipedia: 2.2 x 10^14 kg)"<<std::endl;

    std::cout<<"halley_comet::(anonymous)::test002: done."<<std::endl;
  }

  void test003()
  {
    std::cout<<"halley_comet::(anonymous)::test003: "<<std::endl;
    std::cout<<"halley_comet::(anonymous)::test003: done."<<std::endl;
  }

  void test004()
  {
    std::cout<<"halley_comet::(anonymous)::test004: "<<std::endl;
    std::cout<<"halley_comet::(anonymous)::test004: done."<<std::endl;
 }

  void test005()
  {
    std::cout<<"halley_comet::(anonymous)::test005: "<<std::endl;
    std::cout<<"halley_comet::(anonymous)::test005: done."<<std::endl;
  }
}; // anonymous namespace

void halley_comet::test()
{
  test001();
  test002();
}

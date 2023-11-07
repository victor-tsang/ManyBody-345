//
// Created by Victor Tsang on 2023/10/27.
//

#include"Particle.hpp"
#include"CollisionDetection.hpp"

#include<iostream>
#include<stdexcept>
#include<cassert>

namespace {

  const int Trace = 1;

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
    Particle::value_type t_;
    Particle::value_type inverse_total_mass_;
    Particle::value_type total_pe_;
    double variation_{0.05}; // +/- percentage.
    size_t tolorance_break_{0};
    double initial_total_energy_{0};
    CollisionDetection::Collision collision_detector_;

    void update_acceleration_()
    {
      // calculate forces
      //
      // force exerted in m1: F1 = F12 + F13
      // force exerted in m2: F2 = F21 + F23 = F23 - F12
      // force exerted in m3: F3 = F31 + F32 = -(F13 + F23)
      //
      // force from i to j:
      //   Fij = G mi mj rij / (rij)^3
      //

      Vector2 r12{m2_.r - m1_.r};
      Vector2 r13{m3_.r - m1_.r};
      Vector2 r23{m3_.r - m2_.r};

      if(Trace && 0)
      {
        std::cout<<"Runner::update_acceleration_: r1="<<m1_.r<<", r2="<<m2_.r<<", r13="<<r12<<", norm="<<r12.norm()<<std::endl;
        std::cout<<"Runner::update_acceleration_: r1="<<m1_.r<<", r3="<<m3_.r<<", r13="<<r13<<", norm="<<r13.norm()<<std::endl;
        std::cout<<"Runner::update_acceleration_: r2="<<m2_.r<<", r3="<<m3_.r<<", r23="<<r23<<", norm="<<r23.norm()<<std::endl;
      }

      r12 *= r12.inverse_norm_cubic()*G_;
      r13 *= r13.inverse_norm_cubic()*G_;
      r23 *= r23.inverse_norm_cubic()*G_;

      m1_.a = r12*m2_.mass + r13*m3_.mass;
      m2_.a = r23*m3_.mass - r12*m1_.mass;
      m3_.a = -(r13*m1_.mass + r23*m2_.mass);
    }

    void update_total_pe_()
    {
      Vector2 r12{m2_.r - m1_.r};
      Vector2 r13{m3_.r - m1_.r};
      Vector2 r23{m3_.r - m2_.r};

      total_pe_=-G_*(m1_.mass*m2_.mass/r12.norm() + m1_.mass*m3_.mass/r13.norm() + m2_.mass*m3_.mass/r23.norm());
    }

    int compare_(Particle const &prev, Particle const &current) noexcept
    {
      int ret;

      ret=compare_(prev.r, current.r,"r");
      ret+=compare_(prev.v, current.v,"v");
      ret+=compare_(prev.a, current.a,"a");

      if(ret)
        ++tolorance_break_;

      return ret;
    }

    int compare_(Vector2 const &prev, Vector2 const &current, const char *msg) const noexcept
    {
      auto diff=current - prev;
      if((prev.x && abs(diff.x/prev.x)>variation_) || (prev.y && abs(diff.y/prev.y)>variation_))
      {
        std::cout<<"****** at t = "<<t_<<", "<<msg<<": "<<prev<<" -> "<<current<<" break the tolorance."<<std::endl;
        return 1;
      }
      return 0;
    }

  public:
    using value_type = Particle::value_type;

    Runner(Particle::value_type dt, Particle::value_type G=1):
    m1_(M1), m2_(M2), m3_(M3),
    G_(G), dt_(dt), half_dt_(0.5*dt),t_(0),
    inverse_total_mass_(1.0/(m1_.mass+m2_.mass+m3_.mass)),
    collision_detector_(dt)
    {
      if(dt<0)
      {
        throw std::invalid_argument("dt is negative.");
      }
      assert(m1_.mass==3);
      assert(m2_.mass==4);
      assert(m3_.mass==5);

      assert(inverse_total_mass_ == 1./(m1_.mass+m2_.mass+m3_.mass));
      assert(inverse_total_mass_ == 1./12.);

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

      //initial_total_energy_=-G_*769./60.;
      initial_total_energy_=total_energy();
    }

    virtual ~Runner()=default;

    void step()
    {
      if(tolorance_break_ > 10)
      {
        t_ += dt_;
        return;
      }

      auto const previous_m1=m1_;
      auto const previous_m2=m2_;
      auto const previous_m3=m3_;

      m1_.v += m1_.a*half_dt_;
      m2_.v += m2_.a*half_dt_;
      m3_.v += m3_.a*half_dt_;

      m1_.r += m1_.v*dt_;
      m2_.r += m2_.v*dt_;
      m3_.r += m3_.v*dt_;

      update_acceleration_();

      m1_.v += m1_.a*half_dt_;
      m2_.v += m2_.a*half_dt_;
      m3_.v += m3_.a*half_dt_;

      t_ += dt_;
      compare_(previous_m1,m1_);
      compare_(previous_m2,m2_);
      compare_(previous_m3,m3_);
    }

    void set_variation(double variation)
    {
      variation_ = variation;
    }

    int step_with_energy_check()
    {
      if(tolorance_break_ > 10)
      {
        t_ += dt_;
        return 1;
      }

      m1_.v += m1_.a*half_dt_;
      m2_.v += m2_.a*half_dt_;
      m3_.v += m3_.a*half_dt_;

      m1_.r += m1_.v*dt_;
      m2_.r += m2_.v*dt_;
      m3_.r += m3_.v*dt_;

      update_acceleration_();

      m1_.v += m1_.a*half_dt_;
      m2_.v += m2_.a*half_dt_;
      m3_.v += m3_.a*half_dt_;

      t_ += dt_;

      update_total_pe_();
      auto energy = total_pe_ + m1_.ke() + m2_.ke() + m3_.ke();
      if(variation_ < (abs((energy-initial_total_energy_)/initial_total_energy_)))
      {
        ++tolorance_break_;

        std::cout<<"\t!!!!!!!! tolorance violation: at time "<<t_<<", energy = "<<energy<<", initial total energy = "<<initial_total_energy_<<std::endl;
        Vector2 r12{m2_.r - m1_.r};
        Vector2 r13{m3_.r - m1_.r};
        Vector2 r23{m3_.r - m2_.r};
        std::cout<<"\t!!!!!!!!  r_12 = "<<r12<<", r_13 = "<<r13<<", r_23 = "<<r23<<std::endl;

        return -1;
      }

      return 0;
    }

    int step_with_collision_check()
    {
      if(tolorance_break_ > 10)
      {
        t_ += dt_;
        return 1;
      }

      m1_.v += m1_.a*half_dt_;
      m2_.v += m2_.a*half_dt_;
      m3_.v += m3_.a*half_dt_;

      m1_.r += m1_.v*dt_;
      m2_.r += m2_.v*dt_;
      m3_.r += m3_.v*dt_;

      update_acceleration_();

      m1_.v += m1_.a*half_dt_;
      m2_.v += m2_.a*half_dt_;
      m3_.v += m3_.a*half_dt_;

      t_ += dt_;

      if(collision_detector_.collide(m1_,m2_)||collision_detector_.collide(m1_,m3_)||collision_detector_.collide(m2_,m3_))
      {
        ++tolorance_break_;

        std::cout<<"\t!!!!!!!! possible collision: at time "<<t_<<std::endl;
        Vector2 r12{m2_.r - m1_.r};
        Vector2 r13{m3_.r - m1_.r};
        Vector2 r23{m3_.r - m2_.r};
        std::cout<<"\t!!!!!!!!  r_12 = "<<r12<<", r_13 = "<<r13<<", r_23 = "<<r23<<std::endl;

        return -1;
      }

      return 0;
    }


    value_type time() const noexcept
    {
      return t_;
    }

    Vector2 centre_of_mass() const noexcept
    {
      Vector2 ret{0,0};
      ret.x = (m1_.r.x*m1_.mass + m2_.r.x*m2_.mass + m3_.r.x*m3_.mass) * inverse_total_mass_;
      ret.y = (m1_.r.y*m1_.mass + m2_.r.y*m2_.mass + m3_.r.y*m3_.mass) * inverse_total_mass_;

      return ret;
    }

    value_type total_pe() noexcept
    {
      update_total_pe_();
      return total_pe_;
    }

    value_type total_energy() noexcept
    {
      update_total_pe_();
      return total_pe_ + m1_.ke() + m2_.ke() + m3_.ke();
    }

    void dump_masses(std::ostream &os) const
    {
      os<<"m1: mass = "<<m1_.mass<<", at "<<m1_.r<<", v = "<<m1_.v<<", speed = "<<m1_.v.norm()<<", KE = "<<m1_.ke()<<std::endl;
      os<<"m2: mass = "<<m2_.mass<<", at "<<m2_.r<<", v = "<<m2_.v<<", speed = "<<m2_.v.norm()<<", KE = "<<m2_.ke()<<std::endl;
      os<<"m3: mass = "<<m3_.mass<<", at "<<m3_.r<<", v = "<<m3_.v<<", speed = "<<m3_.v.norm()<<", KE = "<<m3_.ke()<<std::endl;
    }
  }; // class Runner

  void test001()
  {
    Particle a{{0,0},{0,0},{0,0},0};

    //Runner runner(1e-4,1.0);
    Runner runner(1e-4,1.0e-10);

    const size_t loops=1000000;

    std::cout<<"test001:\n\ttime is "<< runner.time()<<std::endl;
    std::cout<<"\tcentre of mass = "<< runner.centre_of_mass()<<", total PE = "<<runner.total_pe()<<", total energy = "<<runner.total_energy()<<std::endl;

    std::cout<<"Running "<<loops<<" iterations"<<std::endl;

    for(size_t i=0;i<loops;++i)
      runner.step();

    std::cout<<"test001:\n\ttime is "<< runner.time()<<std::endl;
    std::cout<<"\tcentre of mass = "<< runner.centre_of_mass()<<", total PE = "<<runner.total_pe()<<", total energy = "<<runner.total_energy()<<std::endl;

    std::cout<<"\n\ntest001: done."<<std::endl;
  }

  void test002()
  {
    const double target_time=100;
    //const double G=1.0;
    const double G=1.0e-10;

    std::cout<<"test002:"<<std::endl;
    size_t loops;

    for(double dt:{1e-6, 1e-7})
    //for(double dt:{1e-4, 1e-5, 1e-6, 1e-7, 1e-8})
    {
      Runner runner(dt,G);

      std::cout<<"dt = "<<dt<<", estimated "<<(target_time/dt)<<" steps."<<std::endl;
      std::cout<<"\ttime is "<< runner.time()<<std::endl;
      std::cout<<"\tcentre of mass = "<< runner.centre_of_mass()<<", total PE = "<<runner.total_pe()<<", total energy = "<<runner.total_energy()<<std::endl;
      runner.dump_masses(std::cout);

      loops=0;
      while(runner.time()<target_time)
      {
        runner.step();
        ++loops;
      }

      std::cout<<"\ttime is "<< runner.time()<<", "<<loops<<" steps."<<std::endl;
      std::cout<<"\tcentre of mass = "<< runner.centre_of_mass()<<", total PE = "<<runner.total_pe()<<", total energy = "<<runner.total_energy()<<std::endl;
      runner.dump_masses(std::cout);
      std::cout<<std::endl;
      std::cout<<std::endl;
    }

    std::cout<<"test002: done."<<std::endl;
  }

  void test003()
  {
    const double target_time=100;
    const double G=1.0;
    //const double G=1.0e-10;

    std::cout<<"test003:"<<std::endl;
    size_t loops;

    for(double dt:{1e-6, 1e-7})
    //for(double dt:{1e-4, 1e-5, 1e-6, 1e-7, 1e-8})
    {
      Runner runner(dt,G);

      runner.set_variation(0.01); // 1%

      std::cout<<"dt = "<<dt<<", estimated "<<(target_time/dt)<<" steps."<<std::endl;
      std::cout<<"\ttime is "<< runner.time()<<std::endl;
      std::cout<<"\tcentre of mass = "<< runner.centre_of_mass()<<", total PE = "<<runner.total_pe()<<", total energy = "<<runner.total_energy()<<std::endl;
      //assert(runner.total_energy() == -G*769./60.);
      std::cout<<"\ttotal energy = PE = -769/60, a.k.a. "<<(-G*769./60.)<<std::endl;
      //runner.dump_masses(std::cout);

      loops=0;
      while(runner.time()<target_time)
      {
        if(runner.step_with_energy_check())
        {
          break;
        }
        ++loops;
      }

      std::cout<<"\ttime is "<< runner.time()<<", "<<loops<<" steps."<<std::endl;
      std::cout<<"\tcentre of mass = "<< runner.centre_of_mass()<<", total PE = "<<runner.total_pe()<<", total energy = "<<runner.total_energy()<<std::endl;
      runner.dump_masses(std::cout);
      std::cout<<std::endl;
      std::cout<<std::endl;
    }

    std::cout<<"test003: done."<<std::endl;
  }

  void test004()
  {
    const double target_time=100;
    const double G=1.0;
    //const double G=1.0e-10;

    std::cout<<"test004:"<<std::endl;
    size_t loops;

    for(double dt:{1e-6, 1e-7})
    //for(double dt:{1e-4, 1e-5, 1e-6, 1e-7, 1e-8})
    {
      Runner runner(dt,G);

      runner.set_variation(0.01); // 1%

      std::cout<<"dt = "<<dt<<", estimated "<<(target_time/dt)<<" steps."<<std::endl;
      std::cout<<"\ttime is "<< runner.time()<<std::endl;
      std::cout<<"\tcentre of mass = "<< runner.centre_of_mass()<<", total PE = "<<runner.total_pe()<<", total energy = "<<runner.total_energy()<<std::endl;
      //assert(runner.total_energy() == -G*769./60.);
      std::cout<<"\ttotal energy = PE = -769/60, a.k.a. "<<(-G*769./60.)<<std::endl;
      //runner.dump_masses(std::cout);

      loops=0;
      while(runner.time()<target_time)
      {
        if(runner.step_with_collision_check())
        {
          break;
        }
        ++loops;
      }

      std::cout<<"\ttime is "<< runner.time()<<", "<<loops<<" steps."<<std::endl;
      std::cout<<"\tcentre of mass = "<< runner.centre_of_mass()<<", total PE = "<<runner.total_pe()<<", total energy = "<<runner.total_energy()<<std::endl;
      runner.dump_masses(std::cout);
      std::cout<<std::endl;
      std::cout<<std::endl;
    }

    std::cout<<"test004: done."<<std::endl;
  }

}; // anonymous namespace

void test_particle_001()
{
  //test001();
  //test002();
  //test003();
  test004();
}

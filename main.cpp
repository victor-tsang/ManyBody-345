//
// Created by Victor Tsang on 2023/10/26.
//

#include<iostream>
#include<chrono>

#include<fmt/format.h>
#include<fmt/chrono.h>
#include<fmt/ostream.h>

#include<Aliz/AlizChrono.hpp>

#include"Particle.hpp"

// Particle.cpp
void test_particle_001();

// HalleyComet.cpp
namespace halley_comet
{
  void test();
};

// for estimating or comparing hardware performance
void test001(double a, double b, double c, size_t array_length)
{
  using Clock=std::chrono::steady_clock;
  using Timer=Aliz::ChronoTimer<Clock>;

  Timer timer;
  timer.start("allocation");

  double *array_A=new double[array_length];
  double *array_B=new double[array_length];
  double *array_C=new double[array_length];

  memset(array_A,0,sizeof(double)*array_length);
  memset(array_B,0,sizeof(double)*array_length);
  memset(array_C,0,sizeof(double)*array_length);

  double dt=1e-10;
  double dx=1e-10;
  double pi = 2.0 * asin(1.0);
  //std::cout<<"Pi = "<<pi<<std::endl;
  double two_pi_dt = pi*2*dt;

  // a += sin(2 pi i dt) - b
  // b += cos(2 pi i dt) + 0.25 * a^2
  // c += sqrt(a) + b

  size_t loops=10;

  timer.stop();
  std::cout<<timer<<std::endl;

  timer.start("computation");
  for(size_t l=0;l<loops;++l)
  {
    for(size_t i=0;i<array_length;++i)
    {
      array_A[i] += sin(two_pi_dt*i) - array_B[i];
      array_B[i] += cos(two_pi_dt*i) + 0.25*array_A[i]*array_A[i];
      array_C[i] += sqrt(array_A[i]) + array_B[i];
    }
  } // loops
  timer.stop();
  std::cout<<timer<<std::endl;

  delete[] array_A;
  delete[] array_B;
  delete[] array_C;
}

namespace test {
  struct vector2 {
    double x, y;

    vector2 operator+(const vector2 &rhs) const {
      return {x + rhs.x, y + rhs.y};
    }

    vector2 operator-(const vector2 &rhs) const {
      return {x - rhs.x, y - rhs.y};
    }

    double norm_square() const {
      return x * x + y * y;
    }

    double norm() const {
      return sqrt(x * x + y * y);
    }

    vector2 operator*(double scalar) const {
      return {scalar * x, scalar * y};
    }

    vector2 operator/(double scalar) const {
      double inv = 1. / scalar;
      return {x * inv, y * inv};
    }

    vector2 operator-() const {
      return {-x, -y};
    }

    double dot(const vector2 &another) const {
      return x * another.x + y * another.y;
    }
  };

  vector2 operator*(const double scalar, const vector2 &vector) {
    return {vector.x * scalar, vector.y * scalar};
  }

  struct particle
  {
    double x, y;
    double vx, vy;
    double ax, ay;
    double mass;
  };

  vector2 displacement(const particle &from, const particle &to) {
    return {to.x - from.x, to.y - from.y};
  }


  // centre of mass of these 3 particles are (0,0).
  const particle M1{1., 3., 0, 0, 0, 0, 3.};
  const particle M2{-2., -1., 0, 0, 0, 0, 4.};
  const particle M3{1., -1., 0, 0, 0, 0, 5.};


  vector2 centre_of_mass(const particle &m1, const particle &m2, const particle &m3)
  {
    double inv_sum_of_masses=1.0/(m1.mass + m2.mass + m3.mass);
    return {  (m1.x*m1.mass + m2.x*m2.mass + m3.x*m3.mass) * inv_sum_of_masses,
              (m1.y*m1.mass + m2.y*m2.mass + m3.y*m3.mass) * inv_sum_of_masses};
  }

}; // namespace test

void test002()
{
  using namespace test;
  particle m1{1.,3.,0,0,0,0,3.};
  particle m2{-2.,-1.,0,0,0,0,4.};
  particle m3{1.,-1.,0,0,0,0,5.};

  auto centre_of_mass=[&m1, &m2, &m3]() -> vector2 {
    vector2 ret{0,0};
    ret.x = (m1.x*m1.mass + m2.x*m2.mass + m3.x*m3.mass) / (m1.mass + m2.mass + m3.mass);
    ret.y = (m1.y*m1.mass + m2.y*m2.mass + m3.y*m3.mass) / (m1.mass + m2.mass + m3.mass);

    return ret;
  };

  auto cm=centre_of_mass();
  std::cout<<fmt::format("center of mass = ({0},{1}).",cm.x, cm.y)<<std::endl;
}


void test003(size_t steps, const double dt, const double G)
{
  using namespace test;

  auto m1=M1;
  auto m2=M2;
  auto m3=M3;
  vector2 cm{0,0};
  vector2 r12{displacement(m1,m2)};
  vector2 r13{displacement(m1,m3)};
  vector2 r23{displacement(m2,m3)};
  double l12sq{r12.norm_square()};
  double l13sq{r13.norm_square()};
  double l23sq{r23.norm_square()};
  double l12{r12.norm()};
  double l13{r13.norm()};
  double l23{r23.norm()};

  auto update_acceleration=[&](){
    // force without masses
    vector2 f12=r12*(G/l12sq/l12);
    vector2 f13=r13*(G/l13sq/l13);
    vector2 f23=r23*(G/l23sq/l23);

    vector2 a1=f12*m2.mass+f13*m3.mass;
    vector2 a2=m3.mass*f23-m1.mass*f12;
    vector2 a3=-m1.mass*f13-m2.mass*f23;

    m1.ax=a1.x;
    m1.ay=a1.y;
    m2.ax=a2.x;
    m2.ay=a2.y;
    m3.ax=a3.x;
    m3.ay=a3.y;
  };

  update_acceleration();

  double t{0};
  double half_dt=0.5*dt;

  for(size_t i=0;i<steps;++i)
  {
    vector2 v1_half{m1.vx + half_dt*m1.ax, m1.vy+half_dt*m1.ay};
    vector2 v2_half{m2.vx + half_dt*m2.ax, m2.vy+half_dt*m2.ay};
    vector2 v3_half{m3.vx + half_dt*m3.ax, m3.vy+half_dt*m3.ay};

    m1.x += v1_half.x*dt;
    m1.y += v1_half.y*dt;
    m2.x += v2_half.x*dt;
    m2.y += v2_half.y*dt;
    m3.x += v3_half.x*dt;
    m3.y += v3_half.y*dt;

    update_acceleration();

    m1.vx = v1_half.x + half_dt * m1.ax;
    m1.vy = v1_half.y + half_dt * m1.ay;
    m2.vx = v2_half.x + half_dt * m2.ax;
    m2.vy = v2_half.y + half_dt * m2.ay;
    m3.vx = v3_half.x + half_dt * m3.ax;
    m3.vy = v3_half.y + half_dt * m3.ay;

  } // for steps

  cm = centre_of_mass(m1,m2,m3);
  std::cout<<fmt::format("After {2} steps, with dt={3}, G={4}, the centre of mass = ({0},{1}).",cm.x, cm.y, steps, dt,G)<<std::endl;
}

namespace test_particle
{
  const Particle M1{{ 1, 3},{0,0},{0,0},3};
  const Particle M2{{-2,-1},{0,0},{0,0},4};
  const Particle M3{{ 1,-1},{0,0},{0,0},5};


  void test001()
  {
    Particle m1{{1,3},{0,0},{0,0},3};

  }
}; // namespace test_particle

int main(int argc,char **argv)
{
  std::cout<<"Hello, world!"<<std::endl;

  //test001(3., 4., 5., 1000000);
  //test002();
  //test003(10000000,1e-5,1.0);
  //test003(10000000,1e-6,1.0);
  //test003(10000000,1e-7,1.0);
  //test003(10000000,1e-8,1.0);
  //test003(10000000,1e-9,1.0);
  //test003(10000000,1e-10,1.0);

  test_particle_001();

  //halley_comet::test();

  return 0;
}
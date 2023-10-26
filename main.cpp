//
// Created by Victor Tsang on 2023/10/26.
//

#include<iostream>
#include<chrono>

#include<fmt/format.h>
#include<fmt/chrono.h>
#include<fmt/ostream.h>

#include<Aliz/AlizChrono.hpp>

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


int main(int argc,char **argv)
{
  std::cout<<"Hello, world!"<<std::endl;

  test001(3., 4., 5., 1000000);

  return 0;
}
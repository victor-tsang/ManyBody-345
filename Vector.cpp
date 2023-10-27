//
// Created by Victor Tsang on 2023/10/27.
//

#include"Vector.hpp"

#include<iostream>


namespace
{

  void test001()
  {
    Vector2 a;
    Vector2 b(1,2);
    Vector2 c=-a;

    std::cout<<"a = "<<a<<";\n"
             <<"b = "<<b<<";\n"
             <<"c = "<<c<<"."<<std::endl;
  }


};
//
// Created by Victor Tsang on 2023/10/27.
//

#ifndef MANYBODY_345_VECTOR_HPP
#define MANYBODY_345_VECTOR_HPP

#include<type_traits>
#include<ostream>

template<typename T>
struct basic_vector2
{
  using value_type=T;

  value_type x;
  value_type y;

  basic_vector2() noexcept
  : basic_vector2(0, 0)
  {
  }

  basic_vector2(const value_type x, const value_type y) noexcept
  : x(x), y(y)
  {}

  ~basic_vector2()=default;

  basic_vector2 operator+(const basic_vector2 &another) const
  {
    return {x+another.x, y+another.y};
  }

  basic_vector2 operator-(const basic_vector2 &another) const
  {
    return {x-another.x, y-another.y};
  }

  basic_vector2 operator-() const
  {
    return {-x,-y};
  }

  basic_vector2 operator*(const value_type &scalar) const
  {
    return {x*scalar, y*scalar};
  }

  basic_vector2 operator+=(const basic_vector2 &another)
  {
    x += another.x;
    y += another.y;
    return *this;
  }

  basic_vector2 operator*=(const value_type &scalar)
  {
    x *= scalar;
    y *= scalar;
    return *this;
  }

  value_type dot(const basic_vector2 &another) const
  {
    return x*another.x + y*another.y;
  }

  value_type norm_square() const
  {
    return x*x+y*y;
  }

  value_type norm() const
  {
    return sqrt(x*x+y*y);
  }

  value_type norm_cubic() const
  {
    auto ret=x*x+y*y;
    return ret*sqrt(ret);
  }

  value_type inverse_norm_cubic() const
  {
    auto ret=x*x+y*y;
    return 1.0/(ret*sqrt(ret));
  }

};

using Vector2 = basic_vector2<double>;

template<typename T>
std::ostream& operator<<(std::ostream &os, const basic_vector2<T> &value) noexcept
{
  os<<'('<<value.x<<','<<value.y<<')';
  return os;
}

#endif //MANYBODY_345_VECTOR_HPP

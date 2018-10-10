#pragma once

#include<iostream>

using std::cout;
using std::endl;
using std::cerr;

//#include<cstdint>

#ifndef TRUE
# define TRUE                          (1U)
#endif

#ifndef FALSE
# define FALSE                         (0U)
#endif

//#include <glm/glm.hpp>

//#include <limits.h>

namespace sv {

#if 0
typedef int int8;
typedef unsigned int uint8;
typedef int int16;
typedef unsigned int uint16;
typedef long int32;
typedef unsigned long uint32;
typedef float real32;
typedef double real64;
#endif

#if 0
typedef glm::vec2 Point;
typedef glm::vec3 Vec3;
#else
//using Vec3 = glm::vec3;
#endif

//class Real3 : public glm::tmat3x3<T, P>
//class Vec3 : public glm::tmat3x3<T, P>
//typedef glm::vec3 Vec3;
//using Real3 = glm::vec<Real, 3>;
//typedef glm::tvec3<float, glm::precision::highp> Real3;

} // namespace sv

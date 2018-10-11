//
//  svpoint.h
//  SphericalVoronoi
//
//  Created by Home on 2014-04-17.
//  Copyright (c) 2014 whenitsdone.org. All rights reserved.
//

#ifndef SphericalVoronoi_svpoint_h
#define SphericalVoronoi_svpoint_h

#include <cassert>
#include <ostream>
#include <sstream>
#include <vector>
#include <tuple>
#include <bitset>
#include <numeric> // for std::accumulate
//#include <static_cast>
#include <math.h>       /* atan2 */

#if 0
#ifndef GLM_FORCE_RADIANS
#define GLM_FORCE_RADIANS
#endif


#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#define GLM_ENABLE_EXPERIMENTAL //NOTE: For using glm::transform, maybe deprecated?
#include <glm/gtx/transform.hpp>

using namespace glm;
#endif

//#include "svPrefix.h"
#include "xs_Float.h"  // http://stereopsis.com/sree/fpu2006.html

#include "svBase.h"
#include "svColor.h"

#define FLOOR_TO_INT(val)      xs_FloorToInt(val)
#define ROUND_TO_INT(val)      xs_RoundToInt(val)
#define CEIL_TO_INT(val)       xs_CeilToInt(val)

#include <limits>


#define FLT_MAX std::numeric_limits<float>::max()
#define FLT_MIN std::numeric_limits<float>::min()
//#std::numeric_limits<float>::infinity();

#include <Magnum/Math/Vector2.h>
#include <Magnum/Math/Vector3.h>
#include <Magnum/Math/Vector4.h>
#include <Magnum/Math/Matrix4.h>
#include <Magnum/Math/Functions.h>
#include <Magnum/Math/Geometry/Distance.h>

namespace sv
{

  typedef Magnum::Math::Vector2<float> Real2;
  typedef Magnum::Math::Vector3<float> Real3;
  typedef Magnum::Math::Vector4<float> Real4;
  typedef Magnum::Math::Matrix4<float> Mat4;

  //typedef Magnum::Math::Vector2<float> F2;
  //typedef Magnum::Math::Vector3<float> F3;
  //typedef Magnum::Math::Vector4<float> F4;

    
  inline std::string stringFromReal3(const Real3& r3)
  {
    std::stringstream ss;
    ss << "<" << r3.x() << "," << r3.y() << "," << r3.z() << ">";
    return ss.str();
  }
    
  //    typedef float Real;
  //    typedef glm::vec3 Real3;
  //    typedef glm::vec4 Real4;
  //    typedef glm::mat4 Mat4;
  //
  //    typedef glm::dvec3 Double3;
  //    typedef glm::dvec4 Double4;
  //    typedef glm::dmat4 DMat4;
    
  class Ray;
  class AABB;
  class Plane;
  class Point;
    
  inline Real3 getTransformPosition(const Mat4& transform) {
    return Real3(transform[3][0], transform[3][1], transform[3][2]);
  }
    
  enum ECubeFace
  {
    CF_POSX = 0,
    CF_FIRST = 0,
    CF_NEGX,
    CF_POSY,
    CF_NEGY,
    CF_POSZ,
    CF_NEGZ,
        
    CF_MAX,
        
    CF_POSX_BITMASK = (1 << CF_POSX),
    CF_NEGX_BITMASK = (1 << CF_NEGX),
    CF_POSY_BITMASK = (1 << CF_POSY),
    CF_NEGY_BITMASK = (1 << CF_NEGY),
    CF_POSZ_BITMASK = (1 << CF_POSZ),
    CF_NEGZ_BITMASK = (1 << CF_NEGZ),
        
    CF_INVALID = 0xFF
  };
    
  typedef std::bitset<CF_MAX> CubeFaceBitSet;
    
  class Point
  {
 public:
 Point() : theta(0), phi(0) { computePosition(); }
        
 Point(float theta_, float phi_)
     : theta(theta_), phi(phi_)
    {
      computePosition();
    }
        
    Point(const Real3& direction)
    {
      assignDirection(direction);
    }
        
    Point(float x, float y, float z)
    {
      assignDirection(Real3(x, y, z));
    }
        

    Magnum::Math::Rad<float> theta;
    Magnum::Math::Rad<float> phi;

    //float phi;
    Real3 position;
        
    std::tuple<Real3, CubeFaceBitSet> cubeCoord() const;
        
    Real3 tangent() const;
    Real3 binormal() const;
        
    void assignDirection(const Real3& direction)
    {
      //float r = Magnum::Math::length(direction);
      float r = direction.length();
      assert(r > 0);
      theta = Magnum::Math::acos(Magnum::Math::clamp(direction.z() /r , -1.0f, 1.0f));

      //TODO: Couldnt make Magnum::Math::atan accept two values and return Rad<float>
      //phi = Magnum::Math::atan(direction.y(), direction.x());
      //phi = Magnum::Math::atan(direction.y());
      phi = Magnum::Math::Rad<float>(atan2(direction.y(), direction.x()));

      position = direction / r;
    }

    float sphericalDistance(const Point& p2) const
    {
      float dot = Magnum::Math::dot(position, p2.position);
      float result = static_cast<float>(Magnum::Math::acos(Magnum::Math::clamp<float>(dot, -1.0, 1.0)));
      return result;
    }

    bool equalWithEps(const Point& p2, float eps) const
    {
      return Magnum::Math::abs(position.x() - p2.position.x()) < eps &&
          Magnum::Math::abs(position.y() - p2.position.y()) < eps &&
          Magnum::Math::abs(position.z() - p2.position.z()) < eps;
    }
        
    bool operator < (const Point& p2)
    {
      return (theta < p2.theta) || (theta == p2.theta && phi < p2.phi);
    }
        
    friend std::ostream& operator << (std::ostream& stream, const Point& p)
    {
      return stream << static_cast<float>(p.theta) << "," << static_cast<float>(p.phi);
    }
        
    template <class Archive>
        void serialize(Archive& ar)
    {
      ar(theta, phi, position.x(), position.y(), position.z());
    }
        
     private:
        void computePosition()
        {
          using namespace Magnum::Math;

          position = Real3{sin(theta) * cos(phi),
                           sin(theta) * sin(phi),
                           cos(theta)};
        }
  };
    
      namespace Util
      {
        extern std::vector<Real3> splitSphericalLineSegment(const Point& start, const Point& end, float deltaAngle = M_PI / 180.0);
        extern float lagrangeInterpolate(float x, const std::vector<float>& xArray, const std::vector<float>& yArray);
        extern float interpolateSphericalSamples(const Point& p0, const std::vector<Point>& points, const std::vector<float>& values);
        extern float computeTriangleArea(const Real3& p0, const Real3& p1, const Real3& p2);
        
        extern void faceAxisDirection(ECubeFace face, Real3& s_dir, Real3& t_dir, Real3& p_dir);

        inline float sqrDistance(Real2 a, Real2 b)
        {
          auto p = a - b;
          return Magnum::Math::dot(p, p);
        }
      }
    
      class SphericalLine
      {
     public:
     SphericalLine() : direction(Real3(0, 0, 1)), xi(0) {}
     SphericalLine(const Real3& direction_, Magnum::Math::Rad<float> xi_) : direction(direction_.normalized()), xi(xi_) {}
        
        Real3 direction;
        Magnum::Math::Rad<float> xi;
      };
    
      class AABB
      {
     public:
        AABB();
        AABB(const Real3& p);
        AABB(const Real3& min, const Real3& max);
        
        void reset();
        bool isValid() const;
        bool isEmpty() const;
        
        const Real3& min() const { assert(isValid()); return m_min; }
        const Real3& max() const { assert(isValid()); return m_max; }
        
        Real3 center() const { assert(isValid()); return (m_min + m_max) * 0.5f; }
        Real3 size() const { assert(isValid()); return m_max - m_min; }
        Real3 extent() const { return size() * 0.5f; }
        
        void getMajorVertices(const Real3& direction, Real3& P, Real3& N) const;
        
        void unionWith(const Real3& p);
        void unionWith(const AABB& aabb);
        bool contains(const Real3& p) const;
        bool contains(const AABB& aabb) const;
        
        bool operator == (const AABB& aabb) const;
     private:
        Real3 m_min, m_max;
      };
    
      class Ray
      {
     public:
        Ray();
        Ray(const Real3& origin, const Real3& direction);
        
        inline const Real3& origin() const { return m_origin; }
        inline void setOrigin(const Real3& origin) { m_origin = origin; }
        inline const Real3& direction() const { return m_direction; }
        inline void setDirection(const Real3& direction) { m_direction = direction.normalized(); }
        inline void setNormalizedDirection(const Real3& direction) { m_direction = direction; }
        
     private:
        Real3 m_origin;
        Real3 m_direction;
      };
    
      class Plane
      {
     public:
        Plane();
        Plane(const Plane& other);
        Plane(const Real3& normal, float distance);
        Plane(const Real4& vec);
        Plane(const Real3& a, const Real3& b, const Real3& c);
        
        const Real3& normal() const { return m_normal; }
        void setNormal(const Real3& normal) { m_normal = normal; }
        float distance() const { return m_distance; }
        void setDistance(float distance) { m_distance = distance; }
        
        Plane normalize() const;
        Plane transform(const Mat4 transform) const;
        float distance(const Real3& point) const;
        bool pointOnSide(const Real3& point) const;
        bool lineIntersection(const Real3& ptA, const Real3& ptB, Real3& resultDestination) const;

     private:
        Real3 m_normal;
        float m_distance;
      };
    
      bool threePlanesIntersection(const Plane& planeA, const Plane& planeB, const Plane& planeC, Real3& result);
    
      bool rayAabbIntersection(const Ray& ray, const AABB& aabb);
    
      template <typename T>
          class PositionT
      {
     public:
        using Vec3 = Magnum::Math::Vector3<T>;
        PositionT();
        PositionT(ECubeFace face, T s, T t, T p);
        PositionT(ECubeFace face, const Vec3& stp);
        
        ECubeFace face() const { return m_face; }
        const Vec3& surfacePoint() const { return m_surfacePoint; }
        Vec3 stpCoords() const;
        const Vec3& spacePosition() const { return m_spacePosition; }
        
     private:
        ECubeFace m_face;
        T m_height;
        Vec3 m_surfacePoint;
        Vec3 m_spacePosition;
      };
    
      //typedef PositionT<float> Position;
      typedef PositionT<float> PositionF;
      typedef PositionT<float> Position;
      typedef PositionT<double> PositionD;
    }

#include "svMath.hpp"

#endif

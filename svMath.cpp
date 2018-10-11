//
//  svMath.cpp
//  SphericalVoronoi
//
//  Created by Home on 2014-06-28.
//  Copyright (c) 2014 whenitsdone.org. All rights reserved.
//

#include "svPrefix.h"

#include "svMath.h"

#include <Magnum/Math/Math.h>
#include <Magnum/Math/Functions.h>

namespace sv
{
    std::tuple<Real3, CubeFaceBitSet> Point::cubeCoord() const
    {
      Real3 absPos = Magnum::Math::abs(position);
      Real3 result;
      CubeFaceBitSet faceSet;
      float ratio;
      if (absPos.z() >= Magnum::Math::max(absPos.x(), absPos.y()))
      {
        ratio = 1.0 / absPos.z();
        faceSet |= (position.z() > 0 ? CF_POSZ_BITMASK : CF_NEGZ_BITMASK);
        if (absPos.z() == absPos.x())
        {
          faceSet |= (position.x() > 0 ? CF_POSX_BITMASK : CF_NEGX_BITMASK);
        }
        if (absPos.z() == absPos.y())
        {
          faceSet |= (position.y() > 0 ? CF_POSY_BITMASK : CF_NEGY_BITMASK);
        }
      }
      else if (absPos.y() >= Magnum::Math::max(absPos.x(), absPos.z()))
       {
         ratio = 1.0 / absPos.y();
         faceSet |= (position.y() > 0 ? CF_POSY_BITMASK : CF_NEGY_BITMASK);
         if (absPos.y() == absPos.x())
         {
           faceSet |= (position.x() > 0 ? CF_POSX_BITMASK : CF_NEGX_BITMASK);
         }
         if (absPos.y() == absPos.z())
         {
           faceSet |= (position.z() > 0 ? CF_POSZ_BITMASK : CF_NEGZ_BITMASK);
         }
       }
       else if (absPos.x() >= Magnum::Math::max(absPos.y(), absPos.z()))
       {
         ratio = 1.0 / absPos.x();
         faceSet |= (position.x() > 0 ? CF_POSX_BITMASK : CF_NEGX_BITMASK);
         if (absPos.x() == absPos.y())
         {
           faceSet |= (position.y() > 0 ? CF_POSY_BITMASK : CF_NEGY_BITMASK);
         }
         if (absPos.x() == absPos.z())
         {
           faceSet |= (position.z() > 0 ? CF_POSZ_BITMASK : CF_NEGZ_BITMASK);
         }
       }
       result = position * ratio;
        return std::tie(result, faceSet);
    }
    
    Real3 Point::tangent() const
    {
        CubeFaceBitSet bitSet;
        
        Real3 binormal(1, 0, 0);
        
        if (bitSet.test(CF_POSX))
        {
            binormal = Real3(0, 1, 0);
        }
        else if (bitSet.test(CF_NEGX))
        {
            binormal = Real3(0, 1, 0);
        }
        else if (bitSet.test(CF_POSY))
        {
            binormal = Real3(1, 0, 0);
        }
        else if (bitSet.test(CF_NEGY))
        {
            binormal = Real3(-1, 0, 0);
        }
        else if (bitSet.test(CF_POSZ))
        {
            binormal = Real3(0, 1, 0);
        }
        else if (bitSet.test(CF_NEGZ))
        {
            binormal = Real3(0, 1, 0);
        }
        else
        {
            assert(false);
        }
        
        Real3 normal = position;
        Real3 result = Magnum::Math::cross(binormal, normal).normalized();
        return result;
    }
    
    Real3 Point::binormal() const
    {
        Real3 t = tangent();
        Real3 normal = position;
        Real3 result = Magnum::Math::cross(normal, t).normalized();
        return result;
    }
    
    // http://www.cgafaq.info/wiki/Intersection_of_three_planes
    bool threePlanesIntersection(const Plane& planeA, const Plane& planeB, const Plane& planeC, Real3& result)
    {
        Real3 bcCross = Magnum::Math::cross(planeB.normal(), planeC.normal());
        float denom = Magnum::Math::dot(planeA.normal(), bcCross);
        
        if (denom == 0) {
          result = Real3(0);
          return false;
        }
        else {
          result = (-planeA.distance() * bcCross
                    - planeB.distance() * Magnum::Math::cross(planeC.normal(), planeA.normal())
                    - planeC.distance() * Magnum::Math::cross(planeA.normal(), planeB.normal())) / denom;
          return true;
        }
    }
    
    // http://tavianator.com/2011/05/fast-branchless-raybounding-box-intersections/
    bool rayAabbIntersection(const Ray& ray, const AABB& aabb)
    {
        Real3 n_inv = Real3(1.0) / ray.direction();
        
        double tx1 = (aabb.min().x() - ray.origin().x())*n_inv.x();
        double tx2 = (aabb.max().x() - ray.origin().x())*n_inv.x();
        
        double tmin = Magnum::Math::min(tx1, tx2);
        double tmax = Magnum::Math::max(tx1, tx2);
        
        double ty1 = (aabb.min().y() - ray.origin().y())*n_inv.y();
        double ty2 = (aabb.max().y() - ray.origin().y())*n_inv.y();
        
        tmin = Magnum::Math::max(tmin, Magnum::Math::min(ty1, ty2));
        tmax = Magnum::Math::min(tmax, Magnum::Math::max(ty1, ty2));
        
        double tz1 = (aabb.min().z() - ray.origin().z())*n_inv.z();
        double tz2 = (aabb.max().z() - ray.origin().z())*n_inv.z();
        
        tmin = Magnum::Math::max(tmin, Magnum::Math::min(tz1, tz2));
        tmax = Magnum::Math::min(tmax, Magnum::Math::max(tz1, tz2));
        
        return tmax >= Magnum::Math::max(tmin, 0.0);
    }
    
    namespace Util
    {
    
    std::vector<Real3> splitSphericalLineSegment(const Point& start, const Point& end, Magnum::Math::Rad<float> deltaAngle)
    {
      std::vector<Real3> result;
            
            assert(start.position != -end.position);
            
            using namespace Magnum::Math;
            auto direction = Magnum::Math::cross(start.position, end.position).normalized();
            Rad<float> distance = Magnum::Math::acos(Magnum::Math::dot(start.position, end.position));
            
             result.push_back(start.position);
            
             for (auto angle=deltaAngle; angle<distance; angle+=deltaAngle)
             {
               //Mat4 rotation = Magnum::Math::rotate(Mat4(1.0), angle, direction);
               Mat4 rotation = Mat4::rotation(angle, direction);

               //Real3 pos = Magnum::Math::normalize(Real3(rotation * Real4(start.position, 1.0)));
               //Real3 pos = Real3{rotation * Real4(start.position, 1.0)};

               Real3 pos = (rotation * Real4(start.position, 1.f)).xyz().normalized(); //NOTE: .xyz() converts Vec4 to Vec3 as per glm semantics.

               result.push_back(pos);
             }
            
             result.push_back(end.position);
            
            return result;
        }
        
    float lagrangeInterpolate(float x, const std::vector<float>& xArray, const std::vector<float>& yArray)
    {
      assert(xArray.size() == yArray.size());
            
      float sum = 0.0;
      for (unsigned int i = 0; i < xArray.size(); ++i)
       {
         float Xi, Yi;
         Xi = xArray[i];
         Yi = yArray[i];
         float factor = 1.0;
         for (unsigned int j = 0; j < xArray.size(); ++j)
         {
           if (i != j)
           {
             float Xj = xArray[j];
             factor *= (x - Xj) / (Xi - Xj);
           }
         }
         sum += factor * Yi;
       }
       return sum;
    }
        
    float interpolateSphericalSamples(const Point& p0, const std::vector<Point>& points, const std::vector<float>& values)
    {
      float totalSqrDistance = std::accumulate(points.begin(), points.end(), 0.0, [p0](float res, const Point& p) {
                                                                                    float d = static_cast<float>(p.sphericalDistance(p0));
                                                                                    return res + d * d;
                                                                                 });
            
      float sum = 0.0;
      float weight = 0.0;
            
      for (size_t i = 0; i < points.size(); ++i)
      {
        const Point& p = points[i];
        float d = static_cast<float>(p.sphericalDistance(p0));
        float w = (totalSqrDistance - d*d) / totalSqrDistance;
        sum += w * values[i];
        weight += w;
      }
      return sum / weight;
    }
        
    float computeTriangleArea(const Real3& p0, const Real3& p1, const Real3& p2)
    {
      Real3 v12 = p2 - p1;
      Real3 v02 = p2 - p0;
      Real3 v12n = v12.normalized();
      float t = Magnum::Math::dot(v02, v12n);
      Real3 c = p2 - v12n * t;
      //float d = glm::distance(p0, c);
      float d = p0.length() - c.length(); //NOTE: should work? https://glm.g-truc.net/0.9.0/api/a00119.html
      float l12 = v12.length();
      return l12 * d * 0.5;
    }

    void faceAxisDirection(ECubeFace face, Real3& s_dir, Real3& t_dir, Real3& p_dir)
    {
      switch (face)
      {
        case CF_POSX:
          p_dir = Real3(1, 0, 0);
          s_dir = Real3(0, 0, -1);
          t_dir = Real3(0, 1, 0);
                    break;
                case CF_NEGX:
                    p_dir = Real3(-1, 0, 0);
                    s_dir = Real3(0, 0, 1);
                    t_dir = Real3(0, 1, 0);
                    break;
                case CF_POSY:
                    p_dir = Real3(0, 1, 0);
                    s_dir = Real3(0, 0, 1);
                    t_dir = Real3(1, 0, 0);
                    break;
                case CF_NEGY:
                    p_dir = Real3(0, -1, 0);
                    s_dir = Real3(0, 0, 1);
                    t_dir = Real3(-1, 0, 0);
                    break;
                case CF_POSZ:
                    p_dir = Real3(0, 0, 1);
                    s_dir = Real3(1, 0, 0);
                    t_dir = Real3(0, 1, 0);
                    break;
                case CF_NEGZ:
                    p_dir = Real3(0, 0, -1);
                    s_dir = Real3(-1, 0, 0);
                    t_dir = Real3(0, 1, 0);
                    break;
                default:
                    assert(0);
                    p_dir = Real3(1, 0, 0);
                    s_dir = Real3(0, 0, -1);
                    t_dir = Real3(0, 1, 0);
            }
        }
    }
    
}

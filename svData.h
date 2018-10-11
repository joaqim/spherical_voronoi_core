//
//  svData.h
//  SphericalVoronoi
//
//  Created by Xiang Wei on 2014-04-18.
//  Copyright (c) 2014 whenitsdone.org. All rights reserved.
//

#ifndef SphericalVoronoi_svdata_h
#define SphericalVoronoi_svdata_h

#include <memory>
#include <vector>
#include <cassert>
#include <set>
#include <ostream>
#include <iterator>
#include <algorithm>
#include "svMath.h"

namespace sv
{
  class HalfEdge;
  class Cell;
  class Vertex;
  class BeachArc;
  class SiteEvent;
  class CircleEvent;

  typedef std::shared_ptr<Vertex> Vertex_ptr;
  typedef std::shared_ptr<HalfEdge> HalfEdge_ptr;
  typedef std::shared_ptr<Cell> Cell_ptr;

  class Cell
  {
 public:
    Cell (uint32_t i, const Point& p) : index(i), point(p)
    {
    }
    uint32_t index;
    Point point;
    std::vector<HalfEdge_ptr> halfEdges;

    void reset()
    {
      halfEdges.clear();
    }

    friend std::ostream& operator<< (std::ostream& stream, const Cell& c)
    {
      return stream << "<" << c.index << "> " << "(" << c.point << ")";
    }
  };

  class Vertex
  {
 public:
    Vertex (const Point& p, std::shared_ptr<Cell> c0, std::shared_ptr<Cell> c1) : point(p)
    {
      cells.insert(c0);
      cells.insert(c1);
    }
    Vertex (const Point& p, std::shared_ptr<Cell> c0, std::shared_ptr<Cell> c1, std::shared_ptr<Cell> c2) : point(p)
    {
      cells.insert(c0);
      cells.insert(c1);
      cells.insert(c2);
    }
    uint32_t index;
    Point point;
    std::vector<HalfEdge_ptr> halfEdges;
    std::set<Cell_ptr> cells;

    void reset()
    {
      halfEdges.clear();
      cells.clear();
    }

    friend std::ostream& operator<< (std::ostream& stream, const Vertex& v)
    {
      stream << "point(" << v.point << ")";
      stream << " cells<";
      for (auto c : v.cells)
      {
        stream << c->index << ",";
      }
      stream << ">";
      return stream;
    }
  };

  class HalfEdge
  {
 public:
 HalfEdge(std::shared_ptr<Vertex> s, std::shared_ptr<Vertex> e)
     : start(s), end(e)
    {
    }
    uint32_t index;
    Cell_ptr pCell;
    Cell_ptr otherCell;
    Vertex_ptr start;
    Vertex_ptr end;
    HalfEdge_ptr prev;
    HalfEdge_ptr next;
    HalfEdge_ptr twin;

    void reset()
    {
      pCell.reset();
      start.reset();
      end.reset();
      prev.reset();
      next.reset();
      twin.reset();
    }

    friend std::ostream& operator<< (std::ostream& stream, const HalfEdge& e)
    {
      stream << "s: " << *e.start << "e: " << *e.end;
      return stream;
    }
  };

  class BeachArc
  {
 public:
 BeachArc(std::shared_ptr<Cell> cell_)
     : pCell(cell_)
    {
    }

    std::shared_ptr<Cell> pCell;

    std::shared_ptr<CircleEvent> circleEvent;      // the related circle event

    std::shared_ptr<Vertex> startVertex;

    bool operator< (const BeachArc& ba) const
    {
      return pCell->point.phi < ba.pCell->point.phi;
    }

    friend std::ostream& operator<< (std::ostream& stream, const BeachArc& arc)
    {
      stream << "cell " << *arc.pCell;
      if (arc.startVertex)
      {
        stream << "startVertex " << *arc.startVertex;
      }
      else
      {
        stream << "startVertex NONE";
      }
      return stream;
    }
  };

  typedef std::shared_ptr<BeachArc> BeachArc_ptr;
  typedef std::vector<BeachArc_ptr> beach_type;

  class SiteEvent
  {
    using Rad = Magnum::Math::Rad<float>;
 public:
 SiteEvent(std::shared_ptr<Cell> cell_)
     : pCell(cell_)
    {
      theta = pCell->point.theta;
      phi = pCell->point.phi;
    }

    std::shared_ptr<Cell> pCell;
    Rad theta;
    Rad phi;

    bool operator< (const SiteEvent& right) const
    {
      return (theta < right.theta) || (theta == right.theta && phi < right.phi);
    }

    friend std::ostream& operator<< (std::ostream& stream, const SiteEvent& e)
    {
      return stream << *e.pCell;
    }
  };

  class CircleEvent
  {
 public:
 CircleEvent(const BeachArc_ptr& arc_i_, const BeachArc_ptr& arc_j_, const BeachArc_ptr& arc_k_)
     : arc_i(arc_i_), arc_j(arc_j_), arc_k(arc_k_)
    {

      auto pij = cell_i()->point.position - cell_j()->point.position;
      auto pkj = cell_k()->point.position - cell_j()->point.position;
      auto direction = cross(pij, pkj);
      circle_center = Point(direction);
      circle_radius = Magnum::Math::acos(dot(circle_center.position, cell_i()->point.position));
      theta = Magnum::Math::acos(circle_center.position.z()) + circle_radius;
    }

    BeachArc_ptr arc_i;
    BeachArc_ptr arc_j;
    BeachArc_ptr arc_k;

    Cell_ptr cell_i() const { return arc_i->pCell; }
    Cell_ptr cell_j() const { return arc_j->pCell; }
    Cell_ptr cell_k() const { return arc_k->pCell; }

    Point circle_center;
    Magnum::Math::Rad<float> circle_radius;

    Magnum::Math::Rad<float> theta;        // the lowest point on circle

    bool operator< (const CircleEvent& ce) const
    {
      return theta < ce.theta;
    }

    friend std::ostream& operator<< (std::ostream& stream, const CircleEvent& e)
    {
      stream << "[" << e.cell_i()->index << "," << e.cell_j()->index << "," << e.cell_k()->index << "] " << "theta " << static_cast<float>(e.theta);
      return stream;
    }
  };

  typedef std::shared_ptr<CircleEvent> CircleEvent_ptr;

  struct compare_circle_event_priority
  {
    bool operator()(const std::shared_ptr<CircleEvent>& left, const std::shared_ptr<CircleEvent>& right) const
    {
      return *left < *right;
    }
    };

}

#endif

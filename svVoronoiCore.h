//
//  svVoronoiCore.h
//  SphericalVoronoi
//
//  Created by Xiang Wei on 2014-05-03.
//  Copyright (c) 2014 whenitsdone.org. All rights reserved.
//

#ifndef SphericalVoronoi_svvoronoi_h
#define SphericalVoronoi_svvoronoi_h

#include "svMath.h"
#include "svData.h"

#include <functional> // for std::function

namespace sv
{
    class SphericalVoronoiCore
    {
    public:
        SphericalVoronoiCore(const std::vector<Real3>& directions);

        void setDebugMode(bool debugMode) { this->debugMode = debugMode; }

        bool isFinished() const;
        void step(float maxDeltaXi);
        void solve(std::function<void(int)> cb = nullptr);       // step until finished

        const std::vector<HalfEdge_ptr>& getHalfEdges() const { return halfEdges; }
        const std::vector<Vertex_ptr>& getVertices() const { return vertices; }
        const std::vector<Cell_ptr>& getCells() const { return cells; }

   protected:
        bool debugMode;

        void dumpBeachState(std::ostream& stream);

        void finializeGraph();
        void cleanupMiddleVertices();
        void duplicateHalfEdges();
        void bindHalfEdgesToCells();

        beach_type::const_iterator getPrevArcOnBeach(beach_type::const_iterator it) const
        {
            if (it != beach.begin())
            {
                return std::prev(it);
            }
            else
            {
                return std::prev(beach.end());
            }
        }

        beach_type::const_iterator getNextArcOnBeach(beach_type::const_iterator it) const
        {
            auto next = std::next(it);
            if (next == beach.end())
            {
                next = beach.begin();
            }
            return next;
        }

        using Rad = Magnum::Math::Rad<float>;

        bool intersectWithNextArc(beach_type::const_iterator itArc, float xi, Point& oPoint) const;
        bool intersectWithPrevArc(beach_type::const_iterator itArc, float xi, Point& oPoint) const;
        bool intersectWithNextArc(beach_type::const_iterator itArc, Rad xi, Point& oPoint) const;
        bool intersectWithPrevArc(beach_type::const_iterator itArc, Rad xi, Point& oPoint) const;
        void handleSiteEvent(SiteEvent& event);
        void handleCircleEvent(const CircleEvent_ptr& event);

        
        static Point thetaToPoint(float theta, bool positive, float xi, float theta1, float phi1);
        static Point thetaToPoint(Rad theta, bool positive, Rad xi, Rad theta1, Rad phi1);
        static Point phiToPoint(float phi, float xi, float theta1, float phi1);
        static Point phiToPoint(Rad phi, Rad xi, Rad theta1, Rad phi1);
        static bool arcsIntersection(const BeachArc& arc1, const BeachArc& arc2, float xi, Point& oPoint);
        static bool arcsIntersection(const BeachArc& arc1, const BeachArc& arc2, Magnum::Math::Rad<float> xi, Point &oPoint);

        int nbSteps;
        SphericalLine scanLine;

        constexpr static float eps = 1e-5;

        std::vector<HalfEdge_ptr> halfEdges;
        std::vector<Vertex_ptr> vertices;
        std::vector<Cell_ptr> cells;

        beach_type beach;

        bool isArcOnBeach(const BeachArc_ptr& arc) const
        {
            return find(beach.begin(), beach.end(), arc) != beach.end();
        }

        std::vector<SiteEvent> siteEventQueue;
        std::vector<CircleEvent_ptr> circleEventQueue;

        void addNewSiteEvent(const SiteEvent& event)
        {
            using namespace std;
            auto it = lower_bound(siteEventQueue.begin(), siteEventQueue.end(), event);
            siteEventQueue.insert(it, event);
        }

        void addNewCircleEvent(const std::shared_ptr<CircleEvent>& event)
        {
            using namespace std;
            auto it = lower_bound(circleEventQueue.begin(), circleEventQueue.end(), event, compare_circle_event_priority());
            circleEventQueue.insert(it, event);
        }

        void removeCircleEvent(const std::shared_ptr<CircleEvent>& event)
        {
            using namespace std;
            auto it = find(circleEventQueue.begin(), circleEventQueue.end(), event);
            assert(it != circleEventQueue.end());
            circleEventQueue.erase(it);
        }
    };
}

#endif

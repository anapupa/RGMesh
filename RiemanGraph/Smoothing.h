//
// Created by pupa on 2022/3/5.
//
#pragma once

#include  "./RiemannGraph.h"


class Smoothing{
public:
    Smoothing(PointCloud& cloud, const UGraph& graph): cloud_(cloud), graph_(graph){}

    void smooth(size_t iters, float alpha = 0.2) {
        std::vector<PointType> new_pts(cloud_.pts.size() );
        for(size_t i_iter = 0; i_iter < iters; i_iter++) {
            for(size_t vi = 0; vi < new_pts.size(); vi++) {
                float w , ww = 0;
                PointType p_total(0, 0, 0), diff;
                for( auto& j: graph_.neigh_vec(vi) ) {
                    diff = cloud_.pts[j] - cloud_.pts[vi] ;
                    w = 1.0f/diff.norm();
                    ww += w;
                    p_total += diff * w;
                }
                new_pts[vi] = cloud_.pts[vi] + p_total / ww * alpha;
            }
            std::swap(new_pts, cloud_.pts);
        }

    }


    PointCloud&     cloud_;
    const UGraph&         graph_;
};
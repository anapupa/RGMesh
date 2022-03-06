//
// Created by pupa on 2021/12/28.
//

#include <string>
#include <RiemanGraph/RiemannGraph.h>
#include <RiemanGraph/UGraph.h>
#include <RiemanGraph/Smoothing.h>

void vertex_relocation(std::string lidar_prefix_name) {
    std::string lidar_points_file = lidar_prefix_name + "_lidars.ply";
    std::string viewpoints_file = lidar_prefix_name + "_viewpoint.ply";

    PointCloud lidar_pts, view_pts;


    lidar_pts.read_from_ply(lidar_points_file);
    view_pts.read_from_ply(viewpoints_file);

    UGraph graph = lidar_pts.delaunay_neighbor();
    {
        Smoothing smooth(lidar_pts, graph);
        smooth.smooth(2, 0.4);
        graph = lidar_pts.delaunay_neighbor();
        export_graph(lidar_prefix_name+"_graph.ply", lidar_pts, graph);
    }

    {
        auto triangles = lidar_pts.local_mesh();
        export_mesh(lidar_prefix_name+"_local_mesh.ply", lidar_pts, triangles);
    }






//    {
//            PointRieGraph point_rg(pointData, 12, false);
//            point_rg.compute_variance();
////            NormalFiltering normalFiltering(point_rg);
////            normalFiltering.flip_normal(viewpoint);
////            normalFiltering.export_ply(lidar_prefix_name );
//    }
}

int main(int argc, char** argv) {
    vertex_relocation(argv[1]);
}

//
// Created by pupa on 2022/3/4.
//

#include "RiemannGraph.h"
#include <happly.h>

#define JC_VORONOI_IMPLEMENTATION
#define JCV_REAL_TYPE double
#include "jc_voronoi.h"
#include <set>



const float eps = std::numeric_limits<double>::epsilon()*20;
bool greater(const jcv_point& lhs, const jcv_point& rhs) {
    return std::abs(lhs.x - rhs.x) < eps ? ( lhs.y - rhs.y > eps): (lhs.x - rhs.x  > eps); // 2078 5954
//    return lhs.x == rhs.x ? (lhs.y > rhs.y): lhs.x > rhs.x; // 4490 5948
}


std::vector<int> delaunay_neighbor_(MatrixX2rd& points_) {
    jcv_diagram diagram_{nullptr};
    jcv_diagram_generate(points_.rows(), (jcv_point*)points_.data(), nullptr, 0, &diagram_);
    std::map<jcv_point, int, decltype(greater)*> V_map(greater);
    const jcv_site* sites = jcv_diagram_get_sites(&diagram_);
    for (int i = 0; i < diagram_.numsites; ++i)   V_map[sites[i].p] = i;

    std::vector<int> neighbor_index;
    for( int i = 0; i < diagram_.numsites; ++i ){
//        std::cout << i << ' ' << sites[i].index << "> " << sites[i].p.x << ' ' << sites[i].p.x << " . " << points_.row(sites[i].index) << std::endl;
        if(sites[i].index != points_.rows()-1 ) continue;
        for(auto graph_edge = sites[i].edges; graph_edge; graph_edge = graph_edge->next){
            if(graph_edge->neighbor == nullptr) continue;
//            std::cout << graph_edge->neighbor->index << ' ';
            neighbor_index.push_back(graph_edge->neighbor->index);
//            auto it = V_map.find(graph_edge->neighbor->p);
//            if(it != V_map.end())  {
//                std::cout << it->second << ' ';
//                neighbor_index.push_back(it->second);
//            }
        }

    }


//    std::cout << " >>>  " << neighbor_index.size() <<  " " <<  points_.rows() << std::endl;

    jcv_diagram_free( &diagram_ );
    return  neighbor_index;
}

Eigen::MatrixXi  PointCloud::KNeighborIndex(size_t k) {
    int pointNum = pts.size();
    Eigen::MatrixXi knn( pts.size(), k);

    // 1. build kd-tree
    typedef nanoflann::L2_Simple_Adaptor<float, PointCloud > AdaptorType;
    typedef nanoflann::KDTreeSingleIndexAdaptor<AdaptorType, PointCloud, 3/*dim*/ > my_kd_tree_t;
    my_kd_tree_t index(3 /*dim*/, *this, nanoflann::KDTreeSingleIndexAdaptorParams(16 ));
    index.buildIndex();

    std::vector<float> sq_dis(k+1);
    std::vector<size_t> adjacent(k+1);

    for (size_t i = 0; i < pointNum; ++i) {
        index.knnSearch(pts[i].data(), k+1, adjacent.data(), sq_dis.data());
        for(size_t j = 1; j <= k; j++)
            knn(i, j-1) = adjacent[j];
    }
    return knn;
}

void PointCloud::read_from_ply(std::string filepath) {
    happly::PLYData ply_in(filepath);
    std::vector<std::array<double, 3>> points = ply_in.getVertexPositions();
//    if(ply_in.getElement("vertex").hasProperty("intensity"))
//        this->intensity = ply_in.getElement("vertex").getProperty<u_char>("intensity");
//    if(ply_in.getElement("vertex").hasProperty("red"))
//        this->colors = ply_in.getVertexColors();
    this->pts.resize(points.size(), Eigen::Vector3f(0,0,0));
    for(int i = 0 ; i < points.size(); i++)
        this->pts[i] = Eigen::Vector3f(points[i][0], points[i][1], points[i][2]);
}

void PointCloud::export_subset_ply(std::string filename, const std::vector<size_t>& v_subset) {
    happly::PLYData ply_out;
    size_t N = v_subset.size();
    N = (N == 0 ? pts.size(): N);
    std::vector<std::array<double, 3>> position( N );
    std::vector<float> normal_x, normal_y, normal_z;
    for (size_t i = 0; i < N; i++) {
        Eigen::Vector3f p = pts[v_subset.size() ? v_subset[i]: i];
        position[i]    = {p.x(), p.y(), p.z()};
    }
    ply_out.addVertexPositions(position);
//    if(colors.size()) {
//        auto _color =  v_subset.size() ? PointCloud::subset(colors, v_subset):colors;
//        ply_out.addVertexColors(_color);
//    }
//    if(intensity.size()) {
//        auto _intensity = v_subset.size() ?  PointCloud::subset(intensity, v_subset): intensity;
//        ply_out.getElement("vertex").addProperty("intensity", _intensity);
//    }
//    if(normals.size()) {
//        normal_x.resize(N);
//        normal_y.resize(N);
//        normal_z.resize(N);
//
//        for (size_t i = 0; i < N; i++) {
//            Eigen::Vector3f p = normals[v_subset.size() ? v_subset[i]: i];
//            normal_x[i] = p.x();
//            normal_y[i] = p.y();
//            normal_z[i] = p.z();
//        }
//
////        std::cout << position.size() << ' ' << normal_x.size() << std::endl;
//        ply_out.getElement("vertex").addProperty("nx", normal_x);
//        ply_out.getElement("vertex").addProperty("ny", normal_y);
//        ply_out.getElement("vertex").addProperty("nz", normal_z);
//    }

    if(filename.find(".ply") == std::string::npos)
        filename += ".ply";
    ply_out.write(filename, happly::DataFormat::Binary);
}

UGraph PointCloud::delaunay_neighbor() {
    UGraph graph;
    graph.add_nodes(pts.size());
    Eigen::MatrixXi knn = KNeighborIndex(22);
    for(size_t vi = 0; vi < pts.size(); vi++) {
        for(auto& j: this->delaunay_neighbor(vi, knn.row(vi)) )
            graph.add_edge(vi, j);
    }
    return graph;
}


std::vector< EigenLet<int,3> >  PointCloud::local_mesh() {
    UGraph graph = delaunay_neighbor();
    graph.sort_edge();
    std::set< EigenLet<int, 3> > triangles;
    for(EdgeType<float>& edge : graph.unique_edges<float>()) {
        for(auto& k: graph.neigh_vec(edge.j)) {
            if(graph.find_edge(edge.i, k))
                triangles.insert(EigenLet<int,3> (edge.i, edge.j, k).sort() );
        }
    }
    return std::vector(triangles.begin(), triangles.end());
}

std::vector<int> PointCloud::delaunay_neighbor(size_t i, const Eigen::RowVectorXi& knn_index) {
    Eigen::MatrixX3f X(knn_index.size()+1, 3);
    for(size_t j = 0; j < knn_index.size(); j++) X.row(j) = pts[knn_index[j]];

    X.bottomRows(1) = pts[i];
    Eigen::Vector3f centroid = X.colwise().mean();
    X.rowwise() -= centroid.transpose();
    Eigen::VectorXf dis = 1.0f/(std::pow(X.rowwise().norm().mean(),2)+X.rowwise().squaredNorm().array());
    dis /= dis.sum();
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigen_solver(  X.transpose() * (dis.asDiagonal() * X) );

    Eigen::Vector3f eig_value  =  eigen_solver.eigenvalues();
    Eigen::Matrix3f eig_vector =  eigen_solver.eigenvectors();

    MatrixX2rd points2(knn_index.size()+1, 2);
    points2.row(0) *= 0;
    for(size_t j = 0; j < knn_index.size(); j++) {
        PointType p_vec = pts[knn_index[j]] - pts[i];
        points2(j, 0) = eig_vector.col(1).dot(p_vec);
        points2(j, 1) = eig_vector.col(2).dot(p_vec);
    }

//    std::ofstream  obj_f(std::to_string(i)+".obj");
//    for(size_t j = 0; j < points2.rows(); j++) {
//        obj_f << "v " << points2(j, 0) <<' ' << points2(j, 1) << " 0 \n";
//    }

    std::vector<int> index = delaunay_neighbor_(points2);
    for (size_t j = 0; j < index.size(); j++)
        index[j] = knn_index[index[j]];

    return index;
}

void export_graph(std::string filename, const PointCloud& cloud,  UGraph& u_graph)  {
    happly::PLYData ply_out;
    size_t N = cloud.pts.size();
    std::vector<std::array<double, 3>> position( N );
    for (size_t i = 0; i < N; i++) {
        Eigen::Vector3f p = cloud.pts[i];
        position[i]    = {p.x(), p.y(), p.z()};
    }
    ply_out.addVertexPositions(position);

    auto edges_ = u_graph.unique_edges<float>();
    std::vector<int> vertex1(edges_.size()), vertex2(edges_.size());
    for(size_t i = 0; i < edges_.size(); i++) {
        vertex1[i] = edges_[i].i;
        vertex2[i] = edges_[i].j;
    }
    ply_out.addElement("edge", edges_.size());
    ply_out.getElement("edge").addProperty("vertex1", vertex1);
    ply_out.getElement("edge").addProperty("vertex2", vertex2);
    ply_out.write(filename, happly::DataFormat::ASCII);
}


void export_mesh(std::string filename, const PointCloud& cloud,  std::vector< EigenLet<int,3> >& triangles)  {
    happly::PLYData ply_out;
    size_t N = cloud.pts.size();
    std::vector<std::array<double, 3>> position( N );
    for (size_t i = 0; i < N; i++) {
        Eigen::Vector3f p = cloud.pts[i];
        position[i]    = {p.x(), p.y(), p.z()};
    }
    ply_out.addVertexPositions(position);

    std::vector< std::vector<int> > polygons;
    polygons.reserve(triangles.size());
    for(auto& triangle: triangles) {
        polygons.push_back({triangle.x(), triangle.y(), triangle.z()});
    }
    ply_out.addFaceIndices(polygons);
    ply_out.write(filename, happly::DataFormat::ASCII);
}

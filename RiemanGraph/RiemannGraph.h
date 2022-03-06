#pragma once

#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Eigenvalues> // header file
#include <vector>
#include <unordered_set>
#include <iostream>
#include <numeric>
#include <nanoflann.hpp>
#include "UGraph.h"

#define PointType Eigen::RowVector3f
typedef Eigen::Matrix<double, -1, 2, Eigen::RowMajorBit> MatrixX2rd;

std::vector<int> delaunay_neighbor_(MatrixX2rd& points_);


template<typename T, int Dim>
class EigenLet: public Eigen::Matrix<T, Dim, 1> {
public:
    using Data = Eigen::Matrix<T, Dim, 1>;
    using Eigen::Matrix<T, Dim, 1>::Matrix;
    virtual ~EigenLet() = default;

    EigenLet<T, Dim>& unique() {
        auto max_it = std::max_element(EigenLet::data(), EigenLet::data()+EigenLet::size());
        std::rotate(EigenLet::data(), max_it, EigenLet::data()+EigenLet::size());
        return *this;
    }

    EigenLet<T, Dim>& sort() {
        auto max_it = std::max_element(EigenLet::data(), EigenLet::data()+EigenLet::size());
        std::sort(EigenLet::data(), EigenLet::data()+EigenLet::size());
        return *this;
    }

    bool operator<(const EigenLet<T, Dim>& other) const {
        for (size_t i=0; i<Dim; i++) {
            if ((*this)[i] < other[i]) return true;
            else if ((*this)[i] > other[i]) return false;
        }
        return false;
    }

    int hash() const {
        static Eigen::RowVector4i p(73856093, 19349663, 83492791, 100663319);
        return p.leftCols(Dim).dot((*this));
    }
};

template<typename KeyType>
struct EigenLetHash {
    inline int operator() (const KeyType& key) const {
        return key.hash();
    }
};


struct PointCloud{
    std::vector<PointType>                  pts;
    PointCloud &operator = (const PointCloud &info) {
        this->pts    = info.pts;
        return *this;
    }
    // Must return the number of data points
    inline size_t kdtree_get_point_count() const { return pts.size(); }

    // Returns the dim'th component of the idx'th point in the class:
    // Since this is inlined and the "dim" argument is typically an immediate value, the
    //  "if/else's" are actually solved at compile time.
    inline double kdtree_get_pt(const size_t idx, int dim) const {
        return dim < pts[idx].size() ? pts[idx][dim]: 0;
    }

    // Optional bounding-box computation: return false to default to a standard bbox computation loop.
    //   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
    //   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /* bb */) const { return false; }

    // IO
    std::vector<std::array<double, 3>> points() const{
        std::vector<std::array<double, 3>> position( pts.size() );
        for (size_t i = 0; i < pts.size(); i++)
            position[i] = {pts[i].x(), pts[i].y(), pts[i].z()};
        return position;
    }

    template<class T>
    std::vector<T> static subset(const std::vector<T>& Src, const std::vector<size_t>& subset) {
        std::vector<T> subset_(subset.size());
        for (size_t i = 0; i < subset.size(); i++)
            subset_[i] = Src[subset[i]];
        return subset_;
    }

    void export_subset_ply(std::string filename, const std::vector<size_t>& v_subset=std::vector<size_t>()) ;

    void read_from_ply(std::string filename) ;


    Eigen::MatrixXi     KNeighborIndex(size_t k) ;


    std::vector<int>    delaunay_neighbor(size_t i, const Eigen::RowVectorXi& knn_index) ;

    UGraph delaunay_neighbor() ;
    std::vector< EigenLet<int,3> > local_mesh();

};


void export_graph(std::string filename, const PointCloud& cloud,  UGraph& u_graph)  ;
void export_mesh(std::string filename, const PointCloud& cloud,  std::vector< EigenLet<int,3> >& triangles);
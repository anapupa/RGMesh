//
// Created by pupa on 2022/1/12.
//
#pragma once
#include <unordered_set>
#include <vector>
#include <set>
#include <unordered_map>
#include <algorithm>



template<class T>
struct EdgeType{
    size_t i, j;
    T value;
    bool operator > (const EdgeType& rhs) const {
        return i == rhs.i ? j > rhs.j : i > rhs.i;
    }
    bool operator < (const EdgeType& rhs) const {
        return i == rhs.i ? j < rhs.j : i < rhs.i;
    }
    bool operator ==(const EdgeType& rhs) const {
        return i == rhs.i && j == rhs.j;
    }
};


struct UGraph{
    void inline add_node() {
        graph_.emplace_back();
        index_map_.emplace_back(index_map_.size());
    }
    void inline add_nodes(size_t n_nodes) {
        graph_.reserve(n_nodes+graph_.size());
        index_map_.reserve(n_nodes+index_map_.size());
        for(size_t i = 0; i < n_nodes; i++) add_node();
    }

    void reserve_degree(size_t k){
        for(size_t i = 0; i < num_nodes(); i++) graph_.reserve(k);
    }

    void unique(){
        for(size_t i = 0; i < num_nodes(); i++) {
            std::sort(graph_[i].begin(), graph_[i].end());
            graph_[i].erase( std::unique(graph_[i].begin(), graph_[i].end()), graph_[i].end() );
        }
    }

    void add_edges(const std::vector<size_t>& source, const std::vector<size_t>& target) {
        assert(source.size() == target.size());
        std::vector<size_t> v_degree(num_nodes(), 0);
        for(size_t i = 0; i < source.size(); i++){
            v_degree[source[i]] ++;
            v_degree[target[i]] ++;
        }
        for(size_t i = 0; i < source.size(); i++)
            graph_[i].reserve(graph_[i].size()+v_degree[i]);

        for(size_t i = 0; i < source.size(); i++) add_edge(source[i], target[i]);
    }

    void add_edges(const Eigen::VectorXi& source, const Eigen::VectorXi& target) {
        assert(source.size() == target.size());
        std::vector<size_t> v_degree(num_nodes(), 0);
        for(size_t i = 0; i < source.size(); i++){
            v_degree[source[i]] ++;
            v_degree[target[i]] ++;
        }
        for(size_t i = 0; i < source.size(); i++)
            graph_[i].reserve(graph_[i].size()+v_degree[i]);

        for(size_t i = 0; i < source.size(); i++) add_edge(source[i], target[i]);
    }

    void add_edge(size_t i, size_t j) {
        if(i == j) return ;
        graph_[i].push_back(j);
        graph_[j].push_back(i);
    }

    std::vector<size_t> neigh_vec(size_t s) const {
        return graph_[s];
    }

    size_t num_nodes() const { return graph_.size(); }

    size_t num_edges() const {
        return std::accumulate(graph_.begin(), graph_.end(), 0,
                               []( size_t s, const std::vector<size_t>& e){ return s + e.size();})/2;
    }

    void sort_edge() {
        for(size_t s = 0; s < graph_.size(); s++)
            std::sort(graph_[s].begin(), graph_[s].end());
    }

    bool find_edge(size_t vi, size_t vj) const {
        return std::binary_search(graph_[vi].begin(), graph_[vi].end(), vj);
    }

    template<class T>
    std::vector<EdgeType<T>> unique_edges() const{
        std::vector<EdgeType<T>> E;
        E.reserve(num_edges());
        for(size_t s = 0; s < graph_.size(); s++) {
            for(size_t t: graph_[s]) {
                if(s > t) continue;
                E.push_back({s, t, T()});
            }
        }
        return E;
    }

    std::vector<size_t>                index_map_;
protected:
    std::vector< std::vector<size_t> > graph_;

};



/*
 * Notice that the list of included headers has
 * expanded a little. As before, you are not allowed
 * to add to this.
 */
#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <array>
#include <list>
#include <forward_list>
#include <deque>
#include <map>
#include <cstddef>
#include <string>
#include <utility>
#include <algorithm>
#include <limits>
#include <optional>
#include <exception>
#include <stdexcept>

#include "directed_graph.hpp"
# define INF 9000000

template <typename vertex>
std::pair<directed_graph<vertex>, std::list<vertex>> kahns_algorithm(const directed_graph<vertex> & d){
  directed_graph<vertex> graph(d);
  std::vector<vertex> zero_indegrees;
  std::list<vertex> topological_order;

  for(auto i : graph)
    if(d.in_degree(i) == 0) zero_indegrees.push_back(i);

  while(!zero_indegrees.empty()){
    vertex v = zero_indegrees.back();
    zero_indegrees.pop_back();
    topological_order.push_back(v);

    for(auto u = d.nbegin(v); u != d.nend(v); ++u){
      graph.remove_edge(v, *u);
      
      if(graph.in_degree(*u) == 0) zero_indegrees.push_back(*u);
    }
  }
  
  return std::make_pair(graph, topological_order);
}

/*
 * Computes whether the input is a Directed Acyclic Graph (DAG).
 * A digraph is a DAG if there is no vertex that has a cycle.
 * A cycle is a non-empty set of [out-]edges that starts at one 
 * vertex, and returns to it.
 */
template <typename vertex>
bool is_dag(const directed_graph<vertex> & d) {
  // If it doesn't have any vertices, it's acyclic.
  if(d.num_vertices()) {
    directed_graph<vertex> g = kahns_algorithm(d).first;

    for(auto i : g) 
      if(g.in_degree(i)) return false;
  }
  
  return true;
}

/*
 * Computes a topological ordering of the vertices.
 * For every vertex u in the order, and any of its
 * neighbours v, v appears later in the order than u.
 */
template <typename vertex>
std::list<vertex> topological_sort(const directed_graph<vertex> & d) {

  if(is_dag(d)) return kahns_algorithm(d).second;
  
  return std::list<vertex>();
}

/*
 * Given a DAG, computes whether there is a Hamiltonian path.
 * a Hamiltonian path is a path that visits every vertex
 * exactly once.
 */
template <typename vertex>
bool is_hamiltonian_dag(const directed_graph<vertex> & d) {  
  std::list<vertex> ts(topological_sort(d));
  
  for(auto i = ts.begin() ; i != ts.end(); ++i){
    auto nx = std::next(i, 1);
    
    if(nx != ts.end() && !d.adjacent(*i, *nx)) return false;
  }
  
  return true;
}

/*
 * Computes the weakly connected components of the graph.
 * A [weak] component is the smallest subset of the vertices
 * such that the in and out neighbourhood of each vertex in
 * the set is also contained in the set.
 */
template <typename vertex>
std::vector<std::vector<vertex>> components(const directed_graph<vertex> & d) {
  
  return std::vector<std::vector<vertex>>();
}

/*
 * Computes the strongly connected components of the graph.
 * A strongly connected component is a subset of the vertices
 * such that for every pair u, v of vertices in the subset,
 * v is reachable from u and u is reachable from v.
 */

template <typename vertex>
std::vector<std::vector<vertex>> strongly_connected_components(const directed_graph<vertex> & d) {
  std::stack<vertex> s;
  int index = 0;
  
  for(auto i : d){
  }
  
  return std::vector<std::vector<vertex>>();
}

template <typename vertex>
vertex minimum_dist(std::unordered_map<vertex, std::size_t> & stp_set, 
                    std::unordered_map<vertex, bool> & visited){
  
  std::size_t min_val = INF;
  vertex min_index;
  
  for(auto& it: stp_set)
    if(visited[it.first] == false && it.second <= min_val)
      min_index = it.first, min_val = it.second;

  return min_index;
}

/*
 * Computes the shortest distance from u to every other vertex
 * in the graph d. The shortest distance is the smallest number
 * of edges in any path from u to the other vertex.
 * If there is no path from u to a vertex, set the distance to
 * be the number of vertices in d plus 1.
 */
template <typename vertex>
std::unordered_map<vertex, std::size_t> shortest_distances(const directed_graph<vertex> & d, const vertex & u) {
  // Dijkstra's algorithm
  
  int size = d.num_vertices();
  std::unordered_map<vertex, std::size_t> stp_set;
  std::unordered_map<vertex, bool> visited;
  //std::priority_queue<vertex, std::vector<vertex>, std::greater<vertex>> pq;

  for(auto& it: d)
    stp_set[it] = INF, visited[it] = false;

  stp_set[u] = 0;
  
  for(auto i = 0; i < size-1 ; i++){
  
    auto v = minimum_dist(stp_set, visited);
    
    visited[v] = true;
    
    for(auto& it: d){
      if(!visited[it] && d.adjacent(v, it) && stp_set[v] != INF && stp_set[v]+1 < stp_set[it]) 
        stp_set[it] = stp_set[v]+1;
      else if(stp_set[it] == INF) stp_set[it] = d.num_vertices()+1;
    }
      
  }
  
  return stp_set;
}

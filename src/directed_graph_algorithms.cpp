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

#define INF 9000000
#define un_map_int std::unordered_map<vertex, int>
#define un_map_bool std::unordered_map<vertex, bool>

template <typename vertex>
std::pair<directed_graph<vertex>, std::list<vertex>> kahns(const directed_graph<vertex> & d){
  directed_graph<vertex> g(d);
  std::queue<vertex> q;
  std::list<vertex> t_order;

  for(auto i : d)
    if(d.in_degree(i) == 0) q.push(i);

  while(!q.empty()){
    vertex v = q.front();
    q.pop();
    t_order.push_back(v);

    for(auto u = d.nbegin(v); u != d.nend(v); ++u){
      g.remove_edge(v, *u);
      
      if(g.in_degree(*u) == 0) q.push(*u);
    }
  }
  
  return std::make_pair(g, t_order);
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
    directed_graph<vertex> g = kahns(d).first;

    for(auto& i : g) 
      if(g.in_degree(i)) 
        return false;
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

  if(is_dag(d)) return kahns(d).second;
  
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
  
  for(auto i = ts.begin() ; i != ts.end(); i++){
    auto nx = std::next(i, 1);
    
    if(nx != ts.end() && !d.adjacent(*i, *nx)) return false;
  }
  
  return true;
}

template <typename vertex>
void dft(const directed_graph<vertex> & d, const vertex& u, 
         std::unordered_map<vertex, bool>& visited, std::vector<vertex>& s_comp){
  
  visited[u] = true;
  s_comp.push_back(u);
  
  for(auto i = d.nbegin(u); i != d.nend(u); i++)
    if(!visited[*i]) 
      dft(d, *i, visited, s_comp);
}

/*
 * Computes the weakly connected components of the graph.
 * A [weak] component is the smallest subset of the vertices
 * such that the in and out neighbourhood of each vertex in
 * the set is also contained in the set.
 */
template <typename vertex>
std::vector<std::vector<vertex>> components(const directed_graph<vertex> & d) {
  directed_graph<vertex> g(d);
  std::vector<std::vector<vertex>> weak_comp;
  std::vector<vertex> s_comp;
  std::unordered_map<vertex, bool> visited;

  for(auto& i: g){
    visited.insert(std::make_pair(i, false));
    for(auto v = g.nbegin(i); v != g.nend(i); v++)
      g.add_edge(*v, i);
  }

  for(auto& i: g){
    if(visited[i] == false){
      dft(g, i, visited, s_comp);
      weak_comp.push_back(s_comp), s_comp.clear();
    }
  }
  
  return weak_comp;
}

template <typename vertex>
void tarjans(const vertex& u, int& g_index, un_map_int& index, 
             un_map_int& low, std::stack<vertex>& s, 
             un_map_bool& on_stack, const directed_graph<vertex> & d,
             std::vector<std::vector<vertex>>& strong_comp){
  
  index[u] = g_index;
  low[u] = g_index;
  ++g_index;
  s.push(u);
  on_stack[u] = true;
  
  for(auto i = d.nbegin(u); i != d.nend(u); i++){
    if(index[*i] == INF){
      tarjans(*i, g_index, index, low, s, on_stack, d, strong_comp);
      low[u] = std::min(low[u], low[*i]);
    } else if(on_stack[*i]){
      low[u] = std::min(low[u], index[*i]);
    } 
  }
  
  if(low[u] == index[u]){
    std::vector<vertex> c;
    vertex w;
    do{
      w = s.top(), s.pop();
      on_stack[w] = false;
      c.push_back(w);
    } while (w != u);
    strong_comp.push_back(c);  
  }
}

/*
 * Computes the strongly connected components of the graph.
 * A strongly connected component is a subset of the vertices
 * such that for every pair u, v of vertices in the subset,
 * v is reachable from u and u is reachable from v.
 */
template <typename vertex>
std::vector<std::vector<vertex>> strongly_connected_components(const directed_graph<vertex> & d) {
  std::vector<std::vector<vertex>> strong_comp;
  un_map_int low, index;
  un_map_bool on_stack;
  std::stack<vertex> s;
  int g_index = 0;
  
  for(auto& i: d){
    low.insert(std::make_pair(i, INF));
    index.insert(std::make_pair(i, INF));
    on_stack.insert(std::make_pair(i, false));
  }
  
  for(auto& v: d)
    if(index[v] == INF)
      tarjans(v, g_index, index, low, s, on_stack, d, strong_comp);
  
  return strong_comp;
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
  std::unordered_map<vertex, std::size_t> sd;
  un_map_bool visited;

  for(auto& it: d)
    sd[it] = INF, visited[it] = false;

  sd[u] = 0;
  
  for(auto i = 0; i < size-1 ; i++){
    std::size_t min_val = INF;
    vertex v;

    for(auto& it: sd)
      if(visited[it.first] == false && it.second <= min_val)
        v = it.first, min_val = it.second;
      
    visited[v] = true;
    
    for(auto& it: d){
      if(!visited[it] && d.adjacent(v, it) && sd[v] != INF && sd[v]+1 < sd[it]) 
        sd[it] = sd[v]+1;
      else if(sd[it] == INF) 
        sd[it] = size+1;
    }
  }
  
  return sd;
}

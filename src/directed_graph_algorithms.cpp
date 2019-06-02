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

// used a large number to define a value
// that can be considered null or infinity in the graph
#define INF 9000000

// pre-define some variable types to keep lines and params short
#define unordered_map_int std::unordered_map<vertex, int>
#define unordered_map_bool std::unordered_map<vertex, bool>
#define unordered_map_size_t std::unordered_map<vertex, std::size_t>

/*
 * Kahn's algorithm for topological sorting looks for vertices
 * with no incoming edges. It iteratively deconstructs the graph
 * from zero in-degree vertices to the out degree vertices.
 */
template <typename vertex>
std::pair<directed_graph<vertex>, std::list<vertex>> kahns(const directed_graph<vertex> & d){
  // make a copy of directed graph d 
  directed_graph<vertex> g(d);
  // and initialise a queue for the zero in-degree vertices
  std::queue<vertex> q;
  // and a list to store topological ordering
  std::list<vertex> t_order;

  // look for vertices with no in-degrees and add them to the queue
  for(auto i : d)
    if(d.in_degree(i) == 0) q.push(i);

  // while the queue isn't empty
  while(!q.empty()){

    // take the vertex from the queue 
    // and add it to the tail of the topological ordering
    vertex v = q.front();
    q.pop();
    t_order.push_back(v);

    // remove edge from neighbouring vertices
    for(auto u = d.nbegin(v); u != d.nend(v); ++u){
      g.remove_edge(v, *u);
      
      // if it the neigbouring vertex no in-degrees, add it to the queue
      if(g.in_degree(*u) == 0) q.push(*u);
    }
  }
  
  // return reconstructed graph and topological order of vertices
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
  // perform kahns algorithm for topological sort on the graph
  if(d.num_vertices()) {
    directed_graph<vertex> g = kahns(d).first;

    // if the resulting graph has any edges, then it's not acyclic
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
  // perform kahns algorithm for topological sort on the graph
  // list will be empty if its not a dag
  return kahns(d).second;
}

/*
 * Given a DAG, computes whether there is a Hamiltonian path.
 * a Hamiltonian path is a path that visits every vertex
 * exactly once.
 */
template <typename vertex>
bool is_hamiltonian_dag(const directed_graph<vertex> & d) {  
  // do a topological sort on the graph
  std::list<vertex> ts(topological_sort(d));
  
  // if the all the vertices in the topological sort order
  // are adjacent to the vertex after it in the list, then it's a hamiltonian dag
  for(auto i = ts.begin() ; i != ts.end(); i++){
    auto nx = std::next(i, 1);
    if(nx != ts.end() && !d.adjacent(*i, *nx)) return false;
  }
  
  return true;
}

/*
 * A normal recursive depth first traversal but this time we
 * have to pass a few more things we use to keep track of connected components.
 * Method is called by components(const directed_graph<vertex> & d).
 */
template <typename vertex>
void dft(const directed_graph<vertex> & d, const vertex& u, 
         unordered_map_bool& visited, std::vector<vertex>& s_comp){
  // set the vertex as visited
  visited[u] = true;
  // add it as part of a single components
  s_comp.push_back(u);
  
  // if the neighbouring vertices of the current vertex 
  // hasn't been visited yet, do a dft on each one too
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
  // make a copy of the graph to manipulate
  directed_graph<vertex> g(d);
  // initialise variables to store single components, 
  // list of all components and visited components
  std::vector<std::vector<vertex>> weak_comp;
  std::vector<vertex> s_comp;
  unordered_map_bool visited;

  // add an edge going backwards to all edges in the graph
  // essentially make it like an undirected graph
  for(auto& i: g){
    visited.insert(std::make_pair(i, false));
    for(auto v = g.nbegin(i); v != g.nend(i); v++)
      g.add_edge(*v, i);
  }

  // then do a normal dft on each vertex on the modified graph
  for(auto& i: g){
    if(visited[i] == false){
      dft(g, i, visited, s_comp);

      // at the end of the dft, add the formed single component
      // to the list of weakly connected components
      weak_comp.push_back(s_comp), s_comp.clear();
    }
  }
  
  return weak_comp;
}

/*
 * Tarjan's algorithm for finding strongly connected components.
 * Method is called by strongly_connected_components(const directed_graph<vertex> & d).
 */
template <typename vertex>
void tarjans(const vertex& u, int& g_index, unordered_map_int& index, 
             unordered_map_int& low, std::stack<vertex>& s, 
             unordered_map_bool& on_stack, const directed_graph<vertex> & d,
             std::vector<std::vector<vertex>>& strong_comp){
  // set the index and low to the current global index 
  index[u] = g_index;
  low[u] = g_index;
  // increment global index and 
  // add + mark that the vertex is on the stack
  ++g_index;
  s.push(u);
  on_stack[u] = true;
  
  // check if all the neighbouring vertices have a low and index value
  for(auto i = d.nbegin(u); i != d.nend(u); i++){

    // if it doesn't the perform tarjan's on the neighbour vertex
    if(index[*i] == INF){
      tarjans(*i, g_index, index, low, s, on_stack, d, strong_comp);

      // set the vertex's low value as either the 
      // low value of the current vertex or the neighbour vertex;
      // whichever one is lower
      low[u] = std::min(low[u], low[*i]);
    } else if(on_stack[*i]){
      // otherwise if the neighbour vertex is already on the stack
      // then set its low value to either the low of the current vertex
      // or the index of the neighbour vertex; whichever one is lower
      low[u] = std::min(low[u], index[*i]);
    } 
  }
  
  // if the low and index values of the current vertex 
  // are the same, then we have a strongly connected component
  if(low[u] == index[u]){
    std::vector<vertex> c;
    vertex w;

    // go through the stack 
    do{
      w = s.top(), s.pop();
      // mark the vertex as no longer being on stack once it's popped of
      on_stack[w] = false;
      // group vertices as part of a single strong component
      c.push_back(w);
    } while (w != u);
    // add the component as a vector to the vector of strong components
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
  unordered_map_int low, index;
  unordered_map_bool on_stack;
  std::stack<vertex> s;
  int g_index = 0;

  // assign a low + index value for each vertex 
  // (we're treating INF here as NULL or at least, a value that it can hopefully never be) 
  // set all the vertices as not being on the stack just yet
  for(auto& i: d){
    low.insert(std::make_pair(i, INF));
    index.insert(std::make_pair(i, INF));
    on_stack.insert(std::make_pair(i, false));
  }
  
  // for each vertex in graph d, if the value isn't null
  // then perform tarjan's algorithm
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
 * 
 * This method uses Dijkstra's algorithm.
 */
template <typename vertex>
unordered_map_size_t shortest_distances(const directed_graph<vertex> & d, const vertex & u) {
  int size = d.num_vertices();
  unordered_map_size_t sd;
  unordered_map_bool visited;

  // set the distance of each vertex to infinity 
  // (it's really just a large number that the distance will hopefully never be...)
  // mark all vertices as unvisited
  for(auto& it: d)
    sd[it] = INF, visited[it] = false;

  // set distance of starting vertex to itself as 0
  sd[u] = 0;
  
  for(auto i = 0; i < size-1 ; i++){
    std::size_t min_val = INF;
    vertex v;

    // select a vertex which hasn't been visited and is of the smallest distance
    // and assign its distance as the minimum value
    for(auto& it: sd)
      if(visited[it.first] == false && it.second <= min_val)
        v = it.first, min_val = it.second;
    
    // mark selected vertex as visited
    visited[v] = true;
    
    // look for the shortest path for all the vertices
    // all edges have a weight of 1 because it's not a weighted graph
    for(auto& it: d){
      // update the distance of the vertex (it) if it hasn't been visited yet,
      // is adjacent to the current selected vertex 
      // and the total weight of the path is smaller than the value of vertex (it).
      if(!visited[it] && d.adjacent(v, it) && sd[v] != INF && sd[v]+1 < sd[it]) 
        sd[it] = sd[v]+1;

      // otherwise if there is no path then set
      // the number of vertices in d plus 1
      else if(sd[it] == INF)
        sd[it] = size+1;
    }
  }
  
  return sd;
}

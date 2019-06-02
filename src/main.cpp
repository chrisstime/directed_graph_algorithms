/*
* Welcome to the playground.
*/
#include <iostream>
#include <climits>

#include "directed_graph_algorithms.cpp"

std::unordered_set<int> generate_random_set(std::size_t size) {
  std::unordered_set<int> set;
  while (set.size() < size) {
    set.insert(rand()%1000);
  }
  return set;
}

template <typename vertex>
std::vector<std::pair<vertex, vertex> > random_tree(const std::vector<vertex>& vertices){
  
  std::vector<std::pair<vertex, vertex> > tree_edges;
  if (!vertices.empty()){
    std::vector<vertex> connected;
    std::vector<vertex> unconnected;
    connected.push_back(vertices[0]);
	
    for (int i = 1; i < vertices.size(); ++i) unconnected.push_back(vertices[i]);
	
    while (connected.size() < vertices.size()){
	
      int index1 = std::rand()%connected.size();
      int index2 = std::rand()%unconnected.size();
      vertex u = connected[index1];
      vertex v = unconnected[index2];
      tree_edges.push_back({u,v});
      unconnected.erase(unconnected.begin() + index2);
      connected.push_back(v);
	
    }
  }
  return tree_edges;

}

int main() {
	auto verts = generate_random_set(5 + rand()%20);
  std::vector<int> ordered_verts(verts.begin(), verts.end());

  directed_graph<int> d;
  for (auto v : ordered_verts) d.add_vertex(v);

  auto tree_edges = random_tree(ordered_verts);

  for (auto e : tree_edges) d.add_edge(e.first, e.second);

  std::cout << is_dag(d) << std::endl;
  for (auto i : d) {
      std::cout << std::string("Vertex: ") << i << std::endl;
      std::cout << std::string("Edges") << std::endl;
      for(auto j = d.nbegin(i); j != d.nend(i); ++j){
        std::cout << *j << std::endl;
      }
      std::cout << std::string("***********") << std::endl;

  }

  std::list<int> ts(topological_sort(d));
  for(auto i = ts.begin() ; i != ts.end(); ++i){
    auto nx = std::next(i, 1);
    bool foo = *i > *nx ;
    std::cout << *nx << std::endl; 
  }
    int foo = 1;
    std::size_t bar = foo;
    std::cout << bar << std::endl;
    std::vector<int> dist(5, 500000);
    std::cout << dist[0] << std::endl;
    
    
    std::unordered_map<std::string, int> thing;
    thing.insert(std::make_pair("5", 6));
    std::cout << thing["5"] << std::endl;
    
}
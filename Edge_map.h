#ifndef Edge_map_h
#define Edge_map_h

struct Polygon;
template <class Triangulation>
struct Adjacent_face {
  std::vector<Polygon>::iterator polygon;
  typename Triangulation::Face_handle face;
  int opposite_vertex;
  
  Adjacent_face(std::vector<Polygon>::iterator polygon, typename Triangulation::Face_handle face, int opposite_vertex) {
    this->polygon = polygon;
    this->face = face;
    this->opposite_vertex = opposite_vertex;
  }
};

template <class Triangulation>
struct Edge {
  std::vector<Adjacent_face<Triangulation>> adjacent_faces;
};

template <class Kernel, class Triangulation>
struct Edge_map {
  std::unordered_map<typename Kernel::Point_2, std::unordered_map<typename Kernel::Point_2, Edge<Triangulation>>> edges;
  
  void insert(typename Kernel::Point_2 &origin, typename Kernel::Point_2 &destination, std::vector<Polygon>::iterator polygon, typename Triangulation::Face_handle face, int opposite_vertex) {
    edges[origin][destination].adjacent_faces.push_back(Adjacent_face<Triangulation>(polygon, face, opposite_vertex));
  }
  
  std::size_t size() const {
    std::size_t number_of_edges = 0;
    for (auto const &edge: edges) {
      number_of_edges += edge.second.size();
    } return number_of_edges;
  }
};


#endif /* Edge_map_h */

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
    CGAL_assertion(face->vertex(face->ccw(opposite_vertex))->point() == origin);
    CGAL_assertion(face->vertex(face->cw(opposite_vertex))->point() == destination);
    edges[origin][destination].adjacent_faces.push_back(Adjacent_face<Triangulation>(polygon, face, opposite_vertex));
  }
  
  void check_consistency() {
    for (auto const &directed_edges_from_source_vertex: edges) {
      for (auto const &adjacent_faces_with_destination_vertex: directed_edges_from_source_vertex.second) {
        for (auto const &adjacent_face: adjacent_faces_with_destination_vertex.second.adjacent_faces) {
          CGAL_assertion(adjacent_face.face->vertex(adjacent_face.face->ccw(adjacent_face.opposite_vertex))->point() == directed_edges_from_source_vertex.first);
          CGAL_assertion(adjacent_face.face->vertex(adjacent_face.face->cw(adjacent_face.opposite_vertex))->point() == adjacent_faces_with_destination_vertex.first);
        }
      }
    }
  }
  
  std::size_t size() const {
    std::size_t number_of_edges = 0;
    for (auto const &directed_edges_from_source_vertex: edges) {
      number_of_edges += directed_edges_from_source_vertex.second.size();
    } return number_of_edges;
  }
  
  void print_info() const {
    std::size_t one_to_zero = 0, many_to_zero = 0, one_to_one = 0, many_to_one = 0, many_to_many = 0;
    for (auto const &directed_edges_from_source_vertex: edges) {
      for (auto const &adjacent_faces_with_destination_vertex: directed_edges_from_source_vertex.second) {
        if (edges.count(adjacent_faces_with_destination_vertex.first) == 0 ||
            edges.at(adjacent_faces_with_destination_vertex.first).count(directed_edges_from_source_vertex.first) == 0) {
          if (adjacent_faces_with_destination_vertex.second.adjacent_faces.size() == 1) {
            ++one_to_zero;
          } else {
            ++many_to_zero;
          }
        } else if (directed_edges_from_source_vertex.first < adjacent_faces_with_destination_vertex.first) {
          if (adjacent_faces_with_destination_vertex.second.adjacent_faces.size() == 1) {
            if (edges.at(adjacent_faces_with_destination_vertex.first).at(directed_edges_from_source_vertex.first).adjacent_faces.size() == 1) {
              ++one_to_one;
            } else {
              ++many_to_one;
            }
          } else {
            if (edges.at(adjacent_faces_with_destination_vertex.first).at(directed_edges_from_source_vertex.first).adjacent_faces.size() == 1) {
              ++many_to_one;
            } else {
              ++many_to_many;
            }
          }
        }
      }
    } std::cout << "Edges: " << std::endl;
    std::cout << "\t" << one_to_zero << "\t1:0 (border)" << std::endl;
    std::cout << "\t" << many_to_zero << "\t2+:0 (border with overlap)" << std::endl;
    std::cout << "\t" << one_to_one << "\t1:1 (perfect match)" << std::endl;
    std::cout << "\t" << many_to_one << "\t2+:1 (overlap on one side)" << std::endl;
    std::cout << "\t" << many_to_many << "\t2+:2+ (overlap on both sides)" << std::endl;
  }
};


#endif /* Edge_map_h */

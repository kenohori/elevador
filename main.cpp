#include <iostream>
#include <fstream>
#include <random>
#include <mach/mach.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include "Enhanced_constrained_triangulation_2.h"
#include <CGAL/Point_set_3.h>
#include <CGAL/random_simplify_point_set.h>
#include <CGAL/wlop_simplify_and_regularize_point_set.h>
#include "Quadtree.h"

#include <ogrsf_frmts.h>

#include <pdal/PointTable.hpp>
#include <pdal/PointView.hpp>
#include <pdal/StageFactory.hpp>
#include <pdal/io/LasReader.hpp>
#include <pdal/io/LasHeader.hpp>
#include <pdal/Options.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Exact_predicates_tag Tag;
struct Vertex_info;
typedef CGAL::Triangulation_vertex_base_with_info_2<Vertex_info, Kernel> Vertex_base;
typedef CGAL::Constrained_triangulation_face_base_2<Kernel> Face_base;
struct FaceInfo;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo, Kernel, Face_base> Face_base_with_info;
typedef CGAL::Triangulation_data_structure_2<Vertex_base, Face_base_with_info> Triangulation_data_structure;
typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, Triangulation_data_structure, Tag> Constrained_delaunay_triangulation;
typedef Enhanced_constrained_triangulation_2<Constrained_delaunay_triangulation> Triangulation;
typedef CGAL::Point_set_3<Kernel::Point_3> Point_cloud;
typedef CGAL::Parallel_if_available_tag Concurrency_tag;
typedef Quadtree<Kernel, Point_cloud> Index;

struct Vertex_info {
  Kernel::FT z;
  Vertex_info() {
    z = 0.0;
  }
};

struct FaceInfo {
  bool processed;
  bool interior;
  FaceInfo() {
    processed = false;
    interior = false;
  }
};

struct Ring {
  std::vector<Kernel::Point_2> points;
};

struct Polygon {
  Ring outer_ring;
  std::vector<Ring> inner_rings;
  Triangulation triangulation;
  std::string cityjson_class;
  Kernel::FT x_min, x_max, y_min, y_max;
};

void printTimer(clock_t &start_time) {
  clock_t stop_time = clock();
  double seconds = (stop_time-start_time)/(double)CLOCKS_PER_SEC;
  std::cout << std::to_string(seconds) + " seconds";
}

void printMemoryUsage() {
  struct task_basic_info t_info;
  mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
  std::string usage;
  if (task_info(mach_task_self(), TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count) == KERN_SUCCESS) {
    
    if (t_info.resident_size > 1024*1024*1024) {
      usage += std::to_string(t_info.resident_size/(1024.0*1024.0*1024.0)) + " GB";
    } else if (t_info.resident_size > 1024*1024) {
      usage += std::to_string(t_info.resident_size/(1024.0*1024.0)) + " MB";
    } else if (t_info.resident_size > 1024) {
      usage += std::to_string(t_info.resident_size/1024.0) + " KB";
    } else {
      usage += std::to_string(t_info.resident_size) + " bytes";
    }
    
  } std::cout << usage << std::endl;
}

int load_map(const char *input_map, std::vector<Polygon> &map_polygons) {
  
  // Prepare input map
  clock_t start_time = clock();
  GDALAllRegister();
  GDALDataset *input_map_dataset = (GDALDataset*) GDALOpenEx(input_map, GDAL_OF_READONLY, NULL, NULL, NULL);
  if (input_map_dataset == NULL) {
    std::cerr << "Error: Could not open input map." << std::endl;
    return 1;
  }
  
  // Load map polygons
  std::cout << "Input map: " << input_map << " type: " << input_map_dataset->GetDriverName() << std::endl;
  for (auto &&input_layer: input_map_dataset->GetLayers()) {
    input_layer->ResetReading();
    OGRFeature *input_feature;
    while ((input_feature = input_layer->GetNextFeature()) != NULL) {
      if (!input_feature->GetGeometryRef()) continue;
      
      std::string cityjson_class = input_feature->GetFieldAsString("citygmlclass");
      if (wkbFlatten(input_feature->GetGeometryRef()->getGeometryType()) == wkbPolygon ||
          wkbFlatten(input_feature->GetGeometryRef()->getGeometryType()) == wkbTriangle) {
        OGRPolygon *input_polygon = input_feature->GetGeometryRef()->toPolygon();
        map_polygons.emplace_back();
        map_polygons.back().cityjson_class = cityjson_class;
        for (int current_vertex = 0; current_vertex < input_polygon->getExteriorRing()->getNumPoints(); ++current_vertex) {
          map_polygons.back().outer_ring.points.emplace_back(input_polygon->getExteriorRing()->getX(current_vertex),
                                                             input_polygon->getExteriorRing()->getY(current_vertex));
        } for (int current_inner_ring = 0; current_inner_ring < input_polygon->getNumInteriorRings(); ++current_inner_ring) {
          map_polygons.back().inner_rings.emplace_back();
          for (int current_vertex = 0; current_vertex < input_polygon->getInteriorRing(current_inner_ring)->getNumPoints(); ++current_vertex) {
            map_polygons.back().inner_rings.back().points.emplace_back(input_polygon->getInteriorRing(current_inner_ring)->getX(current_vertex),
                                                                       input_polygon->getInteriorRing(current_inner_ring)->getY(current_vertex));
          }
        }
      }
      
      else if (wkbFlatten(input_feature->GetGeometryRef()->getGeometryType()) == wkbMultiPolygon) {
        OGRMultiPolygon *input_multipolygon = input_feature->GetGeometryRef()->toMultiPolygon();
        for (int current_polygon = 0; current_polygon < input_multipolygon->getNumGeometries(); ++current_polygon) {
          OGRPolygon *input_polygon = input_multipolygon->getGeometryRef(current_polygon);
          map_polygons.emplace_back();
          map_polygons.back().cityjson_class = cityjson_class;
          for (int current_vertex = 0; current_vertex < input_polygon->getExteriorRing()->getNumPoints(); ++current_vertex) {
            map_polygons.back().outer_ring.points.emplace_back(input_polygon->getExteriorRing()->getX(current_vertex),
                                                               input_polygon->getExteriorRing()->getY(current_vertex));
          } for (int current_inner_ring = 0; current_inner_ring < input_polygon->getNumInteriorRings(); ++current_inner_ring) {
            map_polygons.back().inner_rings.emplace_back();
            for (int current_vertex = 0; current_vertex < input_polygon->getInteriorRing(current_inner_ring)->getNumPoints(); ++current_vertex) {
              map_polygons.back().inner_rings.back().points.emplace_back(input_polygon->getInteriorRing(current_inner_ring)->getX(current_vertex),
                                                                         input_polygon->getInteriorRing(current_inner_ring)->getY(current_vertex));
            }
          }
        }
      }
      
      else {
//        std::cout << "Ignoring non-areal type..." << std::endl;
      }
    }
    
  } GDALClose(input_map_dataset);
  std::cout << "Loaded " << map_polygons.size() << " polygons in ";
  printTimer(start_time);
  std::cout << " using ";
  printMemoryUsage();
  return 0;
}

int triangulate_polygons(std::vector<Polygon> &map_polygons) {
  
  // Basic polygon repair (pre-requisite for triangulation)
  clock_t start_time = clock();
  auto current_polygon = map_polygons.begin();
  while (current_polygon != map_polygons.end()) {
    // Close polygons and rings
    if (current_polygon->outer_ring.points.back() != current_polygon->outer_ring.points.front()) {
      std::cout << "Warning: Last point != first. Adding it again at the end..." << std::endl;
      current_polygon->outer_ring.points.push_back(current_polygon->outer_ring.points.front());
    } for (auto &ring: current_polygon->inner_rings) {
      if (ring.points.back() != ring.points.front()) {
        std::cout << "Warning: Last point != first. Adding it again at the end..." << std::endl;
        ring.points.push_back(ring.points.front());
      }
    }
    
    // Delete degenerate polygons and rings
    if (current_polygon->outer_ring.points.size() < 4) {
      std::cout << "Deleting polygon with < 3 vertices..." << std::endl;
      auto polygon_to_erase = current_polygon;
      ++current_polygon;
      map_polygons.erase(polygon_to_erase);
      continue;
    } auto current_ring = current_polygon->inner_rings.begin();
    while (current_ring != current_polygon->inner_rings.end()) {
      if (current_ring->points.size() < 4) {
        std::cout << "Deleting ring with < 3 vertices..." << std::endl;
        auto ring_to_erase = current_ring;
        ++current_ring;
        current_polygon->inner_rings.erase(ring_to_erase);
        continue;
      } ++current_ring;
    }
    
    ++current_polygon;
  }
  
  // Triangulate polygons
  for (auto &polygon: map_polygons) {
    
    // Triangle (optimisation)
    if (polygon.outer_ring.points.size() == 4 && polygon.inner_rings.size() == 0) {
      if (polygon.outer_ring.points[0] != polygon.outer_ring.points[1] &&
          polygon.outer_ring.points[1] != polygon.outer_ring.points[2] &&
          polygon.outer_ring.points[2] != polygon.outer_ring.points[0]) {
        polygon.triangulation.insert(polygon.outer_ring.points[0]);
        polygon.triangulation.insert(polygon.outer_ring.points[1]);
        polygon.triangulation.insert(polygon.outer_ring.points[2]);
        for (Triangulation::Finite_faces_iterator current_face = polygon.triangulation.finite_faces_begin();
             current_face != polygon.triangulation.finite_faces_end();
             ++current_face) current_face->info().interior = true;
      } else {
        std::cout << "Warning: skipping degenerate triangle..." << std::endl;
      }
    }
    
    // Polygons
    else {
      
      // Triangulate the edges of the polygon
      std::vector<Kernel::Point_2>::const_iterator current_point = polygon.outer_ring.points.begin();
      Triangulation::Vertex_handle current_vertex = polygon.triangulation.insert(*current_point);
      ++current_point;
      Triangulation::Vertex_handle previous_vertex;
      while (current_point != polygon.outer_ring.points.end()) {
        previous_vertex = current_vertex;
        current_vertex = polygon.triangulation.insert(*current_point);
        if (previous_vertex != current_vertex) polygon.triangulation.odd_even_insert_constraint(previous_vertex, current_vertex);
        ++current_point;
      } for (auto const &ring: polygon.inner_rings) {
        current_point = ring.points.begin();
        current_vertex = polygon.triangulation.insert(*current_point);
        while (current_point != ring.points.end()) {
          previous_vertex = current_vertex;
          current_vertex = polygon.triangulation.insert(*current_point);
          if (previous_vertex != current_vertex) polygon.triangulation.odd_even_insert_constraint(previous_vertex, current_vertex);
          ++current_point;
        }
      }
      
      // Label the triangles to find out interior/exterior
      if (polygon.triangulation.number_of_faces() == 0) {
        std::cout << "Degenerate face produced no triangles. Skipping..." << std::endl;
        continue;
      } std::list<Triangulation::Face_handle> to_check;
      polygon.triangulation.infinite_face()->info().processed = true;
      CGAL_assertion(polygon.triangulation.infinite_face()->info().processed == true);
      CGAL_assertion(polygon.triangulation.infinite_face()->info().interior == false);
      to_check.push_back(polygon.triangulation.infinite_face());
      while (!to_check.empty()) {
        CGAL_assertion(to_check.front()->info().processed == true);
        for (int neighbour = 0; neighbour < 3; ++neighbour) {
          if (to_check.front()->neighbor(neighbour)->info().processed == true) {
            // Note: validation code.
//            if (triangulation.is_constrained(Triangulation::Edge(toCheck.front(), neighbour))) CGAL_assertion(toCheck.front()->neighbor(neighbour)->info().interior != toCheck.front()->info().interior);
//            else CGAL_assertion(toCheck.front()->neighbor(neighbour)->info().interior == toCheck.front()->info().interior);
          } else {
            to_check.front()->neighbor(neighbour)->info().processed = true;
            CGAL_assertion(to_check.front()->neighbor(neighbour)->info().processed == true);
            if (polygon.triangulation.is_constrained(Triangulation::Edge(to_check.front(), neighbour))) {
              to_check.front()->neighbor(neighbour)->info().interior = !to_check.front()->info().interior;
              to_check.push_back(to_check.front()->neighbor(neighbour));
            } else {
              to_check.front()->neighbor(neighbour)->info().interior = to_check.front()->info().interior;
              to_check.push_back(to_check.front()->neighbor(neighbour));
            }
          }
        } to_check.pop_front();
      }
    }
      
  } std::cout << "Repaired and triangulated " << map_polygons.size() << " polygons in ";
  printTimer(start_time);
  std::cout << " using ";
  printMemoryUsage();
  return 0;
}

int index_polygons(std::vector<Polygon> &map_polygons, std::unordered_map<Kernel::Point_2, std::set<std::size_t>> &points_index) {
  
  // Compute bounds and index boundary vertices of repaired polygons
  clock_t start_time = clock();
  for (std::size_t current_polygon = 0; current_polygon < map_polygons.size(); ++current_polygon) {
    bool first_boundary_point = true;
    for (auto vertex: map_polygons[current_polygon].triangulation.finite_vertex_handles()) {
      bool incident_to_interior = false;
      bool incident_to_exterior = false;
      Triangulation::Face_circulator first_face = map_polygons[current_polygon].triangulation.incident_faces(vertex);
      Triangulation::Face_circulator current_face = first_face;
      do {
        if (current_face->info().interior) incident_to_interior = true;
        else incident_to_exterior = true;
        ++current_face;
      } while (current_face != first_face);
      if (incident_to_interior && incident_to_exterior) {
        points_index[vertex->point()].insert(current_polygon);
        if (first_boundary_point) {
          map_polygons[current_polygon].x_min = vertex->point().x();
          map_polygons[current_polygon].x_max = vertex->point().x();
          map_polygons[current_polygon].y_min = vertex->point().y();
          map_polygons[current_polygon].y_max = vertex->point().y();
          first_boundary_point = false;
        } else {
          if (vertex->point().x() < map_polygons[current_polygon].x_min) map_polygons[current_polygon].x_min = vertex->point().x();
          if (vertex->point().x() > map_polygons[current_polygon].x_max) map_polygons[current_polygon].x_max = vertex->point().x();
          if (vertex->point().y() < map_polygons[current_polygon].y_min) map_polygons[current_polygon].y_min = vertex->point().y();
          if (vertex->point().y() > map_polygons[current_polygon].y_max) map_polygons[current_polygon].y_max = vertex->point().y();
        }
      }
    }
  }
  
  // Print extent
  Kernel::FT x_min = map_polygons.front().x_min;
  Kernel::FT x_max = map_polygons.front().x_max;
  Kernel::FT y_min = map_polygons.front().y_min;
  Kernel::FT y_max = map_polygons.front().y_max;
  for (auto const &polygon: map_polygons) {
    if (polygon.x_min < x_min) x_min = polygon.x_min;
    if (polygon.x_max > x_max) x_max = polygon.x_max;
    if (polygon.y_min < y_min) y_min = polygon.y_min;
    if (polygon.y_max > y_max) y_max = polygon.y_max;
  } std::cout << "Map extent: X = [" << x_min << ", " << x_max << "] Y = [" << y_min << ", " << y_max << "]" << std::endl;
  
  std::cout << "Indexed " << points_index.size() << " polygon triangulation boundary points in ";
  printTimer(start_time);
  std::cout << " using ";
  printMemoryUsage();
  return 0;
}

int load_point_cloud(const char *input_point_cloud, Point_cloud &point_cloud) {
  clock_t start_time = clock();
  std::cout << "Input point cloud: " << input_point_cloud << std::endl;
  pdal::Options input_opts;
  input_opts.add("filename", input_point_cloud);
  pdal::PointTable input_table;
  pdal::LasReader las_reader;
  las_reader.setOptions(input_opts);
  las_reader.prepare(input_table);
  pdal::PointViewSet point_view_set = las_reader.execute(input_table);
  pdal::PointViewPtr point_view = *point_view_set.begin();
  pdal::Dimension::IdList dims = point_view->dims();
  pdal::LasHeader las_header = las_reader.header();
  
  for (pdal::PointId idx = 0; idx < point_view->size(); ++idx) {
    point_cloud.insert(Kernel::Point_3(point_view->getFieldAs<double>(pdal::Dimension::Id::X, idx),
                                       point_view->getFieldAs<double>(pdal::Dimension::Id::Y, idx),
                                       point_view->getFieldAs<double>(pdal::Dimension::Id::Z, idx)));
    
  }
  
  // Print extent
  std::cout << "Point cloud extent: X = [" << las_header.minX() << ", " << las_header.maxX() << "] Y = [" << las_header.minY() << ", " << las_header.maxY() << "] Z = [" << las_header.minZ() << ", " << las_header.maxZ() << "]" << std::endl;
  
  std::cout << "Loaded " << las_header.pointCount() << " points in ";
  printTimer(start_time);
  std::cout << " using ";
  printMemoryUsage();
  return 0;
}

int index_point_cloud(Point_cloud &point_cloud, Index &index) {
  
  const int bucket_size = 100;
  const int maximum_depth = 10;
  
  // Index point cloud
  clock_t start_time = clock();
  index.compute_extent(point_cloud);
  for (Point_cloud::const_iterator point_index = point_cloud.begin(); point_index != point_cloud.end(); ++point_index) index.insert_point(point_cloud, *point_index);
  index.optimise(point_cloud, bucket_size, maximum_depth);
  
  // Print quadtree info
  index.print_info();
  std::cout << "Indexed " << point_cloud.size() << " point cloud points in ";
  printTimer(start_time);
  std::cout << " using ";
  printMemoryUsage();
  return 0;
}

int create_terrain(std::vector<Polygon> &map_polygons, Point_cloud &point_cloud, Index &index) {
  
//  const double retain_percentage = 1.0;
//  const double neighbour_radius = 0.5;
  
  // Remove points from undesirable classes
  clock_t start_time = clock();
  Point_cloud terrain_point_cloud = point_cloud;
//  for (auto &polygon: map_polygons) {
//    if (polygon.cityjson_class == "Building" ||
//        polygon.cityjson_class == "WaterBody") {
//      std::vector<Index *> intersected_nodes;
//      index.find_intersections(intersected_nodes, polygon.x_min, polygon.x_max, polygon.y_min, polygon.y_max);
//      for (auto const &node: intersected_nodes) {
//        for (auto const &point_index: node->points) {
//          Triangulation::Face_handle face = polygon.triangulation.locate(Kernel::Point_2(point_cloud.point(point_index).x(),
//                                                                                         point_cloud.point(point_index).y()));
//          if (!polygon.triangulation.is_infinite(face) && face->info().interior) terrain_point_cloud.remove(point_index);
//        }
//      }
//    }
//  } std::cout << "Removed " << terrain_point_cloud.garbage_size() << " points from terrain point cloud" << std::endl;
  
  // Thinning
  Point_cloud::iterator first_point_to_remove = CGAL::random_simplify_point_set(terrain_point_cloud, 99.0);
  terrain_point_cloud.remove(first_point_to_remove, terrain_point_cloud.end());
  std::cout << "Thinned point cloud to " << terrain_point_cloud.number_of_points() << " points" << std::endl;
  
  // Simplify
//  CGAL::wlop_simplify_and_regularize_point_set<Concurrency_tag>(point_cloud, terrain_point_cloud.point_back_inserter(),
//                                                                CGAL::parameters::select_percentage(retain_percentage).
//                                                                neighbor_radius(neighbour_radius));
//  std::cout << terrain_point_cloud.size() << " points in terrain point cloud" << std::endl;
  
  // Create TIN
  Triangulation terrain_triangulation;
  for (auto const &point: terrain_point_cloud.points()) {
    Triangulation::Vertex_handle vertex = terrain_triangulation.insert(Kernel::Point_2(point.x(), point.y()));
    vertex->info().z = point.z();
  }
  
  std::cout << "Created terrain TIN with " << terrain_triangulation.number_of_faces() << " triangles in ";
  printTimer(start_time);
  std::cout << " using ";
  printMemoryUsage();
  return 0;
}

int create_buildings(std::vector<Polygon> &map_polygons, Point_cloud &point_cloud, Index &index) {
  clock_t start_time = clock();
  std::size_t n_buildings = 0;
  for (auto &polygon: map_polygons) {
    if (polygon.cityjson_class == "Building") {
      std::vector<Index *> intersected_nodes;
      index.find_intersections(intersected_nodes, polygon.x_min, polygon.x_max, polygon.y_min, polygon.y_max);
      std::vector<Point_cloud::Index> points_in_polygon;
      for (auto const &node: intersected_nodes) {
        for (auto const &point_index: node->points) {
          Triangulation::Face_handle face = polygon.triangulation.locate(Kernel::Point_2(point_cloud.point(point_index).x(),
                                                                                         point_cloud.point(point_index).y()));
          if (!polygon.triangulation.is_infinite(face) && face->info().interior) points_in_polygon.push_back(point_index);
        }
      }
      
      if (!points_in_polygon.empty()) {
        Kernel::FT sum_of_elevations = 0.0;
        for (auto const &point_index: points_in_polygon) sum_of_elevations += point_cloud.point(point_index).z();
        Kernel::FT building_elevation = sum_of_elevations/points_in_polygon.size();
        for (Triangulation::Finite_vertices_iterator current_vertex = polygon.triangulation.finite_vertices_begin();
             current_vertex != polygon.triangulation.finite_vertices_end();
             ++current_vertex) {
          current_vertex->info().z = building_elevation;
        }
      } else {
        // TODO
      }
      
      ++n_buildings;
    }
  }
  
  std::cout << "Created " << n_buildings << " buildings in ";
  printTimer(start_time);
  std::cout << " using ";
  printMemoryUsage();
  return 0;
}

int write_obj(const char *output_3dcm, std::vector<Polygon> &map_polygons) {
  clock_t start_time = clock();
  std::ofstream output_stream;
  output_stream.open(output_3dcm);
  output_stream << "mtllib ./elevador.mtl" << std::endl;
  std::string output_objects;
  std::size_t num_faces = 0;
  for (std::size_t current_polygon = 0; current_polygon < map_polygons.size(); ++current_polygon) {
    output_objects += "o " + std::to_string(current_polygon) + "\n";
    output_objects += "usemtl " + map_polygons[current_polygon].cityjson_class + "\n";
    for (Triangulation::Finite_faces_iterator current_face = map_polygons[current_polygon].triangulation.finite_faces_begin();
         current_face != map_polygons[current_polygon].triangulation.finite_faces_end();
         ++current_face) {
      if (current_face->info().interior == false) continue;
      for (int v = 0; v < 3; ++v) output_stream << "v " << current_face->vertex(v)->point().x() << " " << current_face->vertex(v)->point().y() << " " << current_face->vertex(v)->info().z << "\n";
      output_objects += "f " + std::to_string(3*num_faces+1) + " " + std::to_string(3*num_faces+2) + " " + std::to_string(3*num_faces+3) + "\n";
      ++num_faces;
    }
  } output_stream << output_objects;
  output_stream.close();
  std::cout << "Wrote 3D city model in ";
  printTimer(start_time);
  std::cout << " using ";
  printMemoryUsage();
  return 0;
}

int main(int argc, const char * argv[]) {
  
  const char *input_map = "/Users/ken/Downloads/3dfier_os/osmm/osmm.gpkg";
  const char *input_point_cloud = "/Users/ken/Downloads/3dfier_os/osmm/Exeter_VOLTAtest.laz";
  const char *output_3dcm = "/Users/ken/Downloads/exeter.obj";
  
  std::vector<Polygon> map_polygons;
  std::unordered_map<Kernel::Point_2, std::set<std::size_t>> points_index;
  Point_cloud point_cloud;
  Index index;
  
  load_map(input_map, map_polygons);
  triangulate_polygons(map_polygons);
  index_polygons(map_polygons, points_index);
  load_point_cloud(input_point_cloud, point_cloud);
  index_point_cloud(point_cloud, index);
  create_terrain(map_polygons, point_cloud, index);
  create_buildings(map_polygons, point_cloud, index);
  write_obj(output_3dcm, map_polygons);
  
  return 0;
}

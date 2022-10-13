#include <iostream>
#include <fstream>
#include <algorithm>
#include <locale>
#include <iomanip>
#include <mach/mach.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include "Enhanced_constrained_triangulation_2.h"
#include <CGAL/Point_set_3.h>
#include <CGAL/Barycentric_coordinates_2/triangle_coordinates_2.h>
#include "Quadtree.h"
#include "Edge_map.h"

#include <ogrsf_frmts.h>

#include <pdal/PointTable.hpp>
#include <pdal/PointView.hpp>
#include <pdal/StageFactory.hpp>
#include <pdal/io/LasReader.hpp>
#include <pdal/io/LasHeader.hpp>
#include <pdal/Options.hpp>

#include <nlohmann/json.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Exact_predicates_tag Tag;
struct Vertex_info;
typedef CGAL::Triangulation_vertex_base_with_info_2<Vertex_info, Kernel> Vertex_base;
typedef CGAL::Constrained_triangulation_face_base_2<Kernel> Face_base;
struct Face_info;
typedef CGAL::Triangulation_face_base_with_info_2<Face_info, Kernel, Face_base> Face_base_with_info;
typedef CGAL::Triangulation_data_structure_2<Vertex_base, Face_base_with_info> Triangulation_data_structure;
typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, Triangulation_data_structure, Tag> Constrained_delaunay_triangulation;
typedef Enhanced_constrained_triangulation_2<Constrained_delaunay_triangulation> Triangulation;
typedef CGAL::Point_set_3<Kernel::Point_3> Point_cloud;
typedef Quadtree_node<Kernel, Point_cloud> Point_index;
typedef Edge_map<Kernel, Triangulation> Edge_index;

struct Vertex_info {
  Kernel::FT z;
  Vertex_info() {
    z = 0.0;
  }
};

struct Face_info {
  bool processed;
  bool interior;
  Face_info() {
    processed = false;
    interior = false;
  }
};

struct Ring {
  std::vector<Kernel::Point_2> points;
};

struct Triangle {
  Kernel::Point_3 p1, p2, p3;
  Triangle(Kernel::Point_3 &p1, Kernel::Point_3 &p2, Kernel::Point_3 &p3) {
    this->p1 = p1;
    this->p2 = p2;
    this->p3 = p3;
  }
};

struct Polygon {
  Ring outer_ring;
  std::vector<Ring> inner_rings;
  Triangulation triangulation;
  std::vector<Triangle> extra_triangles;
  std::string cityjson_class;
  std::string id;
  std::unordered_map<std::string, std::string> attributes;
  Kernel::FT x_min, x_max, y_min, y_max;
};

void print_timer(clock_t start_time) {
  clock_t stop_time = clock();
  double seconds = (stop_time-start_time)/(double)CLOCKS_PER_SEC;
  std::cout << seconds << " seconds";
}

void print_memory_usage() {
  struct task_basic_info t_info;
  mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
  if (task_info(mach_task_self(), TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count) == KERN_SUCCESS) {
    if (t_info.resident_size > 1024*1024*1024) {
      std::cout << t_info.resident_size/(1024.0*1024.0*1024.0) << " GB";
    } else if (t_info.resident_size > 1024*1024) {
      std::cout << t_info.resident_size/(1024.0*1024.0) << " MB";
    } else if (t_info.resident_size > 1024) {
      std::cout << t_info.resident_size/1024.0 << " KB";
    } else {
      std::cout << t_info.resident_size << " bytes";
    } std::cout << std::endl;
  }
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
      std::string cityjson_class;
      std::unordered_map<std::string, std::string> attributes;
      
      for (int current_field = 0; current_field < input_feature->GetFieldCount(); ++current_field) {
        if (strcmp(input_feature->GetFieldDefnRef(current_field)->GetNameRef(), "citygmlclass") == 0) {
          cityjson_class = input_feature->GetFieldAsString("citygmlclass");
          continue;
        } switch (input_feature->GetFieldDefnRef(current_field)->GetType()) {
          case OFTReal:
            attributes[input_feature->GetFieldDefnRef(current_field)->GetNameRef()] = std::to_string(input_feature->GetFieldAsDouble(current_field));
            break;
          case OFTString:
            attributes[input_feature->GetFieldDefnRef(current_field)->GetNameRef()] = input_feature->GetFieldAsString(current_field);
            break;
          case OFTInteger:
            attributes[input_feature->GetFieldDefnRef(current_field)->GetNameRef()] = std::to_string(input_feature->GetFieldAsInteger(current_field));
            break;
          case OFTInteger64:
            attributes[input_feature->GetFieldDefnRef(current_field)->GetNameRef()] = std::to_string(input_feature->GetFieldAsInteger64(current_field));
            break;
          case OFTDate:
            int year, month, day, hour, minute, second, timezone;
            input_feature->GetFieldAsDateTime(current_field, &year, &month, &day, &hour, &minute, &second, &timezone);
            attributes[input_feature->GetFieldDefnRef(current_field)->GetNameRef()] = std::to_string(year) + "-" + std::to_string(month) + "-" + std::to_string(day);
            break;
          default:
            std::cout << "Unsupported field type: " << input_feature->GetFieldDefnRef(current_field)->GetType() << std::endl;
            break;
        }
      }
      
      
      
      if (wkbFlatten(input_feature->GetGeometryRef()->getGeometryType()) == wkbPolygon ||
          wkbFlatten(input_feature->GetGeometryRef()->getGeometryType()) == wkbTriangle) {
        OGRPolygon *input_polygon = input_feature->GetGeometryRef()->toPolygon();
        map_polygons.emplace_back();
        map_polygons.back().id = std::to_string(input_feature->GetFID());
        map_polygons.back().cityjson_class = cityjson_class;
        map_polygons.back().attributes = attributes;
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
          map_polygons.back().id = std::to_string(input_feature->GetFID()) + "-" + std::to_string(current_polygon);
          map_polygons.back().cityjson_class = cityjson_class;
          map_polygons.back().attributes = attributes;
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
  print_timer(start_time);
  std::cout << " using ";
  print_memory_usage();
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
      
    // Insert the edges of the polygon as constraints
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
    } if (polygon.triangulation.number_of_faces() == 0) {
      std::cout << "Warning: degenerate polygon (no triangles after insertion of constraints)..." << std::endl;
      continue;
    }
    
    // Label the triangles to find out interior/exterior
    std::list<Triangulation::Face_handle> to_check;
    polygon.triangulation.infinite_face()->info().processed = true;
    CGAL_assertion(polygon.triangulation.infinite_face()->info().processed == true);
    CGAL_assertion(polygon.triangulation.infinite_face()->info().interior == false);
    to_check.push_back(polygon.triangulation.infinite_face());
    while (!to_check.empty()) {
      CGAL_assertion(to_check.front()->info().processed == true);
      for (int neighbour = 0; neighbour < 3; ++neighbour) {
        if (to_check.front()->neighbor(neighbour)->info().processed == true) {
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
    
    // Check if the result isn't degenerate
    std::size_t interior_triangles = 0;
    for (auto const &current_face: polygon.triangulation.finite_face_handles()) {
      if (current_face->info().interior) ++interior_triangles;
    } if (interior_triangles == 0) {
      std::cout << "Warning: degenerate polygon (no interior triangles)..." << std::endl;
      polygon.triangulation.clear();
    }
  }
  
  // Compute bounds of repaired polygons
  for (auto &polygon: map_polygons) {
    bool first_boundary_point = true;
    for (auto current_vertex: polygon.triangulation.finite_vertex_handles()) {
      bool incident_to_interior = false;
      bool incident_to_exterior = false;
      Triangulation::Face_circulator first_face = polygon.triangulation.incident_faces(current_vertex);
      Triangulation::Face_circulator current_face = first_face;
      do {
        if (current_face->info().interior) incident_to_interior = true;
        else incident_to_exterior = true;
        ++current_face;
      } while (current_face != first_face);
      if (incident_to_interior && incident_to_exterior) {
        if (first_boundary_point) {
          polygon.x_min = current_vertex->point().x();
          polygon.x_max = current_vertex->point().x();
          polygon.y_min = current_vertex->point().y();
          polygon.y_max = current_vertex->point().y();
          first_boundary_point = false;
        } else {
          if (current_vertex->point().x() < polygon.x_min) polygon.x_min = current_vertex->point().x();
          if (current_vertex->point().x() > polygon.x_max) polygon.x_max = current_vertex->point().x();
          if (current_vertex->point().y() < polygon.y_min) polygon.y_min = current_vertex->point().y();
          if (current_vertex->point().y() > polygon.y_max) polygon.y_max = current_vertex->point().y();
        }
      }
    }
  }
  
  std::cout << "Repaired and triangulated " << map_polygons.size() << " polygons in ";
  print_timer(start_time);
  std::cout << " using ";
  print_memory_usage();
  return 0;
}

int index_map(std::vector<Polygon> &map_polygons, Edge_index &edges_index) {
  clock_t start_time = clock();
  std::unordered_map<Kernel::Point_2, std::unordered_map<Kernel::Point_2, std::tuple<std::vector<Polygon>::iterator, Triangulation::Face_handle, int>>> open_edges;
  
  // Index edges
  for (std::vector<Polygon>::iterator current_polygon = map_polygons.begin(); current_polygon != map_polygons.end(); ++current_polygon) {
    for (auto const &face: current_polygon->triangulation.finite_face_handles()) {
      if (face->info().interior) {
        for (int opposite_vertex = 0; opposite_vertex < 3; ++opposite_vertex) {
          if (face->neighbor(opposite_vertex) == current_polygon->triangulation.infinite_face() ||
              !face->neighbor(opposite_vertex)->info().interior) {
            Triangulation::Vertex_handle origin = face->vertex(face->ccw(opposite_vertex));
            Triangulation::Vertex_handle destination = face->vertex(face->cw(opposite_vertex));
            CGAL_assertion(CGAL::orientation(origin->point(), destination->point(), face->vertex(opposite_vertex)->point()) == CGAL::COUNTERCLOCKWISE);
            edges_index.insert(origin->point(), destination->point(), current_polygon, face, opposite_vertex);
          }
        }
      }
    }
  }
  
  // Print index stats
//  edges_index.print_info();
  
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
  } // std::cout << "Map extent: X = [" << x_min << ", " << x_max << "] Y = [" << y_min << ", " << y_max << "]" << std::endl;
  
  std::cout << "Indexed " << edges_index.size() << " half-edges edges of " << map_polygons.size() << " polygons in ";
  print_timer(start_time);
  std::cout << " using ";
  print_memory_usage();
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
  // std::cout << "Point cloud extent: X = [" << las_header.minX() << ", " << las_header.maxX() << "] Y = [" << las_header.minY() << ", " << las_header.maxY() << "] Z = [" << las_header.minZ() << ", " << las_header.maxZ() << "]" << std::endl;
  
  std::cout << "Loaded " << las_header.pointCount() << " points in ";
  print_timer(start_time);
  std::cout << " using ";
  print_memory_usage();
  return 0;
}

int index_point_cloud(Point_cloud &point_cloud, Point_index &index) {
  
  const int bucket_size = 100;
  const int maximum_depth = 10;
  
  // Index point cloud
  clock_t start_time = clock();
  index.compute_extent(point_cloud);
  for (Point_cloud::const_iterator point_index = point_cloud.begin(); point_index != point_cloud.end(); ++point_index) index.insert_point(point_cloud, *point_index);
  index.optimise(point_cloud, bucket_size, maximum_depth);
  
  // Print quadtree info
//  index.print_info();
  std::cout << "Indexed " << point_cloud.size() << " point cloud points in ";
  print_timer(start_time);
  std::cout << " using ";
  print_memory_usage();
  return 0;
}

int create_terrain_tin(std::vector<Polygon> &map_polygons, Point_cloud &point_cloud, Point_index &index, Triangulation &terrain) {
  clock_t start_time = clock();
  
  // Remove points from undesirable classes
  Point_cloud terrain_point_cloud = point_cloud;
  for (auto &polygon: map_polygons) {
    if (polygon.cityjson_class == "Building" ||
        polygon.cityjson_class == "WaterBody") {
      std::vector<Point_index *> intersected_nodes;
      index.find_intersections(intersected_nodes, polygon.x_min, polygon.x_max, polygon.y_min, polygon.y_max);
      for (auto const &node: intersected_nodes) {
        for (auto const &point_index: node->points) {
          Triangulation::Face_handle face = polygon.triangulation.locate(Kernel::Point_2(point_cloud.point(point_index).x(),
                                                                                         point_cloud.point(point_index).y()));
          if (!polygon.triangulation.is_infinite(face) && face->info().interior) terrain_point_cloud.remove(point_index);
        }
      }
    }
  } // std::cout << "Used map polygons to remove " << terrain_point_cloud.garbage_size() << " points from terrain point cloud" << std::endl;
  
  // Create grid DTM using nearby low points
  Point_cloud dtm_point_cloud;
  const Kernel::FT dtm_cell_size = 5.0;
  const Kernel::FT dtm_search_radius = 10.0;
  const Kernel::FT dtm_ratio_to_use = 0.01;
  const std::size_t minimum_points_in_grid = 100;
  Kernel::FT squared_search_radius = dtm_search_radius*dtm_search_radius;
  for (Kernel::FT x = index.x_min-0.5*dtm_search_radius; x < index.x_max+0.5*dtm_search_radius; x += dtm_cell_size) {
    for (Kernel::FT y = index.y_min-0.5*dtm_search_radius; y < index.y_max+0.5*dtm_search_radius; y += dtm_cell_size) {
      std::vector<Point_index *> intersected_nodes;
      index.find_intersections(intersected_nodes, x-0.5*dtm_search_radius, x+0.5*dtm_search_radius, y-0.5*dtm_search_radius, y+0.5*dtm_search_radius);
      std::vector<Kernel::FT> elevations;
      for (auto const &node: intersected_nodes) {
        for (auto const &point_index: node->points) {
          Kernel::FT squared_2d_distance = (point_cloud.point(point_index).x()-x)*(point_cloud.point(point_index).x()-x) +
                                           (point_cloud.point(point_index).y()-y)*(point_cloud.point(point_index).y()-y);
          if (!terrain_point_cloud.is_removed(point_index) &&
              squared_2d_distance < squared_search_radius) {
            elevations.push_back(point_cloud.point(point_index).z());
          }
        }
      } if (elevations.size() > minimum_points_in_grid) {
        std::sort(elevations.begin(), elevations.end());
        dtm_point_cloud.insert(Kernel::Point_3(x, y, elevations[std::floor(dtm_ratio_to_use*elevations.size())]));
      }
    }
  } // std::cout << "Created rough DTM with " << dtm_point_cloud.number_of_points() << " points" << std::endl;
  
  // Index DTM
  Point_index dtm_index;
  index_point_cloud(dtm_point_cloud, dtm_index);
  
  // Smoothen
  Point_cloud smooth_dtm_point_cloud;
  const Kernel::FT smoothing_radius = 50.0;
  Kernel::FT squared_smoothing_radius = smoothing_radius*smoothing_radius;
  for (Kernel::FT x = index.x_min-0.5*smoothing_radius; x < index.x_max+0.5*smoothing_radius; x += dtm_cell_size) {
    for (Kernel::FT y = index.y_min-0.5*smoothing_radius; y < index.y_max+0.5*smoothing_radius; y += dtm_cell_size) {

      std::vector<Point_index *> intersected_nodes;
      dtm_index.find_intersections(intersected_nodes, x-smoothing_radius, x+smoothing_radius, y-smoothing_radius, y+smoothing_radius);
      Kernel::FT sum_of_elevations = 0.0;
      std::size_t number_of_points = 0.0;
      std::vector<Kernel::Point_2> points_within_radius;
      for (auto const &node: intersected_nodes) {
        for (auto const &point_index: node->points) {
          Kernel::FT squared_2d_distance = (dtm_point_cloud.point(point_index).x()-x)*(dtm_point_cloud.point(point_index).x()-x) +
                                           (dtm_point_cloud.point(point_index).y()-y)*(dtm_point_cloud.point(point_index).y()-y);
          if (squared_2d_distance <= squared_smoothing_radius) {
            sum_of_elevations += dtm_point_cloud.point(point_index).z();
            ++number_of_points;
          }
        }
      } if (number_of_points > 0) {
        smooth_dtm_point_cloud.insert(Kernel::Point_3(x, y, sum_of_elevations/number_of_points));
      }
    }
  }
  
  // Create TIN from grid
  terrain.clear();
  for (auto const &point: smooth_dtm_point_cloud.points()) {
    Triangulation::Vertex_handle vertex = terrain.insert(Kernel::Point_2(point.x(), point.y()));
    vertex->info().z = point.z();
  }
  
  // TODO: simplify

  std::cout << "Created terrain TIN with " << terrain.number_of_faces() << " triangles in ";
  print_timer(start_time);
  std::cout << " using ";
  print_memory_usage();
  
  return 0;
}

int lift_flat_polygons(std::vector<Polygon> &map_polygons, const char *cityjson_class, Point_cloud &point_cloud, Point_index &index, Kernel::FT ratio_to_use) {
  clock_t start_time = clock();
  std::size_t n_polygons = 0;
  for (auto &polygon: map_polygons) {
    if (polygon.cityjson_class == cityjson_class) {
      
      // Find index nodes overlapping the polygon bbox
      std::vector<Point_index *> intersected_nodes;
      index.find_intersections(intersected_nodes, polygon.x_min, polygon.x_max, polygon.y_min, polygon.y_max);
      
      // Find PC points overlapping the polygon
      std::vector<Point_cloud::Index> points_in_polygon;
      for (auto const &node: intersected_nodes) {
        for (auto const &point_index: node->points) {
          Triangulation::Face_handle face = polygon.triangulation.locate(Kernel::Point_2(point_cloud.point(point_index).x(),
                                                                                         point_cloud.point(point_index).y()));
          if (!polygon.triangulation.is_infinite(face) && face->info().interior) points_in_polygon.push_back(point_index);
        }
      }
      
      // If there are points, use those
      if (!points_in_polygon.empty()) {
        
        // Sort elevations to obtain elevation
        std::vector<Kernel::FT> elevations;
        for (auto const &point_index: points_in_polygon) elevations.push_back(point_cloud.point(point_index).z());
        std::sort(elevations.begin(), elevations.end());
        Kernel::FT polygon_elevation = elevations[std::floor(ratio_to_use*elevations.size())];
        
        // Set elevation of polygon points to calculated elevation
        for (Triangulation::Finite_vertices_iterator current_vertex = polygon.triangulation.finite_vertices_begin();
             current_vertex != polygon.triangulation.finite_vertices_end();
             ++current_vertex) {
          current_vertex->info().z = polygon_elevation;
        }
        
      }
      
      // TODO: If there are no points
      else {
        std::cout << "No points in polygon!" << std::endl;
      }
      
      ++n_polygons;
    }
  }
  
  std::cout << "Lifted " << n_polygons << " flat " << cityjson_class << " polygons in ";
  print_timer(start_time);
  std::cout << " using ";
  print_memory_usage();
  return 0;
}

int lift_polygon_vertices(std::vector<Polygon> &map_polygons, const char *cityjson_class, Triangulation &terrain) {
  clock_t start_time = clock();
  std::size_t n_vertices = 0, n_polygons = 0;
  for (auto &polygon: map_polygons) {
    if (polygon.cityjson_class == cityjson_class) {
      for (Triangulation::Finite_vertices_iterator current_vertex = polygon.triangulation.finite_vertices_begin();
           current_vertex != polygon.triangulation.finite_vertices_end();
           ++current_vertex) {
        
        Triangulation::Face_handle face_of_point = terrain.locate(current_vertex->point());
        std::vector<Kernel::FT> barycentric_coordinates;
        CGAL::Barycentric_coordinates::triangle_coordinates_2(face_of_point->vertex(0)->point(),
                                                              face_of_point->vertex(1)->point(),
                                                              face_of_point->vertex(2)->point(),
                                                              current_vertex->point(),
                                                              std::back_inserter(barycentric_coordinates));
        current_vertex->info().z = barycentric_coordinates[0]*face_of_point->vertex(0)->info().z +
                                   barycentric_coordinates[1]*face_of_point->vertex(1)->info().z +
                                   barycentric_coordinates[2]*face_of_point->vertex(2)->info().z;
        
        ++n_vertices;
      } ++n_polygons;
    }
  }
  
  std::cout << "Lifted " << n_vertices << " vertices of " << n_polygons << " " << cityjson_class << " polygons for in ";
  print_timer(start_time);
  std::cout << " using ";
  print_memory_usage();
  return 0;
}

int lift_polygons(std::vector<Polygon> &map_polygons, const char *cityjson_class, Triangulation &terrain) {
  clock_t start_time = clock();
  std::size_t n_polygons = 0;
  for (auto &polygon: map_polygons) {
    if (polygon.cityjson_class == cityjson_class) {
      
      // Lift polygon vertices
      for (Triangulation::Finite_vertices_iterator current_vertex = polygon.triangulation.finite_vertices_begin();
           current_vertex != polygon.triangulation.finite_vertices_end();
           ++current_vertex) {
        
        
        Triangulation::Face_handle face_of_point = terrain.locate(current_vertex->point());
        std::vector<Kernel::FT> barycentric_coordinates;
        CGAL::Barycentric_coordinates::triangle_coordinates_2(face_of_point->vertex(0)->point(),
                                                              face_of_point->vertex(1)->point(),
                                                              face_of_point->vertex(2)->point(),
                                                              current_vertex->point(),
                                                              std::back_inserter(barycentric_coordinates));
        current_vertex->info().z = barycentric_coordinates[0]*face_of_point->vertex(0)->info().z +
                                   barycentric_coordinates[1]*face_of_point->vertex(1)->info().z +
                                   barycentric_coordinates[2]*face_of_point->vertex(2)->info().z;
      }
      
      // Add terrain vertices in polygon
      std::vector<Triangulation::Vertex_handle> terrain_vertices_in_polygon;
      for (auto const &current_vertex: terrain.finite_vertex_handles()) {
        if (current_vertex->point().x() > polygon.x_min && current_vertex->point().x() < polygon.x_max &&
            current_vertex->point().y() > polygon.y_min && current_vertex->point().y() < polygon.y_max) {
          Triangulation::Locate_type locate_type;
          int vertex_index;
          Triangulation::Face_handle face_of_point = polygon.triangulation.locate(current_vertex->point(), locate_type, vertex_index);
          if (locate_type == Triangulation::FACE && face_of_point->info().interior == true) {
            terrain_vertices_in_polygon.push_back(current_vertex);
          }
        }
      } for (auto const &terrain_vertex: terrain_vertices_in_polygon) {
        Triangulation::Vertex_handle polygon_vertex = polygon.triangulation.insert(terrain_vertex->point());
        polygon_vertex->info().z = terrain_vertex->info().z;
      }
      
      // Relabel triangles as interior/exterior
      for (auto const &current_face: polygon.triangulation.finite_face_handles()) {
        current_face->info().processed = false;
        current_face->info().interior = false;
      } std::list<Triangulation::Face_handle> to_check;
      polygon.triangulation.infinite_face()->info().processed = true;
      polygon.triangulation.infinite_face()->info().interior = false;
      to_check.push_back(polygon.triangulation.infinite_face());
      while (!to_check.empty()) {
        CGAL_assertion(to_check.front()->info().processed == true);
        for (int neighbour = 0; neighbour < 3; ++neighbour) {
          if (to_check.front()->neighbor(neighbour)->info().processed == true) {
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
      
      ++n_polygons;
    }
  }
  std::cout << "Lifted " << n_polygons << " " << cityjson_class << " polygons in ";
  print_timer(start_time);
  std::cout << " using ";
  print_memory_usage();
  return 0;
}

int create_vertical_walls(std::vector<Polygon> &map_polygons, Edge_index &edge_index) {
  clock_t start_time = clock();
  
  // For every border edge in map polygon
  for (std::vector<Polygon>::iterator current_polygon = map_polygons.begin(); current_polygon != map_polygons.end(); ++current_polygon) {
    for (auto const &face: current_polygon->triangulation.finite_face_handles()) {
      if (face->info().interior) {
        for (int opposite_vertex = 0; opposite_vertex < 3; ++opposite_vertex) {
          if (face->neighbor(opposite_vertex) == current_polygon->triangulation.infinite_face() ||
              !face->neighbor(opposite_vertex)->info().interior) {
            
            // Get elevations at origin and destination
            Triangulation::Vertex_handle origin = face->vertex(face->ccw(opposite_vertex));
            Triangulation::Vertex_handle destination = face->vertex(face->cw(opposite_vertex));
            CGAL_assertion(CGAL::orientation(origin->point(), destination->point(), face->vertex(opposite_vertex)->point()) == CGAL::COUNTERCLOCKWISE);
            std::set<Kernel::FT> origin_elevations, destination_elevations;
            for (auto const &same_side_face: edge_index.edges[origin->point()][destination->point()].adjacent_faces) {
              Triangulation::Vertex_handle same_side_origin = same_side_face.face->vertex(same_side_face.face->ccw(same_side_face.opposite_vertex));
              Triangulation::Vertex_handle same_side_destination = same_side_face.face->vertex(same_side_face.face->cw(same_side_face.opposite_vertex));
              CGAL_assertion(same_side_origin->point() == origin->point());
              CGAL_assertion(same_side_destination->point() == destination->point());
              origin_elevations.insert(same_side_origin->info().z);
              destination_elevations.insert(same_side_destination->info().z);
            } for (auto const &opposite_side_face: edge_index.edges[destination->point()][origin->point()].adjacent_faces) {
              Triangulation::Vertex_handle opposite_side_origin = opposite_side_face.face->vertex(opposite_side_face.face->ccw(opposite_side_face.opposite_vertex));
              Triangulation::Vertex_handle opposite_side_destination = opposite_side_face.face->vertex(opposite_side_face.face->cw(opposite_side_face.opposite_vertex));
              CGAL_assertion(opposite_side_origin->point() == destination->point());
              CGAL_assertion(opposite_side_destination->point() == origin->point());
              origin_elevations.insert(opposite_side_destination->info().z);
              destination_elevations.insert(opposite_side_origin->info().z);
            } CGAL_assertion(origin_elevations.count(origin->info().z) == 1);
            CGAL_assertion(destination_elevations.count(destination->info().z) == 1);
            
            // Vertical wall is not needed
            if (origin_elevations.size() == 1 && destination_elevations.size() == 1) continue;

            // Quad-like
            if (origin_elevations.size() == 2 && destination_elevations.size() == 2) {
              
              // Quad (from top to avoid duplicates)
              if (origin->info().z == *origin_elevations.rbegin() && destination->info().z == *destination_elevations.rbegin()) {
                Kernel::Point_3 origin_top(origin->point().x(), origin->point().y(), *origin_elevations.rbegin());
                Kernel::Point_3 destination_top(destination->point().x(), destination->point().y(), *destination_elevations.rbegin());
                Kernel::Point_3 origin_bottom(origin->point().x(), origin->point().y(), *origin_elevations.begin());
                Kernel::Point_3 destination_bottom(destination->point().x(), destination->point().y(), *destination_elevations.begin());
                current_polygon->extra_triangles.push_back(Triangle(destination_top, origin_top, origin_bottom));
                current_polygon->extra_triangles.push_back(Triangle(origin_bottom, destination_bottom, destination_top));
              }
              
              // Bowtie (from origin at top and destination at bottom to avoid duplicates)
              else if (origin->info().z == *origin_elevations.rbegin() && destination->info().z == *destination_elevations.begin()) {
                std::cout << "Unsupported case: bowtie" << std::endl;
                Kernel::Point_3 origin_top(origin->point().x(), origin->point().y(), *origin_elevations.rbegin());
                Kernel::Point_3 destination_top(destination->point().x(), destination->point().y(), *destination_elevations.rbegin());
                Kernel::Point_3 origin_bottom(origin->point().x(), origin->point().y(), *origin_elevations.begin());
                Kernel::Point_3 destination_bottom(destination->point().x(), destination->point().y(), *destination_elevations.begin());
                Kernel::Segment_3 top_to_bottom(origin_top, destination_bottom);
                Kernel::Segment_3 bottom_to_top(origin_bottom, destination_top);
                auto result = CGAL::intersection(top_to_bottom, bottom_to_top);
                if (result) {
                  Kernel::Point_3 *intersection_point = boost::get<Kernel::Point_3>(&*result);
                  
                } else {
                  std::cout << "Error: bowtie intersection point not found" << std::endl;
                }
//                current_polygon->extra_triangles.push_back(Triangle(destination_top, origin_top, origin_bottom));
//                current_polygon->extra_triangles.push_back(Triangle(origin_bottom, destination_bottom, destination_top));
              }
            }
              
            // Triangle (from 2 origin elevations and 1 destination elevation to avoid duplicates)
            else if (origin_elevations.size() == 2 && destination_elevations.size() == 1) {
              Kernel::Point_3 origin_top(origin->point().x(), origin->point().y(), *origin_elevations.rbegin());
              Kernel::Point_3 origin_bottom(origin->point().x(), origin->point().y(), *origin_elevations.begin());
              Kernel::Point_3 destination_only(destination->point().x(), destination->point().y(), destination->info().z);
              current_polygon->extra_triangles.push_back(Triangle(destination_only, origin_top, origin_bottom));
            }
            
            // Skip other triangle case
            else if (origin_elevations.size() == 1 && destination_elevations.size() == 2);
            
            // TODO: Warn for other cases (overlapping polygons)
            else {
              std::cout << "Unsupported case: " << origin_elevations.size() << " origin elevations and " << destination_elevations.size() << " destination elevations" << std::endl;
            }
          }
        }
      }
    }
  }
  
  std::cout << "Created vertical walls in ";
  print_timer(start_time);
  std::cout << " using ";
  print_memory_usage();
  return 0;
}

int write_3dcm_obj(const char *output_3dcm, std::vector<Polygon> &map_polygons) {
  clock_t start_time = clock();
  
  const int decimal_digits = 2;
  
  std::ofstream output_stream;
  output_stream.open(output_3dcm);
  output_stream << std::fixed;
  output_stream << std::setprecision(decimal_digits);
  output_stream << "mtllib ./elevador.mtl" << std::endl;
  std::unordered_map<Kernel::Point_3, std::size_t> output_vertices;
  std::size_t num_polygons = 0;
  
  // Vertices
  for (std::vector<Polygon>::iterator current_polygon = map_polygons.begin(); current_polygon != map_polygons.end(); ++current_polygon) {
    for (Triangulation::Finite_faces_iterator current_face = current_polygon->triangulation.finite_faces_begin();
         current_face != current_polygon->triangulation.finite_faces_end();
         ++current_face) {
      if (current_face->info().interior == false) continue;
      for (int v = 0; v < 3; ++v) {
        Kernel::Point_3 point(current_face->vertex(v)->point().x(), current_face->vertex(v)->point().y(), current_face->vertex(v)->info().z);
        if (output_vertices.count(point) == 0) {
          output_stream << "v " << point.x() << " " << point.y() << " " << point.z() << "\n";
          output_vertices[point] = output_vertices.size()+1;
        }
      }
    } for (auto const &triangle: current_polygon->extra_triangles) {
      if (output_vertices.count(triangle.p1) == 0) {
        output_stream << "v " << triangle.p1.x() << " " << triangle.p1.y() << " " << triangle.p1.z() << "\n";
        output_vertices[triangle.p1] = output_vertices.size()+1;
      } if (output_vertices.count(triangle.p2) == 0) {
        output_stream << "v " << triangle.p2.x() << " " << triangle.p2.y() << " " << triangle.p2.z() << "\n";
        output_vertices[triangle.p2] = output_vertices.size()+1;
      } if (output_vertices.count(triangle.p3) == 0) {
        output_stream << "v " << triangle.p3.x() << " " << triangle.p3.y() << " " << triangle.p3.z() << "\n";
        output_vertices[triangle.p3] = output_vertices.size()+1;
      }
    }
  }
  
  // Faces
  for (std::vector<Polygon>::iterator current_polygon = map_polygons.begin(); current_polygon != map_polygons.end(); ++current_polygon) {
    
    // Triangles in polygon triangulation
    output_stream << "o " << std::to_string(num_polygons) << "\n";
    output_stream << "usemtl " << current_polygon->cityjson_class << "\n";
    for (Triangulation::Finite_faces_iterator current_face = current_polygon->triangulation.finite_faces_begin();
         current_face != current_polygon->triangulation.finite_faces_end();
         ++current_face) {
      if (current_face->info().interior == false) continue;
      output_stream << "f";
      for (int v = 0; v < 3; ++v) {
        output_stream << " " << output_vertices[Kernel::Point_3(current_face->vertex(v)->point().x(), current_face->vertex(v)->point().y(), current_face->vertex(v)->info().z)];
      } output_stream << "\n";
    }
    
    // Additional triangles
    bool first_wall = true;
    for (auto const &triangle: current_polygon->extra_triangles) {
      if (first_wall) {
        if (current_polygon->cityjson_class == "Building") output_stream << "usemtl BuildingWall\n";
        first_wall = false;
      } std::size_t v1 = output_vertices[triangle.p1];
      std::size_t v2 = output_vertices[triangle.p2];
      std::size_t v3 = output_vertices[triangle.p3];
      output_stream << "f " << v1 << " " << v2 << " " << v3 << "\n";
    }
    
  ++num_polygons;
  }
  
  output_stream.close();
  std::cout << "Wrote 3D city model in ";
  print_timer(start_time);
  std::cout << " using ";
  print_memory_usage();
  return 0;
}

int write_3dcm_cityjson(const char *output_3dcm, std::vector<Polygon> &map_polygons) {
  clock_t start_time = clock();
  
  const int decimal_digits = 2;
  
  Kernel::FT scale_factor = 1.0;
  for (int digit = 0; digit < decimal_digits; ++digit) scale_factor *= 0.1;
  
  std::ofstream output_stream;
  output_stream.open(output_3dcm);
  output_stream << std::fixed;
  output_stream << std::setprecision(2);
  std::unordered_map<Kernel::Point_3, std::size_t> output_vertices;
  std::size_t num_polygons = 0;
  
  // Compute extent
  Kernel::FT x_min = map_polygons.front().x_min;
  Kernel::FT y_min = map_polygons.front().y_min;
  Kernel::FT z_min = 10000.0;
  for (std::vector<Polygon>::iterator current_polygon = map_polygons.begin(); current_polygon != map_polygons.end(); ++current_polygon) {
    for (auto const &current_vertex: current_polygon->triangulation.finite_vertex_handles()) {
      if (current_vertex->point().x() < x_min) x_min = current_vertex->point().x();
      if (current_vertex->point().y() < y_min) y_min = current_vertex->point().y();
      if (current_vertex->info().z < z_min) z_min = current_vertex->info().z;
    }
  }
  
  // Prepare CityJSON
  nlohmann::json cityjson;
  cityjson["type"] = "CityJSON";
  cityjson["version"] = "1.1";
  cityjson["transform"] = nlohmann::json::object();
  cityjson["transform"]["scale"] = {1.0, 1.0, 1.0};
  cityjson["transform"]["translate"] = {x_min, y_min, z_min};
  cityjson["CityObjects"] = nlohmann::json::object();
  cityjson["vertices"] = nlohmann::json::array();
  
  // Vertices
  for (std::vector<Polygon>::iterator current_polygon = map_polygons.begin(); current_polygon != map_polygons.end(); ++current_polygon) {
    for (Triangulation::Finite_faces_iterator current_face = current_polygon->triangulation.finite_faces_begin();
         current_face != current_polygon->triangulation.finite_faces_end();
         ++current_face) {
      if (current_face->info().interior == false) continue;
      for (int v = 0; v < 3; ++v) {
        Kernel::Point_3 point(current_face->vertex(v)->point().x(), current_face->vertex(v)->point().y(), current_face->vertex(v)->info().z);
        if (output_vertices.count(point) == 0) {
          cityjson["vertices"].push_back({int((point.x()-x_min)/scale_factor), int((point.y()-y_min)/scale_factor), int((point.z()-z_min)/scale_factor)});
          output_vertices[point] = output_vertices.size();
        }
      }
    } for (auto const &triangle: current_polygon->extra_triangles) {
      if (output_vertices.count(triangle.p1) == 0) {
        cityjson["vertices"].push_back({int((triangle.p1.x()-x_min)/scale_factor), int((triangle.p1.y()-y_min)/scale_factor), int((triangle.p1.z()-z_min)/scale_factor)});
        output_vertices[triangle.p1] = output_vertices.size();
      } if (output_vertices.count(triangle.p2) == 0) {
        cityjson["vertices"].push_back({int((triangle.p2.x()-x_min)/scale_factor), int((triangle.p2.y()-y_min)/scale_factor), int((triangle.p2.z()-z_min)/scale_factor)});
        output_vertices[triangle.p2] = output_vertices.size();
      } if (output_vertices.count(triangle.p3) == 0) {
        cityjson["vertices"].push_back({int((triangle.p3.x()-x_min)/scale_factor), int((triangle.p3.y()-y_min)/scale_factor), int((triangle.p3.z()-z_min)/scale_factor)});
        output_vertices[triangle.p3] = output_vertices.size();
      }
    }
  }
  
  // City objects
  for (std::vector<Polygon>::iterator current_polygon = map_polygons.begin(); current_polygon != map_polygons.end(); ++current_polygon) {
    
    cityjson["CityObjects"][current_polygon->id] = nlohmann::json::object();
    cityjson["CityObjects"][current_polygon->id]["type"] = current_polygon->cityjson_class;
    cityjson["CityObjects"][current_polygon->id]["geometry"] = nlohmann::json::array();
    cityjson["CityObjects"][current_polygon->id]["geometry"].push_back(nlohmann::json::object());
    cityjson["CityObjects"][current_polygon->id]["geometry"].back()["lod"] = "1.2";
    cityjson["CityObjects"][current_polygon->id]["geometry"].back()["boundaries"] = nlohmann::json::array();
    cityjson["CityObjects"][current_polygon->id]["geometry"].back()["type"] = "MultiSurface";
    cityjson["CityObjects"][current_polygon->id]["attributes"] = nlohmann::json::object();
    for (auto const &attribute: current_polygon->attributes) {
      cityjson["CityObjects"][current_polygon->id]["attributes"][attribute.first] = attribute.second;
    }
    
    // Triangles in polygon triangulation
    for (Triangulation::Finite_faces_iterator current_face = current_polygon->triangulation.finite_faces_begin();
         current_face != current_polygon->triangulation.finite_faces_end();
         ++current_face) {
      if (current_face->info().interior == false) continue;
      std::vector<std::size_t> vertices_of_face;
      for (int v = 0; v < 3; ++v) {
        vertices_of_face.push_back(output_vertices[Kernel::Point_3(current_face->vertex(v)->point().x(), current_face->vertex(v)->point().y(), current_face->vertex(v)->info().z)]);
      } cityjson["CityObjects"][current_polygon->id]["geometry"].back()["boundaries"].push_back({{vertices_of_face[0], vertices_of_face[1], vertices_of_face[2]}});
    }
    
    // Additional triangles
    for (auto const &triangle: current_polygon->extra_triangles) {
      std::size_t v1 = output_vertices[triangle.p1];
      std::size_t v2 = output_vertices[triangle.p2];
      std::size_t v3 = output_vertices[triangle.p3];
      cityjson["CityObjects"][current_polygon->id]["geometry"].back()["boundaries"].push_back({{v1, v2, v3}});
    }
    
  ++num_polygons;
  }
  
  output_stream << cityjson.dump() << std::endl;
  output_stream.close();
  std::cout << "Wrote 3D city model in ";
  print_timer(start_time);
  std::cout << " using ";
  print_memory_usage();
  return 0;
}

int write_terrain_obj(const char *output_terrain, Triangulation &terrain) {
  clock_t start_time = clock();
  std::ofstream output_stream;
  output_stream.open(output_terrain);
  output_stream << std::fixed;
  output_stream << std::setprecision(2);
  std::unordered_map<Kernel::Point_3, std::size_t> output_vertices;
  
  // Vertices
  for (auto const &vertex: terrain.finite_vertex_handles()) {
    Kernel::Point_3 point(vertex->point().x(), vertex->point().y(), vertex->info().z);
    output_stream << "v " << point.x() << " " << point.y() << " " << point.z() << "\n";
    output_vertices[point] = output_vertices.size()+1;
  }
  
  // Faces
  for (Triangulation::Finite_faces_iterator current_face = terrain.finite_faces_begin();
       current_face != terrain.finite_faces_end();
       ++current_face) {
    output_stream << "f";
    for (int v = 0; v < 3; ++v) {
      output_stream << " " << output_vertices[Kernel::Point_3(current_face->vertex(v)->point().x(),
                                                              current_face->vertex(v)->point().y(),
                                                              current_face->vertex(v)->info().z)];
    } output_stream << "\n";
  } output_stream.close();
  std::cout << "Wrote terrain in ";
  print_timer(start_time);
  std::cout << " using ";
  print_memory_usage();
  return 0;
}

int main(int argc, const char * argv[]) {
  
  const char *input_map = "/Users/ken/Library/Mobile Documents/com~apple~CloudDocs/Teaching/volta/data/vector/osmm.gpkg";
  const char *input_point_cloud = "/Users/ken/Library/Mobile Documents/com~apple~CloudDocs/Teaching/volta/data/point cloud/Exeter_VOLTAtest.laz";
  const char *output_terrain = "/Users/ken/Downloads/terrain.obj";
  const char *output_obj = "/Users/ken/Downloads/exeter.obj";
  const char *output_cityjson = "/Users/ken/Downloads/exeter.json";
  
  std::vector<Polygon> map_polygons;
  Edge_index edge_index;
  Point_cloud point_cloud;
  Point_index point_cloud_index;
  Triangulation terrain;
  
  std::cout.imbue(std::locale("en_US"));
  std::cout << std::fixed;
  std::cout << std::setprecision(2);
  
  load_map(input_map, map_polygons);
  triangulate_polygons(map_polygons);
  load_point_cloud(input_point_cloud, point_cloud);
  index_point_cloud(point_cloud, point_cloud_index);
  create_terrain_tin(map_polygons, point_cloud, point_cloud_index, terrain);
  write_terrain_obj(output_terrain, terrain);
  lift_flat_polygons(map_polygons, "Building", point_cloud, point_cloud_index, 0.7);
  lift_flat_polygons(map_polygons, "WaterBody", point_cloud, point_cloud_index, 0.1);
  lift_polygon_vertices(map_polygons, "Road", terrain);
  lift_polygon_vertices(map_polygons, "Railway", terrain);
  lift_polygon_vertices(map_polygons, "Bridge", terrain);
  lift_polygons(map_polygons, "PlantCover", terrain);
  lift_polygons(map_polygons, "LandUse", terrain);
  index_map(map_polygons, edge_index);
  create_vertical_walls(map_polygons, edge_index);
  write_3dcm_obj(output_obj, map_polygons);
  write_3dcm_cityjson(output_cityjson, map_polygons);
  
  return 0;
}

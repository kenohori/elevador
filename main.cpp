#include <iostream>
#include <mach/mach.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include "Enhanced_constrained_triangulation_2.h"

#include <ogrsf_frmts.h>

#include <pdal/PointTable.hpp>
#include <pdal/PointView.hpp>
#include <pdal/StageFactory.hpp>
#include <pdal/io/LasReader.hpp>
#include <pdal/io/LasHeader.hpp>
#include <pdal/Options.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Exact_predicates_tag Tag;

struct VertexInfo {
  Kernel::FT z;
};

struct FaceInfo {
  bool processed;
  bool interior;
  FaceInfo() {
    processed = false;
    interior = false;
  }
};

typedef CGAL::Triangulation_vertex_base_with_info_2<VertexInfo, Kernel> VertexBase;
typedef CGAL::Constrained_triangulation_face_base_2<Kernel> FaceBase;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo, Kernel, FaceBase> FaceBaseWithInfo;
typedef CGAL::Triangulation_data_structure_2<VertexBase, FaceBaseWithInfo> TriangulationDataStructure;
typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, TriangulationDataStructure, Tag> ConstrainedDelaunayTriangulation;
typedef Enhanced_constrained_triangulation_2<ConstrainedDelaunayTriangulation> Triangulation;

struct Ring {
  std::vector<Kernel::Point_2> vertices;
};

struct Polygon {
  Ring outer_ring;
  std::vector<Ring> inner_rings;
  Triangulation triangulation;
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

int main(int argc, const char * argv[]) {
  
  const char *input_map = "/Users/ken/Downloads/3dfier_os/osmm/osmm.gpkg";
//  const char *input_point_cloud = "/Users/ken/Downloads/3dfier_os/osmm/Exeter_VOLTAtest.laz";
//  const char *output_3dcm = "/Users/ken/Downloads/exeter.obj";
  
  std::vector<Polygon> map_polygons;
  
  // Prepare input map
  clock_t startTime = clock();
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
      
      if (wkbFlatten(input_feature->GetGeometryRef()->getGeometryType()) == wkbPolygon ||
          wkbFlatten(input_feature->GetGeometryRef()->getGeometryType()) == wkbTriangle) {
        OGRPolygon *input_polygon = input_feature->GetGeometryRef()->toPolygon();
        map_polygons.emplace_back();
        for (int current_vertex = 0; current_vertex < input_polygon->getExteriorRing()->getNumPoints(); ++current_vertex) {
          map_polygons.back().outer_ring.vertices.emplace_back(input_polygon->getExteriorRing()->getX(current_vertex),
                                                               input_polygon->getExteriorRing()->getY(current_vertex));
        } for (int current_inner_ring = 0; current_inner_ring < input_polygon->getNumInteriorRings(); ++current_inner_ring) {
          map_polygons.back().inner_rings.emplace_back();
          for (int current_vertex = 0; current_vertex < input_polygon->getInteriorRing(current_inner_ring)->getNumPoints(); ++current_vertex) {
            map_polygons.back().inner_rings.back().vertices.emplace_back(input_polygon->getInteriorRing(current_inner_ring)->getX(current_vertex),
                                                                         input_polygon->getInteriorRing(current_inner_ring)->getY(current_vertex));
          }
        }
      }
      
      else if (wkbFlatten(input_feature->GetGeometryRef()->getGeometryType()) == wkbMultiPolygon) {
        OGRMultiPolygon *input_multipolygon = input_feature->GetGeometryRef()->toMultiPolygon();
        for (int current_polygon = 0; current_polygon < input_multipolygon->getNumGeometries(); ++current_polygon) {
          OGRPolygon *input_polygon = input_multipolygon->getGeometryRef(current_polygon);
          map_polygons.emplace_back();
          for (int current_vertex = 0; current_vertex < input_polygon->getExteriorRing()->getNumPoints(); ++current_vertex) {
            map_polygons.back().outer_ring.vertices.emplace_back(input_polygon->getExteriorRing()->getX(current_vertex),
                                                                 input_polygon->getExteriorRing()->getY(current_vertex));
          } for (int current_inner_ring = 0; current_inner_ring < input_polygon->getNumInteriorRings(); ++current_inner_ring) {
            map_polygons.back().inner_rings.emplace_back();
            for (int current_vertex = 0; current_vertex < input_polygon->getInteriorRing(current_inner_ring)->getNumPoints(); ++current_vertex) {
              map_polygons.back().inner_rings.back().vertices.emplace_back(input_polygon->getInteriorRing(current_inner_ring)->getX(current_vertex),
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
  printTimer(startTime);
  std::cout << " using ";
  printMemoryUsage();
  
  // Basic polygon repair (pre-requisite for triangulation)
  auto current_polygon = map_polygons.begin();
  while (current_polygon != map_polygons.end()) {
    // Close polygons and rings
    if (current_polygon->outer_ring.vertices.back() != current_polygon->outer_ring.vertices.front()) {
      std::cout << "Warning: Last point != first. Adding it again at the end..." << std::endl;
      current_polygon->outer_ring.vertices.push_back(current_polygon->outer_ring.vertices.front());
    } for (auto &ring: current_polygon->inner_rings) {
      if (ring.vertices.back() != ring.vertices.front()) {
        std::cout << "Warning: Last point != first. Adding it again at the end..." << std::endl;
        ring.vertices.push_back(ring.vertices.front());
      }
    }
    
    // Delete degenerate polygons and rings
    if (current_polygon->outer_ring.vertices.size() < 4) {
      std::cout << "Deleting polygon with < 3 vertices..." << std::endl;
      auto polygon_to_erase = current_polygon;
      ++current_polygon;
      map_polygons.erase(polygon_to_erase);
      continue;
    } auto current_ring = current_polygon->inner_rings.begin();
    while (current_ring != current_polygon->inner_rings.end()) {
      if (current_ring->vertices.size() < 4) {
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
  startTime = clock();
  for (auto &polygon: map_polygons) {
    
    // Triangle (optimisation)
    if (polygon.outer_ring.vertices.size() == 4 && polygon.inner_rings.size() == 0) {
      if (polygon.outer_ring.vertices[0] != polygon.outer_ring.vertices[1] &&
          polygon.outer_ring.vertices[1] != polygon.outer_ring.vertices[2] &&
          polygon.outer_ring.vertices[2] != polygon.outer_ring.vertices[0]) {
        polygon.triangulation.insert(polygon.outer_ring.vertices[0]);
        polygon.triangulation.insert(polygon.outer_ring.vertices[1]);
        polygon.triangulation.insert(polygon.outer_ring.vertices[2]);
        for (Triangulation::Finite_faces_iterator current_face = polygon.triangulation.finite_faces_begin();
             current_face != polygon.triangulation.finite_faces_end();
             ++current_face) current_face->info().interior = true;
      } else {
        std::cout << "Warning: skipping degenerate triangle..." << std::endl;
      }
    }
    
    // Polygons
    
  } std::cout << "Triangulated polygons in ";
  printTimer(startTime);
  std::cout << " using ";
  printMemoryUsage();
  
  // Load point cloud
//  pdal::Options input_opts;
//  input_opts.add("filename", input_point_cloud);
//  pdal::PointTable input_table;
//  pdal::LasReader las_reader;
//  las_reader.setOptions(input_opts);
//  las_reader.prepare(input_table);
//  pdal::PointViewSet point_view_set = las_reader.execute(input_table);
//  pdal::PointViewPtr point_view = *point_view_set.begin();
//  pdal::Dimension::IdList dims = point_view->dims();
//  pdal::LasHeader las_header = las_reader.header();
  
  return 0;
}

#include <iostream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include "Enhanced_constrained_triangulation_2.h"

#include <ogrsf_frmts.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Exact_predicates_tag Tag;

struct FaceInfo {
  bool processed;
  bool interior;
  FaceInfo() {
    processed = false;
    interior = false;
  }
};

typedef CGAL::Triangulation_vertex_base_2<Kernel> VertexBase;
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

int main(int argc, const char * argv[]) {
  
  const char *input_map = "/Users/ken/Downloads/3dfier_os/osmm/osmm.gpkg";
  const char *input_point_cloud = "/Users/ken/Downloads/3dfier_os/osmm/Exeter_VOLTAtest.laz";
  
  std::vector<Polygon> polygons;
  
  // Prepare input map
  GDALAllRegister();
  GDALDataset *input_map_dataset = (GDALDataset*) GDALOpenEx(input_map, GDAL_OF_READONLY, NULL, NULL, NULL);
  if (input_map_dataset == NULL) {
    std::cerr << "Error: Could not open input file." << std::endl;
    return 1;
  }
  
  std::cout << "Input file: " << input_map << " type: " << input_map_dataset->GetDriverName() << std::endl;
  
  for (auto &&input_layer: input_map_dataset->GetLayers()) {
    input_layer->ResetReading();
    OGRFeature *input_feature;
    while ((input_feature = input_layer->GetNextFeature()) != NULL) {
      if (!input_feature->GetGeometryRef()) continue;
      
      if (wkbFlatten(input_feature->GetGeometryRef()->getGeometryType()) == wkbPolygon ||
          wkbFlatten(input_feature->GetGeometryRef()->getGeometryType()) == wkbTriangle) {
        OGRPolygon *input_polygon = input_feature->GetGeometryRef()->toPolygon();
        polygons.emplace_back();
        for (int current_vertex = 0; current_vertex < input_polygon->getExteriorRing()->getNumPoints(); ++current_vertex) {
          polygons.back().outer_ring.vertices.emplace_back(input_polygon->getExteriorRing()->getX(current_vertex),
                                                           input_polygon->getExteriorRing()->getY(current_vertex));
        } for (int current_inner_ring = 0; current_inner_ring < input_polygon->getNumInteriorRings(); ++current_inner_ring) {
          polygons.back().inner_rings.emplace_back();
          for (int current_vertex = 0; current_vertex < input_polygon->getInteriorRing(current_inner_ring)->getNumPoints(); ++current_vertex) {
            polygons.back().inner_rings.back().vertices.emplace_back(input_polygon->getInteriorRing(current_inner_ring)->getX(current_vertex),
                                                                     input_polygon->getInteriorRing(current_inner_ring)->getY(current_vertex));
          }
        }
      }
      
      else if (wkbFlatten(input_feature->GetGeometryRef()->getGeometryType()) == wkbMultiPolygon) {
        OGRMultiPolygon *input_multipolygon = input_feature->GetGeometryRef()->toMultiPolygon();
        for (int current_polygon = 0; current_polygon < input_multipolygon->getNumGeometries(); ++current_polygon) {
          OGRPolygon *input_polygon = input_multipolygon->getGeometryRef(current_polygon);
          polygons.emplace_back();
          for (int current_vertex = 0; current_vertex < input_polygon->getExteriorRing()->getNumPoints(); ++current_vertex) {
            polygons.back().outer_ring.vertices.emplace_back(input_polygon->getExteriorRing()->getX(current_vertex),
                                                             input_polygon->getExteriorRing()->getY(current_vertex));
          } for (int current_inner_ring = 0; current_inner_ring < input_polygon->getNumInteriorRings(); ++current_inner_ring) {
            polygons.back().inner_rings.emplace_back();
            for (int current_vertex = 0; current_vertex < input_polygon->getInteriorRing(current_inner_ring)->getNumPoints(); ++current_vertex) {
              polygons.back().inner_rings.back().vertices.emplace_back(input_polygon->getInteriorRing(current_inner_ring)->getX(current_vertex),
                                                                       input_polygon->getInteriorRing(current_inner_ring)->getY(current_vertex));
            }
          }
        }
      }
      
      else {
//        std::cout << "Ignoring non-areal type..." << std::endl;
      }
    }
    
    std::cout << polygons.size() << " polygons in map" << std::endl;
  }
  
  GDALClose(input_map_dataset);
  
  return 0;
}

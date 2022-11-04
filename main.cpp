// main
//
// Copyright © 2022,
// Ken Arroyo Ohori    ken@ken.mx
// All rights reserved.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#include "Elevador.h"

int main(int argc, const char * argv[]) {
  
  if (argc != 2) {
    std::cerr << "usage: elevador input_json_file" << std::endl;
    std::cerr << "for an example, see https://github.com/kenohori/elevador/blob/main/osmm.json" << std::endl;
    return 1;
  } std::ifstream input_stream(argv[1]);
  nlohmann::json input_params = nlohmann::json::parse(input_stream);
  
  std::string class_attribute = input_params["class_attribute"].get<std::string>();
  
  Elevador elevador;

  std::cout.imbue(std::locale("en_US"));
  std::cout << std::fixed;
  std::cout << std::setprecision(2);

  if (input_params.count("input_map") == 1) {
    std::string input_map = input_params["input_map"].get<std::string>();
    elevador.load_map(input_map, class_attribute);
    elevador.triangulate_polygons();
  } if (input_params.count("input_point_cloud") == 1) {
    std::string input_point_cloud = input_params["input_point_cloud"].get<std::string>();
    elevador.load_point_cloud(input_point_cloud);
    elevador.index_point_cloud();
  }
  
  std::set<std::string> ignore_classes;
  if (input_params.count("terrain_ignore_classes") == 1) {
    for (auto const &ignore_class: input_params["terrain_ignore_classes"]) {
      ignore_classes.insert(ignore_class.get<std::string>());
    }
  } elevador.create_terrain_tin(ignore_classes);
  
  if (input_params.count("output_terrain") == 1) {
    std::string output_terrain = input_params["output_terrain"].get<std::string>();
    elevador.write_terrain_obj(output_terrain);
  }

  for (auto &operation: input_params["operations"]) {
    std::string operation_type = operation["operation"].get<std::string>();
    if (operation_type == "lift_flat_polygons") {
      elevador.lift_flat_polygons(operation[class_attribute].get<std::string>().c_str(), operation["ratio"].get<double>());
    } else if (operation_type == "lift_polygon_vertices") {
      elevador.lift_polygon_vertices(operation[class_attribute].get<std::string>().c_str());
    } else if (operation_type == "lift_polygons") {
      elevador.lift_polygons(operation[class_attribute].get<std::string>().c_str());
    } else {
      std::cerr << "Unknown operation: " << operation_type << std::endl;
    }
  }
  
//  compute_height_stats(map, point_cloud, point_cloud_index, terrain);

  elevador.index_map();
  elevador.create_vertical_walls();
  
  if (input_params.count("output_obj") == 1) {
    std::string output_obj = input_params["output_obj"].get<std::string>();
    elevador.write_3dcm_obj(output_obj);
  } if (input_params.count("output_cityjson") == 1) {
    std::string output_cityjson = input_params["output_cityjson"].get<std::string>();
    elevador.write_3dcm_cityjson(output_cityjson);
  }
  
  return 0;
}

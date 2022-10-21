This is the code of a prototype that creates 3D city models from a combination of a 2D map and a point cloud.
It was made for my secondment at the [Ordnance Survey](https://www.ordnancesurvey.co.uk), Southampton, UK.
So far, it has only been tested with sample data around Exeter, UK, which was supplied by OS.
However, you can get similar data from England's Environmental Agency [here](https://www.arcgis.com/apps/webappviewer/index.html?id=f765c2a97d644f08927d5cd5abe58d87) and there's a free sample of the [OS MasterMap](https://www.ordnancesurvey.co.uk/business-government/products/mastermap-topography), which you can reclassify with [this code](https://github.com/kenohori/osclassifier).

In order to run this, you need C++20 or higher, as well as the following libraries: [GDAL](https://gdal.org) (vector input), [CGAL](https://www.cgal.org) (geometric operations), Niels Lohmann's [JSON for Modern C++](https://github.com/nlohmann/json) (JSON input/output) and [PDAL](https://pdal.io/) (point cloud input).
A [CMake](https://cmake.org) file is provided for compilation.

It achieves similar results to [3dfier](https://github.com/tudelft3d/3dfier) and [City4CFD](https://github.com/tudelft3d/city4cfd) but the ultimate aim is to have a more general approach that works in more countries.

More information to come.

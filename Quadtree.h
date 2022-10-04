#ifndef Quadtree_h
#define Quadtree_h

template <class Kernel, class Point_cloud>
struct Quadtree {
  typename Kernel::FT x_min, x_max, y_min, y_max;
  unsigned int depth;
  std::vector<typename Point_cloud::Index> points;
  Quadtree *upper_left, *upper_right, *lower_left, *lower_right;
  
  Quadtree() {
    depth = 0;
    upper_left = NULL;
    upper_right = NULL;
    lower_left = NULL;
    lower_right = NULL;
  }
  
  Quadtree(unsigned int depth, typename Kernel::FT x_min, typename Kernel::FT x_max, typename Kernel::FT y_min, typename Kernel::FT y_max) {
    this->depth = depth;
    upper_left = NULL;
    upper_right = NULL;
    lower_left = NULL;
    lower_right = NULL;
    this->x_min = x_min;
    this->x_max = x_max;
    this->y_min = y_min;
    this->y_max = y_max;
  }
  
  void compute_extent(Point_cloud &point_cloud) {
    typename Kernel::FT x_min = point_cloud.point(*point_cloud.begin()).x();
    typename Kernel::FT x_max = point_cloud.point(*point_cloud.begin()).x();
    typename Kernel::FT y_min = point_cloud.point(*point_cloud.begin()).y();
    typename Kernel::FT y_max = point_cloud.point(*point_cloud.begin()).y();
    for (auto const &point: point_cloud.points()) {
      if (point.x() < x_min) x_min = point.x();
      if (point.x() > x_max) x_max = point.x();
      if (point.y() < y_min) y_min = point.y();
      if (point.y() > y_max) y_max = point.y();
    } this->x_min = x_min;
    this->x_max = x_max;
    this->y_min = y_min;
    this->y_max = y_max;
  }
  
  void insert_point(Point_cloud &point_cloud, typename Point_cloud::Index point) {
      points.push_back(point);
  }
  
  void optimise(Point_cloud &point_cloud, std::size_t bucket_size, unsigned int maximum_depth) {
    if (points.size() < bucket_size || depth+2 >= maximum_depth) return;
    typename Kernel::FT x_split = (x_min+x_max)/2.0;
    typename Kernel::FT y_split = (y_min+y_max)/2.0;
    upper_left = new Quadtree(depth+1, x_min, x_split, y_split, y_max);
    upper_right = new Quadtree(depth+1, x_split, x_max, y_split, y_max);
    lower_left = new Quadtree(depth+1, x_min, x_split, y_min, y_split);
    lower_right = new Quadtree(depth+1, x_split, x_max, y_min, y_split);
    
    for (auto point: points) {
      if (point_cloud.point(point).x() < x_split) {
        if (point_cloud.point(point).y() < y_split) lower_left->insert_point(point_cloud, point);
        else upper_left->insert_point(point_cloud, point);
      } else {
        if (point_cloud.point(point).y() < y_split) lower_right->insert_point(point_cloud, point);
        else upper_right->insert_point(point_cloud, point);
      }
    } points.clear();
    upper_left->optimise(point_cloud, bucket_size, maximum_depth);
    upper_right->optimise(point_cloud, bucket_size, maximum_depth);
    lower_left->optimise(point_cloud, bucket_size, maximum_depth);
    lower_right->optimise(point_cloud, bucket_size, maximum_depth);
  }
  
  void print_info() {
    std::cout << "Quadtree extent: X = [" << x_min << ", " << x_max << "] Y = [" << y_min << ", " << y_max << "]" << std::endl;
    std::vector<std::size_t> nodes_by_depth;
    Quadtree *biggest_node = this;
    std::list<Quadtree *> to_visit;
    to_visit.push_back(this);
    while (!to_visit.empty()) {
      if (to_visit.front()->points.size() > biggest_node->points.size()) biggest_node = to_visit.front();
      if (nodes_by_depth.size() < to_visit.front()->depth+1) nodes_by_depth.push_back(1);
      else nodes_by_depth[to_visit.front()->depth]++;
      if (to_visit.front()->upper_left != NULL) to_visit.push_back(to_visit.front()->upper_left);
      if (to_visit.front()->upper_right != NULL) to_visit.push_back(to_visit.front()->upper_right);
      if (to_visit.front()->lower_left != NULL) to_visit.push_back(to_visit.front()->upper_left);
      if (to_visit.front()->lower_right != NULL) to_visit.push_back(to_visit.front()->lower_right);
      to_visit.pop_front();
    } std::cout << "Quadtree nodes by depth:";
    for (auto const &n_nodes: nodes_by_depth) std::cout << " " << n_nodes;
    std::cout << std::endl;
    std::cout << "Biggest node: " << biggest_node->points.size() << " points at depth " << biggest_node->depth << " for X = [" << biggest_node->x_min << ", " << biggest_node->x_max << "] Y = [" << biggest_node->y_min << ", " << biggest_node->y_max << "]" << std::endl;
  }
  
  ~Quadtree() {
    if (upper_left != NULL) delete upper_left;
    if (upper_right != NULL) delete upper_right;
    if (lower_left != NULL) delete lower_left;
    if (lower_right != NULL) delete lower_right;
  }
};

#endif /* Quadtree_h */

/*
 * Adaptive Mesh Class
 *
 * This class will define the adaptive mesh DE solver used to solve the Stars
 * PDEs. An implementation of RK4 will also be included. 
 *
 * This will provide the numerical component of the Stars project
 */

#ifndef _AdaptMesh_
#define _AdaptMesh_

#include <cmath>

using namespace std;

class AdaptMesh{

 public:
  AdaptMesh();
  ~AdaptMesh();
  void init();

 private:
};

#endif

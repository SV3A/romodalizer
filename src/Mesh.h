#ifndef MESH_H
#define MESH_H

#include <iostream>
#include "Eigen/Dense"

class Mesh
{
  private:
    unsigned int numEl;
    void setGeometry(Eigen::MatrixXd & mesh);

  public:

    // 5xm matrix containing all info related to each element (column)
    Eigen::Matrix<double, 5, Eigen::Dynamic> elements;

    Mesh(Eigen::MatrixXd &mesh);

    void setDensity(double rho, unsigned int element = 0);
    void setEmod(double eMod, unsigned int element = 0);
};

#endif

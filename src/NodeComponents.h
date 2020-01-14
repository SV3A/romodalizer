#ifndef NODECOMPONENTS_H
#define NODECOMPONENTS_H

#include <iostream>
#include <array>
#include "Eigen/Dense"

class Disc
{
  public:
    float m;
    float iD;
    float iP;

    Disc(float mass, float momentInert, float momentPolar);
};


class Bearing
{
  public:
    Eigen::Matrix4d localK;

    Bearing(Eigen::Matrix4d& k);
};

#endif

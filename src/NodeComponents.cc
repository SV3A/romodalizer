#include "NodeComponents.h"

Disc::Disc(float mass, float momentInert, float momentPolar)
{
  m  = mass;
  iD = momentInert;
  iP = momentPolar;
};

Bearing::Bearing(Eigen::Matrix4d& k)
{
  localK = k;
};


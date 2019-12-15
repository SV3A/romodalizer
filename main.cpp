#include <iostream>
#include <cmath>
#include "ModalAnalysis.h"
#include "Eigen/Dense"

class RotorSystem
{
  private:
    // The angular velocity of the shaft
    float omega;

    // Struct containing the shaft properties
    struct Shaft{
      float d;   // Diameter [m]
      float E;   // E-module [N/m^2]
      float rho; // Density  [kg/m^3]
    };

  public:
    // Internal shaft struct
    Shaft shaft;

    RotorSystem(){
      // Setup shaft
      shaft.d   = 0.025;
      shaft.E   = 69e9;
      shaft.rho = 2770;
    }

    void setOmega(float omega)
    {
      this->omega = omega*2.0*M_PI;
    }

    float getOmega(std::string unit){
      // Return the angular velocity 
      if (unit == "rad")
        return omega;
      else if (unit == "rpm")
        return omega*60.0/2.0/M_PI;
      else if (unit == "hz")
        return omega/(2.0 * M_PI);
      else
        // TODO: Make proper exception
        std::cout << "\n\nWrong unit in getOmega" << std::endl;
        exit(1);
    }
};

int main(int argc, char *argv[])
{
  RotorSystem *rotSys = new RotorSystem();

  rotSys->setOmega(1.0);

  // Define mesh
  double sd = rotSys->shaft.d;
  Eigen::MatrixXf mesh(3, 5); mesh << 162,  88, 120,  52,  63,
                                       sd,  sd,  sd,  sd,  sd,
                                        3,   2,   3,   1,   2;

  ModalAnalysis *modAnalysis = new ModalAnalysis(mesh);

  modAnalysis->printInfo();

  // Cleanup
  delete rotSys;
  delete modAnalysis;
}

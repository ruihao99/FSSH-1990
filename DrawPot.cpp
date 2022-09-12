#include <iostream>
#include <fstream>
#include <string>

#include <math.h>
#include <array>
#include <vector>
#include <numeric> 

#include <Eigen/Dense>

#include "./include/utils.hpp"
#include "./include/pot2d.hpp"

template<class T>
void evaluate_pot(T &pot, std::string fname){
    int n1 = 10;
    int n2 = 15;
    std::ofstream fV(fname+"-V.dat"); 
    std::vector<double> X = linspace(-10., 10., 1000);
    fV << std::setw(n1) << "# x" << std::setw(n2) << "E1" << std::setw(n2) << "E2" 
       << std::setw(n2) << "dE1" << std::setw(n2) << "dE2" << std::setw(n2) << "d12" 
       << std::setw(n2) << "d21" << std::endl;


    //std::ofstream fd(fname+"-dij.dat");
    //fd << std::setw(n1) << "# x" << std::setw(n2) << "d11" << std::setw(n2) << "d12" 
    //   << std::setw(n2) << "d21" << std::setw(n2) << "d22" << std::endl;

    for (auto iter = X.begin(); iter<X.end(); iter++)
    {
        pot.eval_V(*iter);    // evaluate V for given x
        pot.eval_dVdx(*iter); // evaluate dVdx for given x
        //pot.diag();
        //pot.eval_d(*iter, h);
        pot.eval_d(*iter);
        fV << std::setw(n1) << *iter << std::setw(n2) << pot.get_Ei(0) 
           << std::setw(n2) << pot.get_Ei(1) << std::setw(n2) << pot.get_Fi(0) 
           << std::setw(n2) << pot.get_Fi(1) << std::setw(n2) << pot.get_dij(0, 1) 
           << std::setw(n2) << pot.get_dij(1, 0) << std::setw(n2) << pot.get_phiij(0, 0) 
           << std::setw(n2) << pot.get_phiij(1, 0) << std::setw(n2) << pot.get_phiij(1, 0)
           << std::setw(n2) << pot.get_phiij(1, 1) << std::endl;

        //fd << std::setw(n1) << *iter << std::setw(n2) << pot.get_dij(0, 0) << std::setw(n2)
        //   << pot.get_dij(0, 1) << std::setw(n2) << pot.get_dij(1, 0) << std::setw(n2) 
        //   << pot.get_dij(1, 1) << std::endl;
    }
    fV.close();
    //fd.close();
}


int main()
{
    // *** POTENTIAL 1: SINGLE AVOIDED CROSSING *** //
    std::array<double, 4> C1 = {0.01, 1.6, 0.005, 1.0};
    P_SAC p_sac(C1);
    evaluate_pot(p_sac, "./SAC");

    // *** POTENTIAL 2: DUAL AVOIDED CROSSING *** //
    std::array<double, 5> C2 = {0.1, 0.28, 0.015, 0.06, 0.05};
    P_DAC p_dac(C2);
    evaluate_pot(p_dac, "./DAC");
   
    // *** POTENTIAL 3: EXTENDED COUPLING WITH REFLECTION *** //
    std::array<double, 3> C3 = {6E-4, 0.10, 0.90};
    P_ECR p_ecr(C3);
    evaluate_pot(p_ecr, "./ECR");
}

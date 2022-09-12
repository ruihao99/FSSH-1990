/*
 * Simulation files for FSSH
 * Author: Rui-Hao Bi
 * email: biruihao@westlake.edu.cn
 */
#include <string>
#include <vector>
#include "./include/integrator.hpp"

// Run FSSH molecular dynamics
// Only one instances
template<class P>
void run_FSSH_MD(P &p2d, std::string fname)
{
    const int    n1     = 15;
    int          nsteps = 20000;
    const double qi     = -10;     // initial positions 
    //const double pi     = 10;    // initial momentum
    double dt           = 10;            // guess dt 
    int           index = 0;

    int N_traj   = 2000;
    //std::vector<double> Px1 = linspace(0., 4., 5);
    //std::vector<double> Px2 = linspace(4.1, 15., 44);
    //std::vector<double> Px3 = linspace(15.5, 35., 20);
    //std::vector<double> Px = Px1 + Px2 + Px3;
    //
    //std::vector<double> Px1 = linspace(0.5, 4.5, 5);
    //std::vector<double> Px2 = linspace(5., 20., 25);
    //std::vector<double> Px3 = linspace(22., 100., 40);
    //std::vector<double> Px = Px1 + Px2 + Px3; 

    std::vector<double> Px1 = linspace(0.5, 2.8, 5);
    std::vector<double> Px2 = linspace(3., 23., 50);
    std::vector<double> Px3 = linspace(23.5, 33.5, 35);
    std::vector<double> Px = Px1 + Px2 + Px3; 

    for (auto iter = Px.begin(); iter<Px.end(); iter++)
    {
        index = iter - Px.begin();
        // Define file name
        std::ofstream f(fname + "-" + std::to_string(index)  + ".dat");
        f << "# " << *iter << std::endl;
        // Define integrators

        f << std::setw(n1) << "# TrajNum" << std::setw(n1) << "final k" 
          << std::setw(n1) << "final q" << std::setw(n1) << "final p"
          << std::endl;
        if (*iter <= 0)
        {
            dt = 1;
            nsteps = 10000;
        }
        else if (dt > 30){
            dt = 1;
            nsteps = 20. / *iter * MyConst::mass / dt;
        }
        else
        {
            nsteps = 10000;
            dt = 20. / *iter * MyConst::mass / nsteps;
        }
        IntegratorFSSH<double, P> int_fssh{qi, *iter, dt};

        for (int j=0; j<N_traj; j++)
        {
            for (int i=0; i<nsteps; i++)
            {
                p2d.eval_V(int_fssh.qp(0));
                p2d.eval_dVdx(int_fssh.qp(0));
                p2d.eval_d(int_fssh.qp(0));

                int_fssh.evolve_nuc(p2d);
                int_fssh.evolve_wc(p2d);
                int_fssh.eval_a(p2d);
            }
            f << std::setw(n1) << j << std::setw(n1) << int_fssh.k 
              << std::setw(n1) << int_fssh.qp(0) << std::setw(n1) 
              << int_fssh.qp(1) << std::endl;
            int_fssh.reset();
        }
        f.close();
    }
}

int main()
{
    // *** POTENTIAL 1: SINGLE AVOIDED CROSSING *** //
    std::array<double, 4> C1 = {0.01, 1.6, 0.005, 1.0};
    P_SAC p_sac(C1);

    // *** POTENTIAL 2: DUAL AVOIDED CROSSING *** //
    std::array<double, 5> C2 = {0.1, 0.28, 0.015, 0.06, 0.05};
    P_DAC p_dac(C2);

    // *** POTENTIAL 3: EXTENDED COUPLING WITH REFLECTION *** //
    std::array<double, 3> C3 = {6E-4, 0.10, 0.90};
    P_ECR p_ecr(C3);

    //run_FSSH_MD(p_sac, "./SAC/SAC");
    //run_FSSH_MD(p_dac, "./DAC/DAC");
    run_FSSH_MD(p_ecr, "./ECR/ECR");

    return 0;
}

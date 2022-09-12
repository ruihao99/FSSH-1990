/*
 * Head file for FSSH integrators
 * Author: Rui-Hao Bi
 * email: biruihao@westlake.edu.cn
 */
#include <math.h>
#include <array>
#include <vector>
#include <complex>
#include <random>
#include <chrono>
#include <fstream>

#include <boost/numeric/odeint.hpp>
namespace odeint = boost::numeric::odeint;
using namespace std::complex_literals;

#include <Eigen/Dense>
#include "pot2d.hpp"
#include "utils.hpp"


// Class for nuclei motion <EOM>
template <class X, class Q>
struct NucleiMotion
{
    X Fk = 0.;
    void update_F(X &F){Fk = F;}
    void operator()(const Q &x, Q &dxdt, double t)
    {
        dxdt(0) = x(1) / MyConst::mass;
        dxdt(1) = - Fk;
    }
};

// Class for wavefunction coefficient time-evolution <enq 8>
template <class X, class C, class P>
struct WaveCoeff
{
    typedef Eigen::Matrix<X, 2, 2> Matrix2X; 
    typedef Eigen::Vector<X, 2> Vector2X; 
    Matrix2X d {{0, 0}, {0, 0}};
    Vector2X E {0, 0};
    X p = 0;

    void update_V_and_d(P &p2d){E = p2d.E; d = p2d.d;}
    void update_p(X &p_in){p = p_in;}

    void operator()(const C &x, C &dxdt, double t)
    {
        dxdt(0) = x(0) * E(0) / 1i - x(1) * p/MyConst::mass * d(0, 1);
        dxdt(1) = x(1) * E(1) / 1i - x(0) * p/MyConst::mass * d(1, 0);
    }
};

// Class for density matrix time evoluition

// Base template for surface hopping integrator
// X is the accuracy of this calculation
// P is the potential class 
// User customize alogorithm for hopping
template <class X, class P>
class IntegratorFSSH
{
    typedef Eigen::Vector<std::complex<X>, 2> Vector2cX;
    typedef Eigen::Vector<X, 2> Vector2X;
    typedef Eigen::Matrix<std::complex<X>, 2, 2> Matrix2cX; 

    public:
        IntegratorFSSH(const X & q, const X & p, const X & dt_in)
        {
            qp << q, p;
            qp0 << q, p;
            dt = dt_in;
        };
        NucleiMotion<X, Vector2X> nm;
        WaveCoeff<X, Vector2cX, P> wc;
        odeint::runge_kutta4<Vector2X>  rk4_qp;
        //odeint::velocity_verlet<Vector2X>  vv_qp;
        odeint::runge_kutta4<Vector2X>  rk4_a;
        odeint::runge_kutta4<Vector2cX> rk4_c;

        void eval_a(P &p2d)
        {

            A = c * c.adjoint();
            //a = A.diagonal().real();

            b(0) = -2. * (std::conj(A(1, 0)) * qp(1)/MyConst::mass * p2d.get_dij(1, 0)).real();
            b(1) = -2. * (std::conj(A(0, 1)) * qp(1)/MyConst::mass * p2d.get_dij(0, 1)).real();
        }
        
        void evolve_nuc(P &p2d)
        {
            this->update_zeta();
            this->eval_g();
            if (g(k) > zeta)
            {
                PEi = p2d.get_Ei(k);
                this->switch_state();       // Make switch
                //std::cout << "swith occur. Now g_kj: " << g(k) << ", and zeta: " << zeta << "." << std::endl;
                DeltaPE = p2d.get_Ei(k) - PEi;
                this->vel_adjust();
                //std::cout << KE << std::endl;
                if (not vel_pass)        // If the current switch doesn't pass the velocity  
                {
                    this->switch_state();    // adjustment test, swith fail, therefore switch back. 
                    //std::cout << "swith rejected. With KE: " << KE << ", and PE change: " << DeltaPE << "." << std::endl;
                }
                else
                    NS += 1;
            }
                            
            nm.update_F(p2d.get_Fi(k));
            rk4_qp.do_step(nm, qp, tn, dt);
            //vv_qp.do_step(nm, qp, tn, dt);
            tn += dt;
        }

        void evolve_wc(P &p2d)
        {
            wc.update_p(qp(1));
            wc.update_V_and_d(p2d);
            rk4_c.do_step(wc, c, tc, dt);
            tc += dt;
        }
        
        void vel_adjust()
        {
            KE  = 1./2 * pow(qp(1), 2) / MyConst::mass;
            if (DeltaPE > KE){vel_pass  = false;}
            else
            {
                KEf = -DeltaPE + KE; 
                pf = sqrt(2.*MyConst::mass * KEf); 
                if (qp(1) > 0){qp(1) = pf;}
                else {qp(1) = -pf;}

                vel_pass = true;
            }
        }

        void update_zeta(){zeta = distro(gen);}
        void eval_g(){g(k) = dt * b(k) / A(k, k).real();}
        void switch_state(){k = (int)(0.5 - ((double)k - 0.5));}
        void conclusion(){
            std::cout << "# Total timesteps: " << tn << std::endl;
            std::cout << "# Number of successful switchs: " << NS << std::endl;
        }
        void reset(){
            qp = qp0;
            tn = tc = k = NS = 0;
            c << 1, 0;    
            A << 1., 0, 
                 0., 0.;  
            //a << 1, 0;    
            b << 0, 0;    
            g << 0, 0;    
        }
        void conclusion(std::ofstream &f){
            f << "# Total timesteps: " << tn << std::endl;
            f << "# Number of successful switchs: " << NS << std::endl;
        }

        X tn = 0.0;                 // Current time for the NucleiMotion 
        X tc = 0.0;                 // Current time for the evolution of WaveFuncCoeff
        X dt = 1.0;                 // Integration time interval
        int k = 0;                  // The system evolves on the PES of state #k.
        int NS = 0;                 // # of switches for current time
        X zeta = 1.;                 // zeta, the random number between 0 and 1
        X PEi, DeltaPE, KE, KEf, pf = 0.;         // Tmp values used for velocity adjustment step
        bool vel_pass = true;         // Tmp values used for velocity adjustment step
        Vector2X  qp, qp0;              // pos and momentum
        Vector2cX c  {1, 0};        // Complex expansion coefficients.
        Matrix2cX A  {{1., 0.}, 
                      {0., 0.}};    // Density matrix coefficients.
        ////Vector2X  a  {1, 0};        // The diagnal of density matrix.
        Vector2X  b  {0, 0};        // The off-diagnal of B, i.e. b_kl defined in eqn14.
        Vector2X  g  {0, 0};        // Switching criterion g_kj.

    private:
        unsigned seed_chrono = std::chrono::system_clock::now().time_since_epoch().count();
        std::default_random_engine gen {seed_chrono};
        //std::default_random_engine gen;
        std::uniform_real_distribution<X> distro{0.0, 1.0};

};

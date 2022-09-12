/*
 * Head file for 2-state problem potential
 * Author: Rui-Hao Bi
 * email: biruihao@westlake.edu.cn
 */
#include <fstream>
#include <string>

#include <math.h>
#include <array>
#include <vector>
#include <numeric> 

#include <Eigen/Dense>


template <class X, class C>
class Potential2D
{
    typedef Eigen::Matrix<X, 2, 2> Matrix2X;
    typedef Eigen::Vector<X, 2> Vector2X;
    public:
        C coeff = 0;
        Potential2D(C coeff_in) : coeff(coeff_in){}; 
        // coeff_in is the input coefficient for the potential 
        
        // 'eval_V' and 'eval_dVdx' defined the
        // 2D-potentials functions, namely V11, V12, V21, and V22,
        // and the derivatives dVdx11, dVdx12, dVdx21, dVdx22.

        virtual void eval_V(const X &x){}
        virtual void eval_dVdx(const X &x){}
        virtual void eval_d(const X &x, const X &h)
        {
            // Numerical method for the evaluation of d
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigensolver(V);
            if (eigensolver.info() != Eigen::Success) abort();
            E = eigensolver.eigenvalues();
            phi = eigensolver.eigenvectors();

            // The two diagonal terms are 0 by difinition
            V0 = V;
            (*this).eval_V(x+h);
            Vp = V;
            (*this).eval_V(x-h);
            Vn = V;
            V = V0;

            Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> esp(Vp);
            if (esp.info() != Eigen::Success) abort();
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> esn(Vn);
            if (esn.info() != Eigen::Success) abort();
            dphi = (esp.eigenvectors() - esn.eigenvectors()) / (2 * h);
            d = phi.transpose() * dphi;
        }

        virtual void eval_d(const X &x)
        {
            // "Exact"-method for evaluating d
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigensolver(V);
            if (eigensolver.info() != Eigen::Success) abort();
            E = eigensolver.eigenvalues();
            phi = eigensolver.eigenvectors();

            if (this->get_Vij(1, 1) > this->get_Vij(0, 0))
                sign = 1;
            else
                sign = -1;

            F = phi.adjoint() * dVdx * phi;
            d(0, 1) = sign * F(0, 1) / (E(1) - E(0));
            d(1, 0) = sign * F(1, 0) / (E(0) - E(1));

            //d(0, 1)  = sign * phi.col(0).adjoint().dot(dVdx * phi.col(1)) / (E(1) - E(0));
            //d(1, 0)  = sign * phi.col(1).adjoint().dot(dVdx * phi.col(0)) / (E(0) - E(1)) ;
            //d(0, 1)   = (phi.col(0).transpose() * dVdx * phi.col(1))(0) / (E(1) - E(0));
            //d(1, 0)   = (phi.col(1).transpose() * dVdx * phi.col(0))(0) / (E(0) - E(1)) ;

            //dE(0) = phi.col(0).transpose().dot((dVdx * phi.col(0))) ;
            //dE(1) = phi.col(1).transpose().dot((dVdx * phi.col(1))) ;
        }

        // For public access of protected data
        double &get_Vij(int i, int j){return V(i, j);}
        double &get_dVdxij(int i, int j){return dVdx(i, j);}
        double &get_Ei(int i){return E(i);}
        double &get_Fi(int i){return F(i, i);}
        double &get_phiij(int i, int j){return phi(i, j);}
        double &get_dij(int i, int j){return d(i, j);}

        Matrix2X V, Vp, Vn, V0 = Matrix2X::Zero();      // Vp and Vn <-> V(x+h) and V(x-h). V0 is tmp.
        Matrix2X dVdx = Matrix2X::Zero();
        Matrix2X d = Matrix2X::Zero(); 
        Matrix2X phi = Matrix2X::Zero();
        Vector2X E = Vector2X::Zero();
    protected:
        /*
         * Here stores the data relating to a 2D-Pot
         * + 2x2 matrix V    is the Hamiltonian matrix
         * + 2x2 matrix dVdx is the matrix of the potential gradient  
         * + 2-vector   E    contains the eigen value for V
         * + 2x2 matrix phi  contains the two eigen vectors for V
         * + 2x2 matrix d    is the coulpling constant matrix
         */
        //Vector2X dE = Vector2X::Zero();
        Matrix2X F = Matrix2X::Zero();
        Matrix2X dphi = Matrix2X::Zero();
        int sign = 1;
};

/*
 * This is a demo for 2-states Single Avoided Crossing (SAC) potential.
 */
class P_SAC : public Potential2D<double, std::array<double, 4>>
{
    public:
        P_SAC(std::array<double, 4> coeff):Potential2D(coeff){}

        void eval_V(const double &x)
        {
            V11 = coeff[0] * (1 - exp(-abs(coeff[1]*x)));
            V12 = coeff[2] * exp(- coeff[3] * pow(x, 2));

            if (x >= 0){
                V(0, 0) =  V11;
                V(1, 1) = -V11;
            }
            else{
                V(0, 0) = -V11;
                V(1, 1) =  V11;
            }
            V(0, 1) = V12;
            V(1, 0) = V12;
        }

        void eval_dVdx(const double &x)
        {
            dVdx11 = coeff[0] * coeff[1] * exp(-abs(coeff[1]*x));
            dVdx12 = -2. * coeff[2] * coeff[3] * x * exp(-coeff[3] * pow(x, 2));

            dVdx(0, 0) =  dVdx11;
            dVdx(0, 1) =  dVdx12;
            dVdx(1, 0) =  dVdx12;
            dVdx(1, 1) = -dVdx11;
        }
        double V11, V12, dVdx11, dVdx12 = 0; 
};

class P_DAC : public Potential2D<double, std::array<double, 5>>
{
    public:
        P_DAC(std::array<double, 5> coeff):Potential2D(coeff){}

        void eval_V(const double &x)
        {
            V12 = coeff[2] * exp(-coeff[3] * pow(x, 2));
            V22 = -coeff[0] * exp(- coeff[1] * pow(x, 2)) + coeff[4];

            V(0, 0) = V11;
            V(1, 1) = V22;
            V(0, 1) = V12;
            V(1, 0) = V12;
        }

        void eval_dVdx(const double &x)
        {
            dVdx22 =  2. * coeff[0] * coeff[1] * x * exp(-coeff[1] * pow(x, 2));
            dVdx12 = -2. * coeff[2] * coeff[3] * x * exp(-coeff[3] * pow(x, 2));

            dVdx(0, 0) = dVdx11;
            dVdx(0, 1) = dVdx12;
            dVdx(1, 0) = dVdx12;
            dVdx(1, 1) = dVdx22;
        }
        double V11, V12, V22, dVdx11, dVdx12, dVdx22 = 0.;
};

class P_ECR : public Potential2D<double, std::array<double, 3>>
{
    public:
        P_ECR(std::array<double, 3> coeff):Potential2D(coeff){}

        void eval_V(const double &x)
        {
            U12 = coeff[1] * exp(- abs(coeff[2] * x));

            if (x >= 0){
                V(0, 1) = coeff[1] * 2. - U12;
                V(1, 0) = coeff[1] * 2. - U12;
            }
            else{
                V(0, 1) = U12;
                V(1, 0) = U12;
            }
            V(0, 0) =  -V11;
            V(1, 1) =  V11;
        }

        void eval_dVdx(const double &x)
        {
            dVdx12 = coeff[1] * coeff[2] * exp(-abs(coeff[2] * x));
            dVdx(0, 0) = dVdx11;
            dVdx(0, 1) = dVdx12;
            dVdx(1, 0) = dVdx12;
            dVdx(1, 1) = dVdx11;
        }

        double V11 = coeff[0];
        double dVdx11, U12, dVdx12 = 0;
};

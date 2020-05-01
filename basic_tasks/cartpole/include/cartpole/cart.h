#ifndef CART_H
#define CART_H

#define N 10
#define d 10.0
#define pi 3.14159265
#define l 0.5
#define m1 0.5
#define m2 0.5
#define g -9.81
#define T 5.0
#define h T/(N-1.0)

#include <iostream>
#include "ifopt/variable_set.h"
#include "ifopt/constraint_set.h"
#include "ifopt/problem.h"
#include "ifopt/solver.h"
#include "ifopt/composite.h"
#include "ifopt/bounds.h"
#include "ifopt/cost_term.h"
#include "Eigen/Eigen"
#include <iostream>
#include <math.h>


namespace cartpole{
    using namespace ifopt;
    using namespace std;
    using VectorXd = Eigen::VectorXd;
    using MatrixXd = Eigen::MatrixXd;

    class myVariables : public VariableSet{

    public:
        // Every variable set has a name, here "var_set1". this allows the constraints
        // and costs to define values and Jacobians specifically w.r.t this variable set.
        myVariables() : myVariables("var_set2"){};
        myVariables(const std::string& name): VariableSet(5*N,name), x_(5,N){
            int c = 0;
            //initial values from where NLP starts iterating
            for (int i = 0; i < N; i++){
                // x stands for state
                x_(0,i) = (double)(i/(N-1.0))*d;
                c += 1;
                x_(1,i) = (double)(i/(N-1.0))*pi;
                c += 1;// x_dot
                x_(2,i) = 0.0;
                c += 1;
                x_(3,i) = 0.0;
                c += 1;//control force: u
                x_(4,i) = 0.0;
                c += 1;
            }
            cout<<x_<<endl; 
            // cout<<"number of variables"<<c<<endl;
        }
        //here is where you start transforming Eigen::Vector into whatever
        //internal representation of your vectors you have(here two doubles, but
        // can also be complex classes such as splines, etc..
        void SetVariables(const VectorXd& x)override{
            // for (int i = 0, j = 0; i < 5*N; i+=5, j+=1){
            //     x_.col(j) = x.segment(i,i+4);
            // }
            VectorXd x1 = x;
            // cout<<"size of vector"<<x.size();
            Eigen::Map<MatrixXd> x_(x1.data(),5,N);
            
        };
        // Here is the reverse transformation from the internal representation to
        // to the Eigen::Vector
        VectorXd GetValues() const override
        {   
            Eigen::Map<const VectorXd> v1(x_.data(), x_.size());
            // cout<<x_<<endl<<"Vector";
            // cout<<v1<<endl;
            return v1;            
        };
        
        // Each variable has an upper and lower bound set here
        VecBound GetBounds() const override
        {
            VecBound bounds(5*N);
            for (int i = 0; i < N;i++){
                bounds.at(i)= Bounds(-10.0,10.0);  
                bounds.at(i+1)= Bounds(-pi,pi);  
                bounds.at(i+2)= Bounds(-20.0, 20.0);
                bounds.at(i+3)= Bounds(-20.0, 20.0);
                bounds.at(i+4)= Bounds(-200.0, 200.0);
            }
            // cout<<"size of vector"<<bounds.size()<<endl; 
            return bounds;
        }

    private:
        MatrixXd x_;

    };

    class myConstraint : public ConstraintSet{
    public:

        VectorXd dynamic(const MatrixXd x) const
        {
            MatrixXd x_dot(4,N);//easier to compute in Matrix form
            for (int i = 0; i<N; i += 1)
            {
                x_dot(0,i) = x(2,i);
                x_dot(1,i) = x(3,i);
                x_dot(2,i) = (l*m2*sin(x(1,i))*pow(x(3,i),2) + x(4,i) + m2*g*cos(x(1,i))*sin(x(1,i)))/(m1 + m2*pow(sin(x(1,i)),2));
                x_dot(3,i) = l*m2*cos(x(1,i))*sin(x(3,i))*pow(x(4,i),2) + x(4,i)*cos(x(1,i)) +(m1+m2)*g*sin(x(1,i));
                x_dot(3,i) = -x_dot(3,i)/(l*m1 + l*m2*pow(sin(x(1,i)),2));
            }    
            Eigen::Map<VectorXd> xd(x_dot.data(),x_dot.size());
            // cout<<"size of vector"<<xd.size()<<endl;
            return xd;
        }

        myConstraint() : myConstraint("constraint2") {}
        // This constraint set just contains 1 constraint, however generally
        // each set can contain multiple related constraints.
        myConstraint(const std::string& name) : ConstraintSet(2, name) {}

        // The constraint value minus the constant value "1", moved to bounds.
        VectorXd GetValues() const override
        {
            // int number = ((4*2)+(2*N));//total number of constraints
            VectorXd c(8);
            VectorXd dump(4*N);
            // VectorXd dumpv(4*(N-1));
            VectorXd x1 = GetVariables()->GetComponent("var_set2")->GetValues();
            Eigen::Map<MatrixXd> x(x1.data(),5,N); //easier to compute
            // constraint for dynamics
            dump = dynamic(x);
            for (int i = 0; i < N-1; i++){
                // cout<<"Not Worked \n";
                dump(i) = (x(0,i+1) - x(0,i)) - (h*(dump(i+4) - dump(i)/2.0));
                dump(1+i) = (x(1,i+1) - x(1,i)) - (h*(dump(i+4+1) - dump(1+i)/2.0));
                dump(2+i) = (x(2,i+1) - x(2,i)) - (h*(dump(i+4+1) - dump(2+i)/2.0));
                dump(3+i) = (x(3,i+1) - x(3,i)) - (h*(dump(i+4+1) - dump(3+i)/2.0));
                // cout<<"Worked \n";
            }
            // constraint for initial and final state
            c(0) = x(0,1);
            c(1) = x(1,1);
            c(2) = x(2,1);
            c(3) = x(3,1);
            c(4) = x(0,N-1) - d;
            c(5) = x(1,N-1) - pi;
            c(6) = x(2,N-1);
            c(7) = x(3,N-1);
            
            // Eigen::Map<VectorXd> dumpv(dump.data(), dump.size());

            VectorXd final(dump.size()+c.size());
            final << dump , c;
            // cout<<"final size "<<final.size()<<" c + dump "<<c.size()+dump.size()<<endl;
            return final;
        };
        // The only constraint in this set is an equality constraint to 1.
        // Constant values should always be put into GetBounds(), not GetValues().
        // For inequality constraints (<,>), use Bounds(x, inf) or Bounds(-inf, x).   
        VecBound GetBounds() const override         
        {
            VecBound b(8+ 4*N);
            for (int i=0;i<(8 + 4*N);i++){
                b.at(i) = Bounds(0.0, 0.0); 
                // cout<<"for loop count"<<GetRows()<<endl;
            }
            return b;
        }
        // This function provides the first derivative of the constraints.
        // In case this is too difficult to write, you can also tell the solvers to
        // approximate the derivatives by finite differences and not overwrite this
        // function, e.g. 
        // use_jacobian_approximation_ = true
        void FillJacobianBlock (std::string var_set, Jacobian& jac_block) const override
        {
        // must fill only that submatrix of the overall Jacobian that relates
        // to this constraint and "var_set1". even if more constraints or variables
        // classes are added, this submatrix will always start at row 0 and column 0,
        // thereby being independent from the overall problem.
        // if (var_set == "var_set2") {
        //     VectorXd x = GetVariables()->GetComponent("var_set2")->GetValues();

        //     jac_block.coeffRef(0, 0) = 2.0*x(0); // derivative of first constraint w.r.t x0
        //     jac_block.coeffRef(0, 1) = 1.0;      // derivative of first constraint w.r.t x1
        // }
        }
    };

    class myCost: public CostTerm {
    public:
        myCost() : myCost("cost_term2"){}
        myCost(const std::string& name) : CostTerm(name) {}
        double GetCost() const override
        {
            VectorXd x = GetVariables()->GetComponent("var_set2")->GetValues();
            
            double result=0;
            for(int i=0 ; i < N; i+=1)
            {
                result = result + (h/2)*(pow(x(i+4),2)+pow(x(i+4+5),2));
            } 
            cout<<"result : "<<result<<endl;
            return result;
        };
                    
        void FillJacobianBlock (std::string var_set, Jacobian& jac) const override
        {
            // if (var_set == "var_set1") {
            // VectorXd x = GetVariables()->GetComponent("var_set1")->GetValues();

            // jac.coeffRef(0, 0) = 0.0;             // derivative of cost w.r.t x0
            // jac.coeffRef(0, 1) = -2.0*(x(1)-2.0); // derivative of cost w.r.t x1
            // }
        }
    };

}

#endif
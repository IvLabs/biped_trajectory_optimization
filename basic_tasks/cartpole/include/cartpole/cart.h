#ifndef CART_H
#define CART_H
#define N 50
#define d 1
#define pi 3.14159265
#define l 0.5
#define m1 10
#define m2 0.3
#define g 9.81
#define T 2

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

using VectorXd = Eigen::VectorXd;
using MatrixXd = Eigen::MatrixXd;

class myVariables : public VariableSet{
public:
    // Every variable set has a name, here "var_set1". this allows the constraints
    // and costs to define values and Jacobians specifically w.r.t this variable set.
    myVariables() : myVariables("var_set2"){};
    myVariables(const std::string& name): VariableSet(5*N,name), x_(5,N){
        //initial values from where NLP starts iterating
        for (int i = 0; i < N; i++){
            // x stands for state
            x_(0,i) = (double)(i/(N-1))*d;
            x_(1,i) = (double)(i/(N-1))*pi;
            // x_dot
            x_(2,i) = 0;
            x_(3,i) = 0;
            //control force: u
            x_(4,i) = 0;
        } 
    }
    //here is where you start transforming Eigen::Vector into whatever
    //internal representation of your vectors you have(here two doubles, but
    // can also be complex classes such as splines, etc..
    void SetVariables(const VectorXd& x) override{
        for (int i = 0;i < 5*N;i+=5){
            for (int j = i; j < 5;j++){ 
                x_(j,i) = x(j);
            }
        }
    };

    VectorXd GetValues() const override
    {    
        return x_;            
    };
    
    // Each variable has an upper and lower bound set here
    VecBound GetBounds() const override
    {
        VecBound bounds(GetRows());
        
        
        bounds.at(0)= Bounds(-10.0,10.0);  
        bounds.at(1)= Bounds(-2*pi,2*pi);  
        bounds.at(2)= NoBound;
        bounds.at(3)= NoBound;
        bounds.at(4)= Bounds(-10, 10);
            
        return bounds;
    }


private:
    MatrixXd x_;
};

class myConstraint : public ConstraintSet{
public:

    MatrixXd dynamic(MatrixXd x)
    {
        MatrixXd x_dot(4,N);
        for (int i=0; i<N; i++)
        {
            
            x_dot(0,i) = x(2,i);
            x_dot(1,i) = x(3,i);
            x_dot(2,i) = (l*m2*sin(x(1,i))*pow(x(3,i),2) + x(4,i) + m2*g*cos(x(1,i))*sin(x(1,i)))/(m1 + m2*pow(sin(x(1,i)),2));
            x_dot(3,i) = l*m2*cos(x(1,i))*sin(x(3,i))*pow(x(4,i),2) + x(4,i)*cos(x(1,i)) +(m1+m2)*g*sin(x(1,i));
            x_dot(3,i) = -x_dot(3,i)/(l*m1 + l*m2*pow(sin(x(1,i)),2));
        }    
        return x_dot;
    }

    myConstraint() : myConstraint("constraint2") {}
    // This constraint set just contains 1 constraint, however generally
    // each set can contain multiple related constraints.
    myConstraint(const std::string& name) : ConstraintSet(2, name) {}

    // The constraint value minus the constant value "1", moved to bounds.
    MatrixXd GetValues()
    {
        int number = ((4*2)+(2*N));//total number of constraints
        MatrixXd c(4,2);
        MatrixXd dummy(4,N);
        MatrixXd x = GetVariables()->GetComponent("var_set2")->GetValues();
        // constraint for dynamics
        dummy = dynamic(x);
        // constraint for initial and final state
        c(1,1) = x(1,1);
        c(2,1) = x(2,1);
        c(3,1) = x(3,1);
        c(4,1) = x(4,1);
        c(1,2) = x(1,N-1) - 5;
        c(2,2) = x(2,N-1) - pi;
        c(3,2) = x(3,N-1);
        c(4,2) = x(4,N-1);
        
        return c;
    }

    VecBound GetBounds() const override
    {
        VecBound b(GetRows());
        b.at(0) = Bounds(1.0, 1.0);
        return b;
    }
    // This function provides the first derivative of the constraints.
    // In case this is too difficult to write, you can also tell the solvers to
    // approximate the derivatives by finite differences and not overwrite this
    // function, e.g. 
    // use_jacobian_approximation_ = true
    // void FillJacobianBlock (std::string var_set, Jacobian& jac_block) const override
    // {
    // // must fill only that submatrix of the overall Jacobian that relates
    // // to this constraint and "var_set1". even if more constraints or variables
    // // classes are added, this submatrix will always start at row 0 and column 0,
    // // thereby being independent from the overall problem.
    // if (var_set == "var_set2") {
    //     VectorXd x = GetVariables()->GetComponent("var_set2")->GetValues();

    //     jac_block.coeffRef(0, 0) = 2.0*x(0); // derivative of first constraint w.r.t x0
    //     jac_block.coeffRef(0, 1) = 1.0;      // derivative of first constraint w.r.t x1
    // }
        
};

class myCost: public CostTerm {
public:
    myCost() : myCost("cost_term2"){}
    myCost(const std::string& name) : CostTerm(name) {}
    double GetCost() const override
    {
        MatrixXd x = GetVariables()->GetComponent("var_set2")->GetValues();
        double h = T/(N-1);
        double result=0;
        for(int i=0 ; i <N-1 ; i++)
        {
            result = result + (h/2)*(pow(x(4,i),2)+pow(x(4,i+1),2));
        } 
        return result;
    };
                
    
};



#endif
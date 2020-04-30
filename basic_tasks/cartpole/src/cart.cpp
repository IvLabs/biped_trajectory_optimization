#include <iostream>
#include <ifopt/ipopt_solver.h>
#include "cartpole/cart.h"

using namespace ifopt;
using namespace cartpole;

int main()
{
    // 1. define the problem
    Problem nlp;
    nlp.AddVariableSet  (std::make_shared<myVariables>());
    nlp.AddConstraintSet(std::make_shared<myConstraint>());
    nlp.AddCostSet      (std::make_shared<myCost>());
    nlp.PrintCurrent();

    // 2. choose solver and options
    IpoptSolver ipopt;
    ipopt.SetOption("linear_solver", "mumps");
    ipopt.SetOption("jacobian_approximation", "finite-difference-values");
    ipopt.SetOption("fixed_variable_treatment","relax_bounds");
    
    // 3 . solve
    ipopt.Solve(nlp);
    Eigen::VectorXd x = nlp.GetOptVariables()->GetValues();
    std::cout << x.transpose() << std::endl;
}
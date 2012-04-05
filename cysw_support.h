#include "stdint.h"
#include <string>
#include <vector>
#include "cassowary/SimplexSolver.h"
#include "cassowary/LinearExpression.h"
#include "cassowary/Constraint.h"


std::vector<size_t> get_cpp_exception_constraint_pointers();
std::string get_cpp_exception_message();
std::string solver_str(SimplexSolver *solver);
P_Constraint *newLinearEquation(const P_LinearExpression &lhs, const P_LinearExpression &rhs, const Strength &strength, double weight);
P_Constraint *newLinearInequality(const P_LinearExpression &lhs, CnRelation op, const P_LinearExpression &rhs, const Strength &strength, double weight);
P_LinearExpression newLinearExpression(double constant);
void delete_P_Constraint(P_Constraint *pcn);
size_t get_P_Constraint_addr(P_Constraint *pcn);

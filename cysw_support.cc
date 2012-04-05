// Support functions for cysw.

#include <sstream>

#include "cysw_support.h"

#include "cassowary/LinearEquation.h"
#include "cassowary/LinearInequality.h"
#include "cassowary/Errors.h"


size_t get_P_Constraint_addr(P_Constraint *pcn) {
    return reinterpret_cast<size_t>(pcn->ptr());
}

std::vector<size_t> get_cpp_exception_constraint_pointers() {
    std::vector<size_t> constraint_pointers;
    try {
        // Re-throw the exception.
        throw;
    }
    catch (const ExCLError &exn) {
        const ConstraintSet* cset = exn.explanation();
        if (cset != NULL) {
            ConstraintSet::const_iterator it = cset->begin();
            for (; it != cset->end(); ++it) {
                constraint_pointers.push_back(reinterpret_cast<size_t>((*it).ptr()));
            }
        }
    }
    return constraint_pointers;
}

std::string get_cpp_exception_message() {
    std::string message;
    try {
        // Re-throw the exception.
        throw;
    }
    catch (const ExCLError &exn) {
        message = exn.description();
    }
    return message;
}

std::string solver_str(SimplexSolver *solver) {
    std::ostringstream ss;
    solver->PrintOn(ss);
    return ss.str();
}


P_Constraint *newLinearEquation(const P_LinearExpression &lhs, const P_LinearExpression &rhs, const Strength &strength, double weight) {
    P_Constraint result(new LinearEquation(lhs, rhs, strength, weight));
    return new P_Constraint(static_cast< const P_Constraint& >(result));
}
P_Constraint *newLinearInequality(const P_LinearExpression &lhs, CnRelation op, const P_LinearExpression &rhs, const Strength &strength, double weight) {
    P_Constraint result(new LinearInequality(lhs, op, rhs, strength, weight));
    return new P_Constraint(static_cast< const P_Constraint& >(result));
}
P_LinearExpression newLinearExpression(double constant) {
    return P_LinearExpression(new LinearExpression(constant));
}

void delete_P_Constraint(P_Constraint *pcn) {
    delete pcn;
}

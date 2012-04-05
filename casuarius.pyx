""" Cython wrapper for Cassowary.
"""

from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc
from libc.math cimport fabs
from libcpp.string cimport string
from libcpp.vector cimport vector

import operator
from collections import defaultdict


cdef extern from "cysw_support.h":
    string get_cpp_exception_message()
    vector[size_t] get_cpp_exception_constraint_pointers()


class CassowaryError(Exception):
    """ An error raised from Cassowary.
    """

class ExplainedCassowaryError(CassowaryError):
    """ An error raised from Cassowary with a list of constraints that
    are involved in the error.

    Use the .explains() method with a list of LinearConstraints to get
    a filtered list of LinearConstraints that explain this error.
    """
    # SORRY: Yes, it's an awkward API. Best I can do right now.

    def explains(self, linear_constraints):
        """ Filter a list of LinearConstraints to return only those
        LinearConstraints that caused the error.
        """
        filtered = []
        for cn in linear_constraints:
            if cn.explains_exception(self):
                filtered.append(cn)
        return filtered


cdef int raise_cassowary_error() except *:
    cdef string message
    cdef vector[size_t] constraint_pointers

    message = get_cpp_exception_message()
    constraint_pointers = get_cpp_exception_constraint_pointers()
    if constraint_pointers.size() > 0:
        e = ExplainedCassowaryError(message.c_str())
        e.explanation_constraint_pointers = []
    else:
        e = CassowaryError(message.c_str())

    cdef vector[size_t].iterator it = constraint_pointers.begin()
    while it != constraint_pointers.end():
        e.explanation_constraint_pointers.append(deref(it))
        inc(it)

    raise e


cdef extern from "cassowary/SymbolicWeight.h":
    cdef cppclass ClSymbolicWeight "SymbolicWeight":
        ClSymbolicWeight (vector[double] weights)
        ClSymbolicWeight operator*(double n)
        ClSymbolicWeight operator/(double n)
        ClSymbolicWeight operator+(ClSymbolicWeight cl)
        double AsDouble()

cdef extern from "cassowary/Strength.h":
    cdef cppclass ClStrength "Strength":
        ClStrength(string Name, ClSymbolicWeight symbolicWeight, bint isRequired)
        ClSymbolicWeight symbolicWeight()

cdef extern from "cassowary/Variable.h":
    cdef cppclass ClVariable "Variable":
          #ClVariable(ClVariable clv_)
          ClVariable(string name, double Value) 
          #ClVariable(double Value=0.0) 
          double Value()
          void SetValue(double value)

cdef extern from "cassowary/LinearExpression.h":
    cdef cppclass ClLinearExpression "LinearExpression":
        ClLinearExpression(double num)
        #ClLinearExpression(ClVariable clv, double value=0.0, double constant=0.0)
        #ClLinearExpression(ClLinearExpression expr)
        #ClLinearExpression Times(double x)
        #ClLinearExpression Divide(double x)
        #ClLinearExpression Plus(ClLinearExpression expr)
        #ClLinearExpression Minus(ClLinearExpression expr)
        ClLinearExpression AddVariable(ClVariable v, double c)

    ctypedef ClLinearExpression* P_LinearExpression

cdef extern from "cassowary/Constraint.h":
    cdef cppclass ClConstraint "Constraint":
        void ChangeStrength(ClStrength strength) except +raise_cassowary_error
        void ChangeWeight(double weight) except +raise_cassowary_error

    ctypedef ClConstraint* P_Constraint

    cdef enum ClCnRelation "CnRelation":
        cnEQ = 0
        cnNEQ = 100
        cnLEQ = 2
        cnGEQ = -2
        cnLT = 3
        cnGT = -3

cdef extern from "cysw_support.h":
    size_t get_P_Constraint_addr(P_Constraint *pcn)

cdef extern from "cassowary/LinearEquation.h":
    cdef cppclass ClLinearEquation "LinearEquation":
        pass
        #ClLinearEquation(ClLinearExpression cle, ClStrength strength, double weight)
        #ClLinearEquation(ClLinearExpression cle1, ClLinearExpression cle2, ClStrength strength, double weight)

    ctypedef ClLinearEquation* P_LinearEquation

cdef extern from "cassowary/LinearInequality.h":
    cdef cppclass ClLinearInequality "LinearInequality":
        pass
        #ClLinearInequality(ClLinearExpression cle, ClStrength strength, double weight=1.0)
        #ClLinearInequality(ClLinearExpression cle1, ClCnRelation op, ClLinearExpression cle2, ClStrength strength, double weight)

    ctypedef ClLinearInequality* P_LinearInequality

cdef extern from "cassowary/SimplexSolver.h":
    cdef cppclass ClSimplexSolver "SimplexSolver":
        ClSimplexSolver()
        # Note: most of these void return types actually should be
        # ClSimplexSolver&, but we don't use the return values, and it causes
        # problems in the generated C++.
        void AddConstraint(P_Constraint pcn) except +raise_cassowary_error
        void RemoveConstraint(P_Constraint pcn) except +raise_cassowary_error
        void AddEditVar(ClVariable v, ClStrength strength, double weight) except +raise_cassowary_error
        void RemoveEditVar(ClVariable v) except +raise_cassowary_error
        void BeginEdit() except +raise_cassowary_error
        void EndEdit() except +raise_cassowary_error
        void RemoveAllEditVars()
        void Solve() except +raise_cassowary_error
        bint SolveNoException()
        void Resolve() except +raise_cassowary_error
        void SetAutosolve(bint f) except +raise_cassowary_error
        void SuggestValue(ClVariable v, double x) except +raise_cassowary_error
        void Reset()
        void SetExplaining(bint f)
        bint FIsExplaining()
        void Solve()
        bint FContainsVariable(ClVariable v) except +raise_cassowary_error

cdef extern from "cysw_support.h":
    string solver_str(ClSimplexSolver *solver)
    P_Constraint *newLinearEquation(P_LinearExpression lhs, P_LinearExpression rhs, ClStrength strength, double weight)
    P_Constraint *newLinearInequality(P_LinearExpression lhs, ClCnRelation op, P_LinearExpression rhs, ClStrength strength, double weight)
    P_LinearExpression newLinearExpression(double constant)
    void delete_P_Constraint(P_Constraint *pcn)

cdef class SymbolicWeight:
    cdef ClSymbolicWeight *symbolic_weight
    cdef readonly tuple weights

    def __cinit__(self, tuple weights):
        cdef vector[double] c_weights
        cdef double c_weight
        cdef int i

        for i from 0 <= i < len(weights):
            c_weights.push_back(<double>weights[i])
        self.symbolic_weight = new ClSymbolicWeight(c_weights)
        self.weights = weights
    
    def __dealloc__(self):
        del self.symbolic_weight

    def __float__(self):
        return self.symbolic_weight.AsDouble()

    def __mul__(self, double other):
        new_weights = tuple(x*other for x in self.weights)
        return SymbolicWeight(new_weights)

    def __rmul__(self, double other):
        return self.__mul__(other)

    def __repr__(self):
        return 'SymbolicWeight({0!r})'.format(self.weights)

    def __hash__(self):
        return id(self)


cdef class Strength:
    cdef ClStrength *strength
    cdef readonly bytes name
    cdef readonly SymbolicWeight symbolic_weight
    cdef readonly bint is_required

    def __cinit__(self, bytes name, SymbolicWeight symbolic_weight, bint is_required=False):
        self.strength = new ClStrength(string(<char*>name), deref(symbolic_weight.symbolic_weight), is_required)
        self.name = name
        self.symbolic_weight = symbolic_weight
        self.is_required = is_required

    def __dealloc__(self):
        del self.strength

    def __repr__(self):
        if self.is_required:
            template = 'Strength({name!r}, {symbolic_weight!r}, is_required={is_required})'
        else:
            template = 'Strength({name!r}, {symbolic_weight!r})'
        return template.format(name=self.name, symbolic_weight=self.symbolic_weight, is_required=self.is_required)

    def __str__(self):
        return repr(self.name)

    def __hash__(self):
        return id(self)


required = Strength('required', SymbolicWeight((1000.0, 1000.0, 1000.0)), is_required=True)
strong = Strength('strong', SymbolicWeight((1.0, 0.0, 0.0)))
medium = Strength('medium', SymbolicWeight((0.0, 1.0, 0.0)))
weak = Strength('weak', SymbolicWeight((0.0, 0.0, 1.0)))

STRENGTH_MAP = dict(
    required=required,
    strong=strong,
    medium=medium,
    weak=weak,
)


cdef bint almost_equal(double a, double b, double eps=1e-8):
    return fabs(a - b) < eps


cdef class LinearSymbolic:

    def nonlinear(self, msg):
        raise TypeError('Non-linear expression: %s' % msg)
    
    def nonlinear_op(self, op):
        raise TypeError('Non-linear operator: `%s`' % op)

    def convert_to_csw(self):
        raise NotImplementedError

    def __repr__(self):
        raise NotImplementedError

    def __add__(self, other):
        raise NotImplementedError
    
    def __mul__(self, other):
        raise NotImplementedError
    
    def __div__(self, other):
        raise NotImplementedError

    def __str__(self):
        return self.__repr__()

    def __sub__(self, other):
        return self + (-1.0 * other)
    
    def __floordiv__(self, other):
        if not isinstance(self, LinearSymbolic):
            self = other
        self.nonlinear_op('//')

    
    def __mod__(self, other):
        if not isinstance(self, LinearSymbolic):
            self = other
        self.nonlinear_op('%')
    
    def __divmod__(self, other):
        if not isinstance(self, LinearSymbolic):
            self = other
        self.nonlinear_op('divmod')
    
    def __pow__(self, other, mod):
        if not isinstance(self, LinearSymbolic):
            self = other
        self.nonlinear_op('**')
    
    def __lshift__(self, other):
        if not isinstance(self, LinearSymbolic):
            self = other
        self.nonlinear_op('<<')

    def __rshift__(self, other):
        if not isinstance(self, LinearSymbolic):
            self = other
        self.nonlinear_op('>>')
    
    def __and__(self, other):
        if not isinstance(self, LinearSymbolic):
            self = other
        self.nonlinear_op('&')
    
    def __or__(self, other):
        if not isinstance(self, LinearSymbolic):
            self = other
        self.nonlinear_op('|')

    def __xor__(self, other):
        if not isinstance(self, LinearSymbolic):
            self = other
        self.nonlinear_op('^')
    
    def __richcmp__(self, other, int op):
        if op == 0 or op == 4 or op == 3:
            raise ValueError('Invalid constraint operation.')
        elif op == 2:
            return self._richcmp_eq(other)
        elif op == 1:
            return self._richcmp_le(other)
        elif op == 5:
            return self._richcmp_ge(other)

    cpdef LinearConstraint _richcmp_le(self, other):
        if isinstance(other, (float, int, long)):
            rhs = LinearExpression([], float(other))
        elif isinstance(other, LinearSymbolic):
            rhs = other
        else:
            msg = 'Invalid type for constraint operation %s' % type(other)
            raise TypeError(msg)
        return LEConstraint(self, rhs)
    
    cpdef LinearConstraint _richcmp_ge(self, other):
        if isinstance(other, (float, int, long)):
            rhs = LinearExpression([], float(other))
        elif isinstance(other, LinearSymbolic):
            rhs = other
        else:
            msg = 'Invalid type for constraint operation %s' % type(other)
            raise TypeError(msg)
        return GEConstraint(self, rhs)

    cpdef LinearConstraint _richcmp_eq(self, other):
        if isinstance(other, (float, int, long)):
            rhs = LinearExpression([], float(other))
        elif isinstance(other, LinearSymbolic):
            rhs = other
        else:
            msg = 'Invalid type for constraint operation %s' % type(other)
            raise TypeError(msg)
        return EQConstraint(self, rhs)

    def __hash__(self):
        return id(self)


cdef class ConstraintVariable(LinearSymbolic):
    cdef ClVariable *variable
    cdef readonly bytes name

    def __cinit__(self, bytes name, double value=0.0):
        self.variable = new ClVariable(string(<char*>name), value)
        self.name = name

    property value:
        def __get__(self):
            return self.variable.Value()

    def __dealloc__(self):
        del self.variable

    def __repr__(self):
        return 'ConstraintVariable({0!r}, {1!r})'.format(self.name, self.variable.Value())

    def __str__(self):
        return '{0}:{1}'.format(self.name, self.value)

    def __add__(self, other):
        if not isinstance(self, LinearSymbolic):
            self, other = other, self
        if isinstance(other, (float, int, long)):
            terms = [Term(self)]
            const = float(other)
            expr = LinearExpression(terms, const)
        elif isinstance(other, Term):
            terms = [Term(self), other]
            expr = LinearExpression(terms)
        elif isinstance(other, ConstraintVariable):
            terms = [Term(self), Term(other)]
            expr = LinearExpression(terms)
        elif isinstance(other, LinearExpression):
            expr = other + self
        else:
            return NotImplemented
        return expr
    
    def __mul__(self, other):
        if not isinstance(self, LinearSymbolic):
            self, other = other, self
        if isinstance(other, (float, int, long)):
            res = Term(self, float(other))
        elif isinstance(other, (Term, ConstraintVariable, LinearExpression)):
            self.nonlinear('[ %s ] * [ %s ]' % (self, other))
        else:
            return NotImplemented
        return res
    
    def __div__(self, other):
        if not isinstance(self, LinearSymbolic):
            other.nonlinear('[ %s ] / [ %s ]' % (self, other))
        if isinstance(other, (float, int)):
            res = (1.0 / float(other)) * self
        elif isinstance(other, (Term, ConstraintVariable, LinearExpression)):
            self.nonlinear('[ %s ] / [ %s ]' % (self, other))
        else:
            return NotImplemented
        return res


cdef class Term(LinearSymbolic):
    cdef public ConstraintVariable var
    cdef public double coeff

    def __cinit__(self, ConstraintVariable var, double coeff=1.0):
        self.var = var
        self.coeff = coeff

    def __repr__(self):
        return 'Term({0!r}, {1!r})'.format(self.var, self.coeff)

    def __str__(self):
        if self.coeff == 1.0:
            template = '{name}:{value}'
        elif self.coeff == -1.0:
            template = '-{name}:{value}'
        else:
            template = '{coeff} * {name}:{value}'
        return template.format(coeff=self.coeff, name=self.var.name, value=self.var.value)

    def __add__(self, other):
        if not isinstance(self, LinearSymbolic):
            self, other = other, self
        if isinstance(other, (float, int, long)):
            terms = [self]
            const = float(other)
            expr = LinearExpression(terms, const)
        elif isinstance(other, Term):
            terms = [self, other]
            expr = LinearExpression(terms)
        elif isinstance(other, ConstraintVariable):
            terms = [self, Term(other)]
            expr = LinearExpression(terms)
        elif isinstance(other, LinearExpression):
            expr = other + self
        else:
            return NotImplemented
        return expr
    
    def __mul__(self, other):
        if not isinstance(self, LinearSymbolic):
            self, other = other, self
        if isinstance(other, (float, int, long)):
            res = Term(self.var, float(other) * self.coeff)
        elif isinstance(other, (Term, ConstraintVariable, LinearExpression)):
            self.nonlinear('[ %s ] * [ %s ]' % (self, other))
        else:
            return NotImplemented
        return res
    
    def __div__(self, other):
        if not isinstance(self, LinearSymbolic):
            other.nonlinear('[ %s ] / [ %s ]' % (self, other))
        if isinstance(other, (float, int, long)):
            res = (1.0 / float(other)) * self
        elif isinstance(other, (Term, ConstraintVariable, LinearExpression)):
            self.nonlinear('[ %s ] / [ %s ]' % (self, other))
        else:
            return NotImplemented
        return res


cdef class LinearExpression(LinearSymbolic):
    cdef public tuple terms
    cdef public double constant

    def reduce_terms(self, terms):
        mapping = defaultdict(float)
        for term in terms:
            mapping[term.var] += term.coeff
        terms = tuple(Term(var, coeff) for (var, coeff) in mapping.iteritems() 
                      if not almost_equal(coeff, 0.0))
        return terms

    def __cinit__(self, terms, constant=0.0):
        self.terms = self.reduce_terms(terms)
        self.constant = constant

    property value:
        def __get__(self):
            cdef double value=self.constant
            for term in self.terms:
                value += term.coeff * term.var.value
            return value

    def __repr__(self):
        if len(self.terms) > 0:
            s = sorted(self.terms, key=operator.attrgetter('var.name'))
            terms = ' + '.join(str(term) for term in s)
            if self.constant > 0.0:
                terms += ' + %s' % self.constant
            elif self.constant < 0.0:
                terms += ' - %s' % -self.constant
        else:
            terms = str(self.constant)
        return terms

    def __add__(self, other):
        if not isinstance(self, LinearSymbolic):
            self, other = other, self
        if isinstance(other, (float, int, long)):
            expr = LinearExpression(self.terms, self.constant + float(other))
        elif isinstance(other, Term):
            terms = list(self.terms) + [other]
            expr = LinearExpression(terms, self.constant)
        elif isinstance(other, ConstraintVariable):
            terms = list(self.terms) + [Term(other)]
            expr = LinearExpression(terms, self.constant)
        elif isinstance(other, LinearExpression):
            terms = list(self.terms) + list(other.terms)
            const = self.constant + other.constant
            expr = LinearExpression(terms, const)
        else:
            return NotImplemented
        return expr
    
    def __mul__(self, other):
        if not isinstance(self, LinearSymbolic):
            self, other = other, self
        if isinstance(other, (float, int, long)):
            terms = [other * term for term in self.terms]
            const = self.constant * other
            res = LinearExpression(terms, const)
        elif isinstance(other, (Term, ConstraintVariable, LinearExpression)):
            self.nonlinear('[ %s ] * [ %s ]' % (self, other))
        else:
            return NotImplemented
        return res
    
    def __div__(self, other):
        if not isinstance(self, LinearSymbolic):
            self, other = other, self
        if isinstance(other, (float, int, long)):
            res = (1.0 / float(other)) * self
        elif isinstance(other, (Term, ConstraintVariable, LinearExpression)):
            self.nonlinear('[ %s ] / [ %s ]' % (self, other))
        else:
            return NotImplemented
        return res

    cdef P_LinearExpression as_cl_linear_expression(self):
        """ Convert to a ClLinearExpression.
        """
        cdef Term term
        cdef P_LinearExpression expr

        expr = newLinearExpression(self.constant)
        for term in self.terms:
            expr.AddVariable(deref(term.var.variable), term.coeff)
        return expr

cdef LinearExpression as_linear_expression(obj):
    """ Cast any appropriate object to a LinearExpression.
    """
    if isinstance(obj, LinearExpression):
        return obj
    elif isinstance(obj, (float, int, long)):
        obj = float(obj)
        return LinearExpression((), obj)
    elif isinstance(obj, ConstraintVariable):
        return LinearExpression((Term(obj),))
    elif isinstance(obj, Term):
        return LinearExpression((obj,))
    else:
        raise TypeError("Cannot cast {0!r}, a {1!r}, to LinearExpression.".format(obj, type(obj)))


cdef class LinearConstraint:
    cdef readonly LinearExpression lhs
    cdef readonly LinearExpression rhs
    cdef Strength _strength
    cdef double _weight
    cdef readonly bytes op
    cdef P_Constraint *cl_linear_constraint

    def __cinit__(self, lhs, rhs, *args, **kwds):
        self.lhs = as_linear_expression(lhs)
        self.rhs = as_linear_expression(rhs)

    def __init__(self, lhs, rhs, Strength strength=required, double weight=1.0):
        self._strength = strength
        self._weight = weight
        self.cl_linear_constraint = self.as_cl_linear_constraint()

    def __dealloc__(self):
        delete_P_Constraint(self.cl_linear_constraint)

    def __repr__(self):
        return '<%s: (%s %s %s) | %s | %s>' % (type(self).__name__, self.lhs,
            self.op, self.rhs, self.strength, self.weight)

    def __str__(self):
        return '%s %s %s' % (self.lhs, self.op, self.rhs)


    property strength:
        def __get__(self):
            return self._strength
        def __set__(self, Strength strength):
            self._strength = strength
            deref(self.cl_linear_constraint).ChangeStrength(deref(strength.strength))

    property weight:
        def __get__(self):
            return self._weight
        def __set__(self, double weight):
            self._weight = weight
            deref(self.cl_linear_constraint).ChangeWeight(weight)

    cpdef bint explains_exception(self, object e):
        """ Returns True if this constraint is in the explanation list
        of the exception.
        """
        if self.cl_linear_constraint == NULL or not isinstance(e, CassowaryError):
            return False
        cdef size_t addr = get_P_Constraint_addr(self.cl_linear_constraint)
        cdef size_t ptr
        for ptr in e.explanation_constraint_pointers:
            if ptr == addr:
                return True
        return False

    def __or__(self, other):
        cdef Strength strength = self.strength
        cdef double weight = self.weight

        if isinstance(other, (float, int, long)):
            weight = other
        elif isinstance(other, basestring):
            if other in STRENGTH_MAP:
                strength = STRENGTH_MAP[other]
            else:
                msg = "Expected a known strength string or a Strength. Got {!r} instead.".format(other)
                raise ValueError(msg)
        elif isinstance(other, Strength):
            strength = other
        else:
            msg = "Expected a known strength string, a Strength or a float weight. Got {!r} instead.".format(other)
            raise ValueError(msg)
        return type(self)(self.lhs, self.rhs, strength=strength, weight=weight)
        
    def __ror__(self, Strength other):
        return self.__or__(other)

    def __hash__(self):
        return id(self)

    cdef P_Constraint *as_cl_linear_constraint(self):
        return NULL

cdef class LEConstraint(LinearConstraint):

    def __init__(self, lhs, rhs, Strength strength=required, double weight=1.0):
        LinearConstraint.__init__(self, lhs, rhs, strength, weight)
        self.op = b'<='

    cdef P_Constraint *as_cl_linear_constraint(self):
        cdef P_Constraint *inequality
        cdef P_LinearExpression lhs_le
        cdef P_LinearExpression rhs_le

        inequality = newLinearInequality(self.lhs.as_cl_linear_expression(), cnLEQ, self.rhs.as_cl_linear_expression(), deref(self._strength.strength), self._weight)
        return inequality

    property error:
        def __get__(self):
            return max(self.lhs.value - self.rhs.value, 0.0)

cdef class GEConstraint(LinearConstraint):

    def __init__(self, lhs, rhs, Strength strength=required, double weight=1.0):
        LinearConstraint.__init__(self, lhs, rhs, strength, weight)
        self.op = b'>='

    cdef P_Constraint *as_cl_linear_constraint(self):
        cdef P_Constraint *inequality
        cdef P_LinearExpression lhs_le
        cdef P_LinearExpression rhs_le

        inequality = newLinearInequality(self.lhs.as_cl_linear_expression(), cnGEQ, self.rhs.as_cl_linear_expression(), deref(self._strength.strength), self._weight)
        return inequality

    property error:
        def __get__(self):
            return max(self.rhs.value - self.lhs.value, 0.0)

cdef class EQConstraint(LinearConstraint):

    def __init__(self, lhs, rhs, Strength strength=required, double weight=1.0):
        LinearConstraint.__init__(self, lhs, rhs, strength, weight)
        self.op = b'=='

    cdef P_Constraint *as_cl_linear_constraint(self):
        cdef P_Constraint *equation
        cdef P_LinearExpression lhs_le
        cdef P_LinearExpression rhs_le

        equation = newLinearEquation(self.lhs.as_cl_linear_expression(), self.rhs.as_cl_linear_expression(), deref(self._strength.strength), self._weight)
        return equation

    property error:
        def __get__(self):
            return abs(self.rhs.value - self.lhs.value)

cdef class Solver:
    cdef ClSimplexSolver *solver
    cdef bint _autosolve
    cdef bint _explaining

    def __cinit__(self, bint autosolve=False, bint explaining=True):
        self.solver = new ClSimplexSolver()
        self.autosolve = autosolve
        self.explaining = explaining

    property autosolve:
        def __get__(self):
            return self._autosolve

        def __set__(self, bint autosolve):
            self._autosolve = autosolve
            self.solver.SetAutosolve(autosolve)

    property explaining:
        def __get__(self):
            return self._explaining

        def __set__(self, bint explaining):
            self._explaining = explaining
            self.solver.SetExplaining(explaining)

    def __dealloc__(self):
        del self.solver

    def __str__(self):
        return solver_str(self.solver).c_str()

    def add_constraint(self, LinearConstraint constraint):
        self.solver.AddConstraint(deref(constraint.cl_linear_constraint))

    def remove_constraint(self, LinearConstraint constraint):
        self.solver.RemoveConstraint(deref(constraint.cl_linear_constraint))

    def suggest_values(self, var_vals, Strength default_strength=strong, double default_weight=1.0):
        return SolverEditContext(self, var_vals, default_strength, default_weight)

    cdef object _begin_edit_suggest_values(self, var_vals, Strength default_strength=strong, double default_weight=1.0):
        cdef ConstraintVariable variable
        cdef double value
        cdef Strength strength
        cdef double weight

        for item in var_vals:
            strength = default_strength
            weight = default_weight
            if len(item) == 2:
                variable, value = item
            elif len(item) == 3:
                variable, value, strength = item
            elif len(item) == 4:
                variable, value, strength, weight = item
            self.solver.AddEditVar(deref(variable.variable), deref(strength.strength), weight)

        self.solver.BeginEdit()
        for item in var_vals:
            strength = default_strength
            weight = default_weight
            if len(item) == 2:
                variable, value = item
            elif len(item) == 3:
                variable, value, strength = item
            elif len(item) == 4:
                variable, value, strength, weight = item
            self.solver.SuggestValue(deref(variable.variable), value)

        self.solver.Resolve()

    cdef _end_edit(self):
        self.solver.EndEdit()


cdef class SolverEditContext:
    """ Context manager for suggesting variables in a solver.
    """
    cdef public Solver solver
    cdef public suggest_var_vals
    cdef public Strength default_strength
    cdef public double default_weight

    def __cinit__(self, solver, suggest_var_vals,
        Strength default_strength=strong, double default_weight=1.0):
        self.solver = solver
        self.suggest_var_vals = suggest_var_vals
        self.default_strength = default_strength
        self.default_weight = default_weight

    def __enter__(self):
        self.solver._begin_edit_suggest_values(self.suggest_var_vals, self.default_strength, self.default_weight)

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.solver._end_edit()

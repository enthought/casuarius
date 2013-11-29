========================================
Casuarius - Cython binding for Cassowary
========================================

The solver source code is derived from the 0.6 release of Cassowary. It has been
modified by Svilen Dobrev to remove memory leaks. His work can be found
separately as SWIG bindings:

    http://pypi.python.org/pypi/cassowarypy

The original Cassowary sources are available from its website:

    http://www.cs.washington.edu/research/constraints/cassowary/

His modifications are available under the BSD license, while Casuarius as
a whole follows the license of Cassowary itself, LGPL v2.1 or (at your option)
a later version of the LGPL. See the files LICENSE and COPYING.LGPL for details.

Please be sure to have Cython 0.15.1 or above in order to compile this extension
module. It has been tested on OS X (using llvm-gcc 4.2) and Windows (using
mingw). Other Windows C++ compilers may or may not work.

Please contact the Enthought mailing list <enthought-dev@enthought.com> for support.

    https://mail.enthought.com/mailman/listinfo/enthought-dev

Simple example::

  import casuarius

  x = casuarius.ConstraintVariable('x')
  y = casuarius.ConstraintVariable('y')
  z = casuarius.ConstraintVariable('z')
  gte = x + y >= z
  eq = z + y == 10
  solver = casuarius.Solver()
  solver.autosolve = True
  solver.add_constraint(gte)
  solver.add_constraint(eq)

  print x, y, z

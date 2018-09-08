This is a utility to compute and explore lower bounds on the minimum distances
of codes from algebraic geometric curves. It computes bounds for AG codes
coming from two-point divisors using various methods.

Reference
=========
I.Duursma, R.Kirov, and S.Park, 'Distance bounds for algebraic geometric
codes', Journal of Pure and Applied Algebra, 2010.

Installation and Usage
======================

The program consists of 3 executable Python scripts:
- build_curve_object.py
- build_data.py
- play.py

Each supports '--help' parameter.

build_curve_object.py
--------------------
This script constructs a .curve.json file containing all the information
necessary for running the minimum bound algorithms. The data in the .curve
files is encoded as JSON. Currently, the script supports the Hermitian and
Suzuki families of curves only.

build_data.py
------------
This script takes a .curve.json file and computes minimum bounds using various
methods. All computations are stored in a .full.gz file.

play.py
-------
When called in interactive mode 'python -i' (or even better 'ipython') this
script takes .full.gz data file and allows for interactive exploration of the
bounds.

The following bounds are loaded and available in the global namespace, as
function taking a divisor C. The points 'P' and 'Q' are available for
arithmetic operations.

ddk - Duursma-Kirov distance bound
ddp - Duursma-Park distance bound
db - Beelen distance bound
dabz - ABZ distance bound
dabzp - ABZP distance bound
dlm - Lundell-McCullough distance bound
dgst - Guneri-Stichtenoth-Taskin distance bound
dgop - Goppa distance bound

The following coset bounds are also available for dk, dp, b, abzp.

cp* - coset bound for C / C + P
cq* - coset bound for C / C + Q

Example
=======
The following commands will generate all bounds for two-point codes on the
Hermitian curve over F16.

./build_curve_object.py -f hermitian 4
./build_data.py hermitian16.curve.json
python -i play.py hermitian16.full.gz

>>> ddp(P + 2 * Q)
>>> cpdp(P + 2 * Q)

Coding Utilities
================
Various coding utility functions are available in codes object.
TODO: add documentation for bestcodes.py.

Requirements
======================
Python (tested on Python 2.7.2)

License
=======

Copyright (C) 2011 Iwan Duursma, Radoslav Kirov

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

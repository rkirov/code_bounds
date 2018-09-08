#!/usr/bin/python
import argparse
import gzip
import json
from curve import Curve
from coding import Coding
from bestcodes import *

def load(inputfile):
    incurve = Curve()
    incurve.load(inputfile)
    incurve.wrap_tables()
    # put curve data into global namespace
    global codes
    codes = Coding(incurve)
    globals().update(incurve._data)
    globals().update(incurve.__dict__)
    globals().update(codes.__dict__)
    global curve, c
    c = curve = incurve
    global best
    best = BestCodes(codes)
    best.load_best_codes()

if __name__ == "__main__":
    desc = '''Load the curve from a given file. Best used with python -i or
        ipython'''
    parser = argparse.ArgumentParser(prog="Interactive", description=desc)
    parser.add_argument('inputfile', type=str,
                       help='Name of the curve file to use.')
    names = vars(parser.parse_args())
    load(names['inputfile'])
    print '''
    The curve data is added to the global namespace.
    The curve is also accessible under variable name 'c' or 'curve'.
  
    Coding theory routines are in 'codes' object.
        codes.print_code(P + 2 * Q)
        codes.best_codes().get_best()
    '''
    print c
    global P, Q
    P = c.P()
    Q = c.Q()

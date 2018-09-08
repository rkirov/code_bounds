#!/usr/bin/python
from curve import Curve
import argparse
import time
import gzip
import os
import json

if __name__ == "__main__":
    starttime = time.time()
    desc = '''Build the bounds for a given curve.'''
    parser = argparse.ArgumentParser(prog="Data Builder", description=desc)
    parser.add_argument('-q', '--quiet', action='store_true', default=False,
                       help='don\'t print progress and debugging messages')
    parser.add_argument("-p", "--plain",
                  action="store_true", default=False,
                  help="Output a plain txt file (otherwise the output is a zipped file)")
    parser.add_argument('inputfile', type=str,
                       help='Name of the curve file to use.')
    names = vars(parser.parse_args())
    curve = Curve()
    curve.load_from_file(names['inputfile'])
    global verbose
    if names['quiet']:
        def verboseprint(*args):
            for arg in args:
                print arg,
            print
        verbose = verboseprint
    else:   
        verbose = lambda *a: None
    curve.build_all()
    curve.save(names['inputfile'].split('.')[0] + '.full', names['plain'])
    print "Done in %f sec" % (time.time() - starttime)

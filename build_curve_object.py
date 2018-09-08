#!/usr/bin/python
from curve_families import * 
import argparse 
import json
    
if __name__ == "__main__": 
    desc = '''Build a .curve.json file containing all data for a curve.'''
    parser = argparse.ArgumentParser(prog="Curve Builder", description=desc)
    parser.add_argument('-f', '--family', required=True, help='family of the curve')
    parser.add_argument('q0', type=int, help='the parameter for the family')
    names = vars(parser.parse_args())
    family = names['family']
    q0 = names['q0']
    if family == 'suzuki':
        c = suzuki(q0)
    elif family == 'hermitian':
        c = hermitian(q0)
    else:
        raise InputError('Family not implemented. Try hermitian or suzuki.')
    c.save(family + str(c.q) + '.curve.json', plain=True)

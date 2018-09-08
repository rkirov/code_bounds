import pickle
from sys import argv

assert len(argv) > 1 

b, w = pickle.load(open(argv[1]+'csstwoptQ'))
b0, w0 = pickle.load(open(argv[1] + 'cssoneptQ'))
b1, w1 = pickle.load(open(argv[1] + 'css+Q'))


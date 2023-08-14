
'''

Definition of functions and variables that are useful to the code in tai.py and
biochem.py.

'''


from functools import reduce
import json


def read_json_file(filename):
    with open(filename, 'r') as f:
        obj = json.load(f)
    return obj

def product(xs):
    return reduce(lambda x,y: x*y,xs)

def geo_mean(xs):
    return (product(xs)) ** (1/len(xs))







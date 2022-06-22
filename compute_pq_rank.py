#!/usr/bin/env sage

import sys
from sage.all import *
from itertools import product
from typing import Tuple, List, Set
from collections import defaultdict

def multiply_monomials(mon1: Tuple[int, ...], mon2: Tuple[int, ...]) -> Tuple[int, ...]:
    """
    Multiplies two monomials of degree two with x^2_i = x_i constraints
    """
    return tuple(sorted(list(set(list(mon1) + list(mon2))))) 

def generate_equations(n: int, p: List[Tuple[int, ...]]) -> Tuple[List[List[Tuple[int, ...]]], List[Tuple[str, Tuple[int, ...]]]]:
    eq_map = defaultdict(list)
    monomials_list = []
    for b in product(range(n), repeat=3):
        if b[0] >= b[1] or b[1] >= b[2]:
            continue
        monomials_list.append(b)
    for b in product(range(n), repeat=2):
        if b[0] >= b[1]:
            continue
        monomials_list.append(b)
    for b in range(n):
        monomials_list.append(tuple([b]))
    for a in p:        
        for b in monomials_list:
            result = multiply_monomials(a, b)
            if (len(result)) > 3:
                eq_map[result].append(('q', b))
            else:
                eq_map[result].append(('q', b))
    for a in eq_map.keys():
        if len(a) <= 3:
            eq_map[a].append(('pq', a))
    return list(eq_map.values())

        
print(generate_equations(5, [(0, 1)]))

A = matrix(GF(2), [[1, 1, 0], [0,1,1], [1,0,1]])
print(A * vector(GF(2), [1,1,1]))
print(matrix(A.kernel().basis()) * vector(GF(2), [1,1,1]))
#print(A.kernel() * vector(GF(2), [1,1,1]).transpose())
A = A.augment(matrix.identity(3), subdivide=True)
#print(A)
print(A.echelon_form())
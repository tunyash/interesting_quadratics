#!/usr/bin/env sage

import sys
from sage.all import *
from itertools import product

import pysat
import pysat.card
import pysat.solvers
from itertools import product, combinations
from typing import Tuple, List, Set
from collections import defaultdict

import random

def upperbound_rank(n: int, p: List[Tuple[int, int]]) -> int:
    pairs_index = [(i,j) for i,j in product(range(n), range(n)) if i < j]
    pairs_map = dict()
    for i, (j, k) in enumerate(pairs_index):
        pairs_map[(j,k)] = i

    equations = []
    pairs_patterns = []
    for i, j in product(range(4), range(4)):
        if i >= j:
            continue
        t = []
        for k in range(4):
            if k not in [i,j]:
                t.append(k)
        pairs_patterns.append(((i, j), tuple(t)))
    for i,j,k,l in product(range(n), repeat=4):
        if i >= j or j >= k or k >= l:
            continue
        indices = [i, j, k, l]
        equation = []
        for ((a, b), (c, d)) in pairs_patterns:
            ra, rb, rc, rd = map(lambda x: indices[x], (a,b,c,d))
            if (ra, rb) in p:
                equation.append((rc, rd))
        equations.append(equation)

    matrix_sparse = dict()
    for i, eq in enumerate(equations):
        for a, b in eq:
            matrix_sparse[i, pairs_map[(a,b)]] = 1

    
    A = matrix(GF(2), len(equations), len(pairs_index), matrix_sparse)
    return len(pairs_index) - A.rank()


def multiply_monomials(mon1: Tuple[int, int], mon2: Tuple[int, int]) -> Tuple[int, ...]:
    """
    Multiplies two monomials of degree two with x^2_i = x_i constraints
    """
    return tuple(sorted(list(set(list(mon1) + list(mon2))))) 

def multiply_polynomials(poly1: List[Tuple[int,int]],
                         poly2: List[Tuple[int,int]]) -> List[Tuple[int, ...]]:
    """
    Multiplies two polynomials over F_2
    """
    result = []
    for a,b in product(poly1, poly2):
        result.append(multiply_monomials(a, b))
    result.sort()
    final = []
    for x in result:
        if not final or final[-1] != x:
            final.append(x)
        elif final[-1] == x:
            final = final[:-1]
    return final

def generate_main_cnf(n: int):
    """
        n: number of possible variables in p and q
    """
    id_pool = pysat.formula.IDPool()
    pairs: List[List[Tuple[int, int]]] = [[], [], []]
    cnf = []
    for poly_id in [1,2]:
        for i, j in product(range(n), repeat=2):
            if i >= j:
                continue
            pairs[poly_id].append((i,j))
            # The variable for the coefficient of x_i * x_j
            id_pool.id((poly_id, i, j))
    # Mapping a monomial of degree <= 4 to list of pairs (monomial in p, monomial in q)
    # that yield the pre-image in multiplication
    product_to_var_pair: Dict[List[int], List[Tuple[int, int]]] = defaultdict(list)
    
    for x, y in product(pairs[1], pairs[2]):
        product_to_var_pair[multiply_monomials(x, y)].append((x, y))
        #print(x, " times ", y, " is ", multiply_monomials(x, y))
        # !(1, x[0], x[1]) or !(2, y[0], y[1]) or (3, x, y)
        cnf.append([-id_pool.id((1, x[0], x[1])),
                    -id_pool.id((2, y[0], y[1])),
                    id_pool.id((3, x, y))])
        cnf.append([-id_pool.id((3, x, y)), id_pool.id((1, x[0], x[1]))])
        cnf.append([-id_pool.id((3, x, y)), id_pool.id((2, y[0], y[1]))])

        
    for result, predecessors in product_to_var_pair.items():
        if len(result) != 4:
            continue
        # Forbid parity 1 for the coefficients of the degree-4 monomials
        for signs in product([0, 1], repeat=len(predecessors)):
            if sum(signs) % 2 == 0:
                continue
            clause = []
            for (x, y), sign in zip(predecessors, signs):
                clause.append((sign * 2 - 1) * id_pool.id((3, x, y)))
            cnf.append(clause)
    return cnf, id_pool

def cnf_with_forced_mathcing(n: int, k :int):
    """
       n: number of variables available in p and q
       k: number of varibles in the matchings
    """
    cnf, id_pool = generate_main_cnf(n)
    matching1 = [(4*i, 4*i+1) for i in range(k)]
    matching2 = [(4*i+2, 4*i+3) for i in range(k)]
    for x, y in matching1:
        cnf.append([id_pool.id((1, x, y))])    
    for x, y in matching2:
        cnf.append([id_pool.id((2, x, y))])
    if k > 0:
        cnf.append([-id_pool.id((2, matching1[0][0], matching1[0][1]))])
        cnf.append([id_pool.id((1, matching1[0][0], matching1[0][1]))])
        
    return cnf, id_pool

def solve(n: int, k :int):
    """
       n: number of variables available in p and q
       k: number of varibles in the matchings
    """
    pairs = [(i, j) for i,j in product(range(n), repeat=2) if i < j]
    cnf, id_pool = cnf_with_forced_mathcing(n, k)
    
    solver = pysat.solvers.Solver(name="m22", bootstrap_with=cnf)
    sat = solver.solve()
    assert sat
    assignment = solver.get_model()
    return [a for a in pairs if id_pool.id((1, *a)) in assignment], \
           [a for a in pairs if id_pool.id((2, *a)) in assignment]

def generate_all(n: int, k: int):
    """
       Enumerates all solutions for the matching-k problem
       n: number of variables available in p and q
       k: number of varibles in the matchings
    """
    pairs = [(i, j) for i,j in product(range(n), repeat=2) if i < j]
    cnf, id_pool = cnf_with_forced_mathcing(n, k)
    
    solver = pysat.solvers.Solver(name="m22", bootstrap_with=cnf)
    sat = solver.solve()
    assert sat
    for assignment in solver.enum_models():
        yield [a for a in pairs if id_pool.id((1, *a)) in assignment], \
              [a for a in pairs if id_pool.id((2, *a)) in assignment]  



def pair_to_letter(p, q, i, j):
    if (i, j) in p and (i, j) in q:
        return 'B'
    elif (i, j) in p:
        return 'P'
    elif (i, j) in q:
        return 'Q'
    else:
        return '.'

def print_tex(n: int, k: int, p: List[Tuple[int, int]], q: List[Tuple[int, int]]):
    # for poly, id in zip([p, q], [1, 2]):
    #     print("p_" + str(id) + " := ")    
    #     for x, y in poly:
    #         print("x_{" + str(x) + "} \cdot x_{" + str(y), end='} + ')
    #     print()        
    
    
    for i in range(n):
        for j in range(n):
            letter = pair_to_letter(p, q, i, j)
            print(letter, end='')
        print()

    # for i in range(n):
    #     for j in range(n):
    #         letter = pair_to_letter(p, q, i, j)
    #         if letter in ['P', 'Q']:
    #             color = 'blue' if letter == 'P' else 'red'
    #             print('\\draw[fill=',color,', color=',color,'](',i,',',j,') circle(1mm);')
    #     print()
    print()

def find_local_patterns(n: int, k: int, p: List[Tuple[int, int]], q: List[Tuple[int, int]]) -> Set[str]:
    local_patterns = set()
    for i1, i2, j1, j2 in product(range(n), repeat=4):
        if i1 >= i2 or i2 >= j1 or j1 >= j2:
            continue
        
        local_patterns.add("".join([pair_to_letter(x, y) for x,y in
                                            [(i1, i2), (i1, j1), (i1, j2), 
                                            (i2, j1), (i2, j2), (j1, j2)]]))
    return local_patterns




exit(0)

for n,k in [(40, 1)]:   
    #try: 
    print(n, k)
    for i, (p, q) in enumerate(generate_all(n, k)):
        print("i=",i)
        print_tex(n, k, p, q)
        print('RANK with fixed p: ', upperbound_rank(n, p))
        print('RANK with fixed q: ', upperbound_rank(n, q))
        
    """for pattern in local_patterns:
        if not ('P' in pattern and 'Q' in pattern):
            continue
        pairs_index = [(i,j) for i, j in product(range(4), repeat=2) if i < j]
        p0 = [(i, j) for k, (i, j) in enumerate(pairs_index) if pattern[k] in ['P', 'B']]
        q0 = [(i, j) for k, (i, j) in enumerate(pairs_index) if pattern[k] in ['Q', 'B']]
        
        print(pattern[:3])
        print('#' + pattern[3:5])
        print('##' + pattern[5])
        print('product=', multiply_polynomials(p0, q0))
        print()"""
    print(len(multiply_polynomials(p, q)))
    print("There are no degree-4 monomials: ", all(lambda x: len(x) <= 3 for x in multiply_polynomials(p, q)))
    print(p)
    #except Exception as e:
    #    print(n, k, " doesn't work")
                    

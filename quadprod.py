import pysat
import pysat.card
import pysat.solvers
from itertools import product, combinations
from typing import Tuple, List
from collections import defaultdict

import random


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


def solve(n: int, k :int):
    """
       n: number of variables available in p and q
       k: number of varibles in the matchings
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

    matching1 = [(4*i, 4*i+1) for i in range(k)]
    matching2 = [(4*i+2, 4*i+3) for i in range(k)]
    for x, y in matching1:
        cnf.append([id_pool.id((1, x, y))])    
    for x, y in matching2:
        cnf.append([id_pool.id((2, x, y))])

    cnf.append([-id_pool.id((2, matching1[0][0], matching1[0][1]))])
    
    solver = pysat.solvers.Solver(name="m22", bootstrap_with=cnf)

    sat = solver.solve()
    assert sat
    assignment = solver.get_model()
    return [a for a in pairs[1] if id_pool.id((1, *a)) in assignment], \
           [a for a in pairs[1] if id_pool.id((2, *a)) in assignment]


for n,k in [(13, 3)]:    
    print(n, k)
    p, q = solve(n, k)

    for poly, id in zip([p, q], [1, 2]):
        print("p_" + str(id) + " := ")    
        for x, y in poly:
            print("x_{" + str(x) + "} \cdot x_{" + str(y), end='} + ')
        print()
    
    for poly in [p, q]:
        for i in range(n):
            for j in range(n):
                if (i, j) in poly:
                    print('#', end='')
                else:
                    print('.', end='')
            print()
        print()

    print("There no degree-4 monomials: ", all(lambda x: len(x) <= 3 for x in multiply_polynomials(p, q)))
    #except Exception as e:
    #    print(str(n), " doesn't work")
                    

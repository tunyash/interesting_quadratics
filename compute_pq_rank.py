#!/usr/bin/env sage

import sys
from sage.all import *
from itertools import product
from typing import Tuple, List, Set, Dict
from collections import defaultdict

def multiply_monomials(mon1: Tuple[int, ...], mon2: Tuple[int, ...]) -> Tuple[int, ...]:
    """
    Multiplies two monomials of degree two with x^2_i = x_i constraints
    """
    return tuple(sorted(list(set(list(mon1) + list(mon2))))) 

def multiply_polynomials(poly1: List[Tuple[int,...]],
                         poly2: List[Tuple[int,...]]) -> List[Tuple[int, ...]]:
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


def num_of_monomials_deg_atmost(n:int, d: int) -> int:
    """
    Calculates the number of monomials over n variables with degree at most d
    """
    num = 0
    for i in range(d + 1):
        num += binomial(n, i)
    return num # without x_i^2=x_i it would have been: (n^(d+1)-1)//(n-1) #sum of 1 + n + n^2 + ... + n^d


# def monomial_to_position(mon: Tuple[int, ...]) -> int:
#     """
#     Calculates the index (for either row or column) corresponding to the input monomial
#     """
#     deg = len(mon)
#     if deg == 0:
#         return 0
#     num_of_smaller_deg = num_of_monomials_deg_atmost(deg - 1)
#     mon.sort()
#     num_of_smaller_vars = 0
#     for i in range(deg):
#         num_of_smaller_vars += binomial(mon[i] - 1, deg - i) # add the number of monomials that agree with mon up to the (i-1)'th variable, but have a smaller i'th variable
#     return  num_of_smaller_deg + num_of_smaller_vars


# def position_to_monimial(i: int) -> tuple:
#     """
#     Converts the index (for either row or column) i into the monomial represented by this position
#     """
#     same_deg_first_i = 0
#     deg = 0
#     next_deg_first_i = 1
#     while next_deg_first_i < i:
#         same_deg_first_i = next_deg_first_i
#         deg += 1
#         next_deg_first_i = num_of_monomials_deg_atmost(deg - 1)
#     #not finished
#     return #TODO

mon_to_i: Dict[Tuple[int, ...], int]
i_to_mon: Dict[int, Tuple[int, ...]]


def create_multiplying_matrix(n: int, eq_map: Dict[tuple, List[tuple]]):
    """
    Creates the matrix A that maps a vector representing a polynomial p to the vector representing the product pq, given the eq_map which lists, for each monomial r in p, the monomials l*r for the all possible monomials l. 
    """
    matrix_dict = {}
    prj_4_dict = {}
    for i, mon in i_to_mon.items():
        prj_4_dict[(i,i)] = 1 if len(mon) >= 4 else 0
    for mon, eq in eq_map.items():
        i = mon_to_i[mon]
        for term in eq:
            if term[0]=='q':
                assert (len(term[1]) <= 3)
                matrix_dict[(i, mon_to_i[term[1]])] = 1
        
    return matrix(GF(2), len(mon_to_i), num_of_monomials_deg_atmost(n, 3), matrix_dict),\
           matrix(GF(2), len(mon_to_i), len(mon_to_i), prj_4_dict)



def generate_equations(n: int, p: List[Tuple[int, ...]]) -> Tuple[List[List[Tuple[int, ...]]], List[Tuple[str, Tuple[int, ...]]]]:
    eq_map = defaultdict(list)
    monomials_list = []
    for b in product(range(n), repeat=3): # cubic monomials
        if b[0] >= b[1] or b[1] >= b[2]:
            continue
        monomials_list.append(b)
    for b in product(range(n), repeat=2): # quadratic monomials
        if b[0] >= b[1]:
            continue
        monomials_list.append(b)
    for b in range(n): # linear monomials
        monomials_list.append(tuple([b]))
    monomials_list.append(()) # the constant monomial
    global mon_to_i
    mon_to_i = dict(zip(monomials_list, range(len(monomials_list))))
    global i_to_mon
    i_to_mon = dict(zip(range(len(monomials_list)), monomials_list))

    for a in p:        
        for b in monomials_list:
            result = multiply_monomials(a, b)
            eq_map[result].append(('q', b))
    for a in eq_map.keys():
        if a not in mon_to_i:
            new_i = len(mon_to_i)
            mon_to_i[a] = new_i
            i_to_mon[new_i] = a
        if len(a) <= 3:
            eq_map[a].append(('pq', a))
    A, prj = create_multiplying_matrix(n, eq_map)
  #  print(A)
   # print(prj)
    # print("A: rows="+str(A.nrows())+" cols="+str(A.ncols()))
    # print("prj: rows="+str(prj.nrows())+" cols="+str(prj.ncols()))
    return A, prj


def calculate_dim_of_prod(n: int, A, prj) -> int:
    return (A * matrix((prj * A).transpose().kernel().basis()).transpose()).rank()

def calculate_dim_of_qs(n: int, A, prj) -> int:
    return matrix((prj * A).transpose().kernel().basis()).rank()

for n in range(11, 40):
    p = [(0, 1, 2), (4, 5, 6), (1, 2), (4,5), (6, 7)]
    A, prj = generate_equations(n, p)
    print('n=', n, 'rank of prod=', calculate_dim_of_prod(n, A, prj), "  rank of qs=", calculate_dim_of_qs(n, A, prj))
    print('p=', p)
    vec = matrix((prj * A).transpose().kernel().basis())[0]
    q = []
    for j in range(len(vec)):
        if vec[j]==1:
            q.append(i_to_mon[j])
    
    print('q in kernel=', q, multiply_polynomials(p, q))

# A = matrix(G3(2), [[1, 1, 0], [0,1,1], [1,0,1]])
# print(A * vector(GF(2), [1,1,1]))
# print(matrix(A.kernel().basis()) * vector(GF(2), [1,1,1]))
# #print(A.kernel() * vector(GF(2), [1,1,1]).transpose())
# A = A.augment(matrix.identity(3), subdivide=True)
# #print(A)
# print(A.echelon_form())

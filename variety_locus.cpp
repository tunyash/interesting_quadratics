#include "bits/stdc++.h"

using namespace std;

struct Monomial {
  int _x, _y;
  Monomial(int x, int y): _x(x), _y(y) {
      if (_x > _y) {
          swap(_x, _y);
      }
  }
  bool intersects(const Monomial& other) {
      return _x == other._x || _x == other._y || _y == other._x || _y == other._y;
  }
  
};

bool operator <(const Monomial& a, const Monomial& b) {
    return a._x < b._x || (a._x == b._x && a._y < b._y);
}

bool operator ==(const Monomial& a, const Monomial& b) {
    return a._x == b._x && a._y == b._y;
}


typedef vector<Monomial> Polynomial;

ostream& operator <<(ostream & ost, const Polynomial& a) {
    for (auto mon: a) {
        ost << mon._x << '*' << mon._y << '+';
    }
    return ost;
}

typedef long long int64;

int evaluate(const Polynomial& a, vector<bool> point) {
    int value = 0;
    for (auto mon : a) {
        value ^= point[mon._x] & point[mon._y]; 
    }
    return value;
}

#define varindex(i,j,r) i*r+j

vector<Polynomial> generate_polynomials(int r) {
    vector<vector<bool>> subsets;
    for (int bit1 = 0; (1<<bit1) < r; ++bit1) {
        for (int bit2 = bit1 + 1; (1<<bit2) < r; ++bit2) {
            for (int bit3 = bit2 + 1; (1<<bit3) < r; ++bit3) {
                for (int assignment = 0; assignment < 8; ++assignment) {
                    vector<bool> subset(r, false);
                    for (int i = 0; i < r; ++i) {
                        int v1 = (i>>bit1) & 1;
                        int v2 = (i>>bit2) & 1;
                        int v3 = (i>>bit3) & 1;
                        if (v1 * 4 + v2 * 2 + v3 == assignment) {
                            subset[i] = 1;
                        }
                    }
                    subsets.push_back(subset);
                }
            }
        }
    }
    //subsets.push_back(vector<bool>(r, true)); // TODO: add other subsets later
    vector<Polynomial> result;
    for (int i1 = 0; i1 < r; ++i1) {
        for (int i2 = 0; i2 < r; ++i2) {
            if (i1 == i2) continue;
            for (int a = 0; a < r; ++a) {
                for (auto subset : subsets) { 
                    Polynomial p = vector<Monomial>();
                    for (int j = 0; j < r; j++) {
                        if (subset[j]) p.push_back(Monomial(varindex(i1, j, r), varindex(i2, (j+a) % r, r)));                        
                    }
                    result.push_back(p);
                }
                for (auto subset : subsets) { // different subsets!!!! 
                    Polynomial p = vector<Monomial>();
                    for (int j = 0; j < r; j++) {
                        if (subset[j]) p.push_back(Monomial(varindex(j, i1, r), varindex((j+a) % r, i2, r)));
                    }
                    result.push_back(p);
                }                
            }
        }
    }
    // Adding quadratics
    for (int i1 = 0; i1 < r; i1++)
       for (int j = 0; j < r; ++j)
         for (int k = j + 1; k < r; ++k) {
            Polynomial p = {Monomial(varindex(i1, j, r), varindex(i1, k, r))};
            result.push_back(p);
            p = {Monomial(varindex(j, i1, r), varindex(k, i1, r))};
            result.push_back(p);
            
         }
    for (auto poly : result) {
        sort(poly.begin(), poly.end());
    }
    return result;
}

bool check_polynomial_construction(int n, const vector<Polynomial>& polys) { 
    set<tuple<int,int,int>> triplets;

    for (auto poly : polys) {
        sort(poly.begin(), poly.end());
        cout << poly << endl;
    }
    for (auto poly1 : polys) {
        for (auto poly2: polys) {
            if (poly1 == poly2) continue;
            int counter = 0;
            int a = -1, b = -1, c = -1;
            for (auto mon1 : poly1) {
                for (auto mon2: poly2) {
                    if (mon1.intersects(mon2)) {
                        counter++;
                        if (counter > 1) {
                            break;
                        }
                        a = mon1._x;
                        b = mon1._y;
                        if (mon2._x != a && mon2._x != b) c = mon2._x;
                        if (mon2._y != a && mon2._y != b) c = mon2._y;                         
                    }
                }
            }
            if (a != -1 && b != -1 && c != -1 && counter == 1 && a != b && a != c && b != c) {
                if (a > b) swap(a,b);
                if (b > c) swap(b,c);
                if (a > b) swap(a,b);
                triplets.insert(make_tuple(a,b,c));
            }
        }
    }
    cout << triplets.size() << endl;
    return (int)triplets.size() == n * (n-1) * (n-2) / 6;
}

pair<bool, vector<bool>> find_nontrivial_point(int n, int k, const vector<Polynomial>& polys, vector<bool> current_vector, int i) {
  if (k < 0) {
      return make_pair(false, vector<bool>());
  }
  if (i == n) {
      if (k != 0) {
          return make_pair(false, vector<bool>());
      }
      for (auto poly: polys) {
          if (evaluate(poly, current_vector) == 1) {
              return make_pair(false, vector<bool>());
          }
      }
      return make_pair(true, current_vector);
  }
  
  current_vector[i] = true;
  auto attempt1 = find_nontrivial_point(n, k - 1, polys, current_vector, i+1);
  if (attempt1.first) return attempt1;
  current_vector[i] = false;
  return find_nontrivial_point(n, k, polys, current_vector, i+1);
}

int main() {
   int r;
   cin >> r;
   int n = r * r;
   auto polys = generate_polynomials(r);
   vector<bool> current_vector(n, false);
   auto result = find_nontrivial_point(n, 4, polys, current_vector, 0);
   cout << (result.first ? "YES" : "NO") << endl;
   if (result.first) {
        for (int i = 0; i < n; ++i) {
            cout << result.second[i];
        }
}
   cout << endl;
}
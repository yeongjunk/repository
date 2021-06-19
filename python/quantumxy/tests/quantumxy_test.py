import quantumxy
from scipy import sparse
from functools import reduce


def test():
    test1 = (quantumxy.ham_xy() == quantumxy.ham_xy_kron()).toarray().all()
    vals, vecs = LA.eigh(quantumxy.ham_xy(h = -1.5).toarray())
    test2 = (vecs[:,0] == [0, 0, 0, 0, 0, 0, 0, 1]).all()
    vals, vecs = LA.eigh(quantumxy.ham_xy(h = 1.5).toarray())
    test3 = (vecs[:,0] == [1, 0, 0, 0, 0, 0, 0, 0]).all()

    print("hamiltonian_kronecker :",test1)
    print("h > 0 :", test2)
    print("h < 0 : ", test3)

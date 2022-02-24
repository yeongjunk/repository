import numpy as np
from scipy import sparse
from functools import reduce

s_p = sparse.csc_matrix([[0., 1.], [0., 0.]], dtype=np.double)
s_m = s_p.transpose()
s_z = sparse.csc_matrix([[1., 0.], [0., -1.]], dtype=np.double)

def bin_arr(value, width):
    string = format(value, '0{}b'.format(width))
    binary_list = np.array([0 if c == '0' else 1 for c in string], dtype = np.bool)
    return binary_list

def basis_vectors(N = 3):
    basis = np.array([bin_arr(i, N) for i in range(2**N-1, -1, -1)])
    return basis

def ham_xy_kron(N = 3, J = 1., h = 1.):
    H = sparse.csc_matrix((2**N, 2**N), dtype=np.double)
    H_offdiag = sparse.csc_matrix((2**N, 2**N), dtype=np.double)

    for i in range(0, N-1):
        I1 = sparse.identity(2**i, dtype = np.double)
        I2 = sparse.identity(2**(N-i-2), dtype = np.double)
        I3 = sparse.identity(2**(N-i-1), dtype = np.double)

        H_offdiag -= J*reduce(sparse.kron, [I1, s_p, s_m, I2])
        H -= h*reduce(sparse.kron, [I1, s_z, I3]) # Diagonal term

    i = N-1
    H -= h*sparse.kron(sparse.identity(2**i, dtype = np.double), s_z)

    H_offdiag -= J*reduce(sparse.kron, [s_p, sparse.identity(2**(N-2), dtype = np.double), s_m])
    H_offdiag += H_offdiag.transpose() # hermitian conjugate
    H += H_offdiag

    return H

def ham_xy(N = 3, J = 1., h = 1.):
    v = basis_vectors(N = N)
    down_up = np.array([False, True])

    row = np.array([], dtype = np.int)
    col = np.array([], dtype = np.int)
    val = np.array([], dtype = np.double)

    for i in range(0, len(v)):
        u = np.copy(v[i])
        onsite = 0.
        for ui in u:
            onsite += -h if ui else h
        row = np.append(row, i)
        col = np.append(col, i)
        val = np.append(val, onsite)


        for j in range(0, N-1):
            if (u[j:j+2] == down_up).all():
                i_new = i - 2**(N-j-2)
                row = np.append(row, np.array([i_new, i]))
                col = np.append(col, np.array([i, i_new]))
                val = np.append(val, -J*np.ones(2))

        j = N-1
        if (u[[j,0]] == down_up).all():
            i_new = i + 2**(N-1) - 1
            row = np.append(row, np.array([i_new, i]))
            col = np.append(col, np.array([i, i_new]))
            val = np.append(val, -J*np.ones(2))

    return sparse.csc_matrix((val, (row, col)), shape=(2**N, 2**N), dtype = np.double)

def basis_magnetizations(N = 3):
    basis = basis_vectors(N = N)
    m = np.zeros(2**N, dtype = np.double)
    for (i, basis_i) in enumerate(basis):
        for s in basis_i:
            m[i] += (1 if s else -1)
    return m/N

def basis_correlations(N = 3, r = 1):
    basis = basis_vectors(N = N)
    corr = np.zeros(2**N, dtype = np.double)
    m = np.zeros(2**N, dtype = np.double)
    for (i, basis_i) in enumerate(basis):
        for s in basis_i:
            m[i] += (1. if s else -1.)

        for n in range(0, N):
            if not(basis_i[n]^basis_i[(n+r)%N]):
                corr[i] += 1
            else:
                corr[i] -= 1

    return corr/N - (m/N)**2

def compute_magnetization(v):
    N = np.int(np.log2(len(v)))
    m_basis = basis_magnetizations(N = N)
    m = 0
    for i in range(0, len(v)):
        m += m_basis[i]*(v[i]**2)
    return m

def compute_correlation(v, r):
    N = np.int(np.log2(len(v)))
    corr_basis = basis_correlations(N = N, r = r)
    corr = 0
    for i in range(0, len(v)):
        corr += corr_basis[i]*(v[i]**2)
    return corr

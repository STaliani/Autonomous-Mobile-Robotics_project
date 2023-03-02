from cdd import Matrix, Polyhedron, RepType
from numpy import array, hstack, zeros

def face_of_span(S):
    """

    Returns the face matrix S^F of the span matrix S,
    that is, a matrix such that

        {x = S z, z >= 0} if and only if {S^F x <= 0}.

    """
    V = hstack([zeros((S.shape[1], 1)), S.T])

    V_cdd = Matrix(V, number_type='float')
    V_cdd.rep_type = RepType.GENERATOR
    P = Polyhedron(V_cdd)
    ineq = P.get_inequalities()
    H = array(ineq)
    if H.shape == (0,):  # H = []
        return H
    # b, A = H[:, 0], -H[:, 1:]  # H matrix is [b, -A]
    A = []
    for i in range(H.shape[0]):
        if H[i, 0] != 0:  # b should be zero
            raise ValueError('The span is not a polytope.')
        elif i not in ineq.lin_set:
            A.append(-H[i, 1:])
    return array(A)
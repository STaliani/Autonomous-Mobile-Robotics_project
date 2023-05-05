from cdd import Matrix, Polyhedron, RepType
from numpy import array, hstack, zeros


def span_of_face(F):
    """

    Compute the span matrix F^S of the face matrix F,
    that is, a matrix such that

        {F x <= 0} if and only if {x = F^S z, z >= 0}.

    """
    b, A = zeros((F.shape[0], 1)), F
    # H-representation: b - A x >= 0
    # ftp://ftp.ifor.math.ethz.ch/pub/fukuda/cdd/cddlibman/node3.html
    # the input to pycddlib is [b, -A]
    F_cdd = Matrix(hstack([b, -A]), number_type='float')
    F_cdd.rep_type = RepType.INEQUALITY
    P = Polyhedron(F_cdd)
    g = P.get_generators()
    V = array(g)
    rays = []
    for i in range(V.shape[0]):
        if V[i, 0] != 0:  # 1 = vertex, 0 = ray
            raise ValueError('The face is not a polytope.')
        elif i not in g.lin_set:
            rays.append(V[i, 1:])
    return array(rays).T

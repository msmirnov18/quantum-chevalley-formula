reset()

# Dynkin type
D = "E6"

# which P_K we quotient out (need not be maximal parabolic)
# maximal parabolic correspond to K's with just one element
K = [2]

print('D =', D, ',', 'K =', K)

G = CartanType(D)

RG = RootSystem(G).root_lattice()

W = RG.weyl_group(prefix="s")
def s(i):
    assert i in G.index_set(), 'the index of S is not allowed'
    return W.gen(i-1)

# # defining minimal length coset representatives (they index Schubert classes)
# WP = []
# for w in W:
#     v = w.coset_representative([j for j in G.index_set() if j not in K])
#     if v not in WP: WP.append(v)
#
# # sorting coset representatives by length
# WP_sorted = []
# for i in range(W.long_element().length()):
#     for v in WP:
#         if v.length() == i: WP_sorted.append(v)

# For complicated examples we load WP and WP_sorted explicitly
# Their definitions above need to be commented out
# load cosetreps-e7-p1.sage
load cosetreps-e6-p2.sage
# load cosetreps-e8-p8.sage

# defining the quantum parameters, i.e. the Novikov ring
# they are indexed by the elements of the index set K that defines the parabolic
qnames=['q%d'%(i) for i in K]
N = PolynomialRing(QQ,qnames)
def q(i):
    assert i in K, 'the index of the q-variable is not allowed'
    return N.gen(K.index(i))

# defining the small quantum cohomology module as a free module over the Novikov ring
# with a basis given by the Schubert classes
QH = CombinatorialFreeModule(N, WP)

# defining a function that takes as input an element w in WP and returns
# the corresponding basis element of QH
def S(w):
    return QH.monomial(w)

# converting a root alpha into the correponding element s_alpha of W
def root_to_reflection(alpha):
    assert alpha in RG.roots(), 'input is not a root'
    return product([W.gen(j-1) for j in list(alpha.associated_reflection())])

# maybe this function needs a better name!
# It is defined in Lemma 3.5 of [FuWo]
def nn(alpha):
    assert alpha in RG.roots(), 'input is not a root'
    nonparabolic_positive_roots = RG.nonparabolic_positive_roots(tuple(j for j in G.index_set() if j not in K))

    return 2*(sum(nonparabolic_positive_roots).to_ambient()).inner_product(alpha.to_ambient())/(alpha.to_ambient()).inner_product(alpha.to_ambient())

# maybe this function needs a better name!
# It is defined right before Lemma 3.1 of [FuWo] and denoted d(\alpha).
# More precisely, this function returns h_{\alpha}(\omega_{\beta}),
# which is a coefficient apperaring in the defintion of d(\alpha) in [FuWo].
# Input:
# alpha can be any root;
# l is the index of the simple root, which is denoted beta in [FuWo].
def dd(alpha,l):
    assert alpha in RG.roots(), 'alpha is not a root'
    assert l in K, 'value of l not allowed'

    beta = RG.simple_root(l)

    return alpha.coefficient(l)*(beta.to_ambient()).inner_product(beta.to_ambient())/(alpha.to_ambient()).inner_product(alpha.to_ambient())

# This is the Quantum Chevalley Formula as in Theorem 10.1 of [FuWo].
# Input:
# l is the index of the simple root giving us a divisor class s_l (corresponds to \sigma_{s_{\beta}} in [FuWo]);
# input is an element in the quantum cohomology module QH
def quantum_chevalley(l,input):

    assert l in K, 'value of l not allowed'
    assert input in QH, 'input is not in QH'

    if len(input.monomials()) == 0:
        return QH.zero()
    if len(input.monomials()) > 1:
        u = next(iter(input.monomial_coefficients()))
        coefficient_of_u = input.coefficient(u)
        first_summand = coefficient_of_u * S(u)
        return quantum_chevalley(l,first_summand) + quantum_chevalley(l, input - first_summand)

    w = next(iter(input.monomial_coefficients()))
    commonfactor = input.coefficient(w)

    nonparabolic_positive_roots = RG.nonparabolic_positive_roots(tuple(j for j in G.index_set() if j not in K))

    result = QH.zero()

    for alpha in nonparabolic_positive_roots:
        v = (w*root_to_reflection(alpha)).coset_representative([j for j in G.index_set() if j not in K])
        if v.length() == w.length() + 1:
           result = result + dd(alpha,l) * S(v)
        if v.length() == w.length() + 1 - nn(alpha):
           result = result + product([q(j)^(dd(alpha,j)) for j in K]) * dd(alpha,l) * S(v)

    return commonfactor*result

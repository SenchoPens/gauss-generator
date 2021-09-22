import sympy

from sympy import sin, cos, Symbol


NEW_LINE = True
PRINT = None
ROWS = 0
COLUMNS = 1


class Matrix(sympy.Matrix):
    def to_latex(self, det=False):
        mtype = ("v" if det else "p") + "matrix"
        l = '\\\\ \n'.join(' & '.join(map(sympy.latex, self.row(i))) for i in range(self.rows))
        return lambda x : f'\\begin{{{mtype}}}\n' + l + f'\n\\end{{{mtype}}}'


class GaussMatrix(Matrix):
    def __init__(self, A, n_of_free=0):
        # if isinstance(A, GaussMatrix):
        if hasattr(A, 'n_of_free'):
            self.n_of_free = A.n_of_free
        else:
            self.n_of_free = n_of_free

    def insertvrule(self, row):
        if self.n_of_free != 0:
            row = row[:-self.n_of_free] + [r'\vrule'] + row[-self.n_of_free:]
        return row

    def to_latex(self):
        l = '\\\\ \n'.join(
            ' & '.join(
                map(
                    lambda x:
                        x if isinstance(x, str)
                        else sympy.latex(x),
                    self.insertvrule(self.row(i))
                )
            )
            for i in range(self.rows)
        )
        def wrap(rowops='', colops=''):
            return (
                '\\begin{gmatrix}[p]\n' + l + '\n'
                + (r'\rowops' '\n' + rowops if rowops else '')
                + (r'\colops' '\n' + colops if colops else '')
                + r'\end{gmatrix}' + '\n' + r'\leadsto'
            )
        return wrap


def eps1_r(A, j, i, k):
    if i == j:
        raise ValueError()
    rowops = r'\add' + f'[{k}]{{{i}}}{{{j}}}'
    for t in range(len(A[j, :])):
        A[j, t] += k * A[i, t]
    return (A, rowops)


def eps1_c(A, j, i, k):
    if i == j:
        raise ValueError()
    colops = r'\add' + f'[{k}]{{{i}}}{{{j}}}'
    for t in range(len(A[j, :])):
        A[t, j] += k * A[t, i]
    return (A, colops)


def eps2_r(A, i, j):
    if i == j:
        raise ValueError()
    rowops = r'\swap' + f'{{{i}}}{{{j}}}'
    A[i, :], A[j, :] = A[j, :], A[i, :]
    return (A, rowops)


def eps2_c(A, i, j):
    if i == j:
        raise ValueError()
    colops = r'\swap' + f'{{{i}}}{{{j}}}'
    A[:, i], A[:, j] = A[:, j], A[:, i]
    return (A, colops)


def eps3_r(A, i, k, div=False, to_int=False):
    if k == 0:
        raise ValueError()
    op = ':' if div else r'\cdot'
    rowops = r'\mult' + f'{{{i}}}{{{op} {k}}}'
    for t in range(len(A[i, :])):
        if div:
            A[i, t] /= k
        else:
            A[i, t] *= k
        if to_int:
            A[i, t] = int(A[i, t])
    return (A, rowops)


def eps3_c(A, i, k, div=False, to_int=False):
    if k == 0:
        raise ValueError()
    op = ':' if div else r'\cdot'
    colops = r'\mult' + f'{{{i}}}{{{op} {k}}}'
    for t in range(len(A[:, i])):
        if div:
            A[t, i] /= k
        else:
            A[t, i] *= k
        if to_int:
            A[t, i] = int(A[t, i])
    return (A, colops)


def gmap_r(A, i, f):
    for t in range(len(A[i, :])):
        try:
            A[i, t] = f(A[i, t])
        except:
            pass
    return (A, '')


def gmap_c(A, i, f):
    for t in range(len(A[:, i])):
        try:
            A[t, i] = f(A[t, i])
        except:
            pass
    return (A, '')


eps1 = (eps1_r, eps1_c)
eps2 = (eps2_r, eps2_c)
eps3 = (eps3_r, eps3_c)
gmap = (gmap_r, gmap_c)


def apply_ops(A, args, i):
    allops = ([], [])
    for t in args:
        if len(t) == 2:
            f, a = t
            A, ops = f[i](A, *a)
        else:
            f, a, kw = t
            A, ops = f[i](A, *a, **kw)
        if ops:
            allops[i].append(ops)
    return A, (('\n'.join(ops) + '\n' if ops else '') for ops in allops)


def mod_print(A, ops, print_=None, i=0):
    if print_ is None:
        print_ = print

    l = A.to_latex()
    A, (ops_r, ops_c) = apply_ops(A, ops, i)
    print_(l(ops_r, ops_c).replace('*', r' \cdot '))
    return A


def mods_print(A, *args, print_=None):
    if print_ is None:
        print_ = print

    i = 0
    cur = []
    cur_orientation = 0
    while i < len(args):
        if args[i] is PRINT:
            A = mod_print(A, cur, print_, i=cur_orientation)
            cur = []
            i += 1
            continue

        if args[i] == COLUMNS or args[i] == ROWS:
            cur_orientation = args[i]
            cur = []
            i += 1
            continue

        if args[i] is NEW_LINE:
            print_(r'$$ $$ \leadsto')
            i += 1
            continue

        if i + 2 < len(args) and isinstance(args[i + 2], dict):
            cur.append((args[i], args[i + 1], args[i + 2]))
            i += 3
        else:
            cur.append((args[i], args[i + 1]))
            i += 2

    if cur:
        A = mod_print(A, cur, print_, i=cur_orientation)
    mod_print(A, [], print_, i=cur_orientation)
    return A

# TODO: https://tex.stackexchange.com/questions/370323/how-to-use-gauss-package-which-displays-as-below


def matan_1_1():
    x = Symbol('x')
    q = Symbol('q')

    # A = sympy.Matrix([
    #     [1 - q * cx, q * sx, q * cx],
    #     [-q * sx, 1 - q * cx, q * sx],
    # ])
    # rref = A.rref()[0]
    # rref[0, 2] = rref[0, 2].simplify()
    # rref[1, 2] = rref[1, 2].simplify()
    # print(sympy.latex(rref))
    # print(A.rref())
    
    A = GaussMatrix([
        [1 - q * cos(x), q * sin(x), q * cos(x)],
        [-q * sin(x), 1 - q * cos(x), q * sin(x)],
    ], n_of_free=1)
    A = mods_print(A,
        ROWS,
        eps3, (0, sin(x)),
        gmap, (0, sympy.expand),
        eps3, (1, cos(x)),
        gmap, (1, sympy.expand),
        PRINT,
        eps1, (0, 1, -1),
        gmap, (0, sympy.simplify),
        eps1, (1, 0, q * cos(x)),
        gmap, (1, sympy.simplify),
        PRINT,
        eps3, (1, (q - cos(x)) * cos(x)), {"div": True},
        gmap, (1, sympy.simplify),
    )


matan_1_1()

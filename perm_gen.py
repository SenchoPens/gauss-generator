import sympy

def to_latex(A):
    l = '\\\\ \n'.join(' & '.join(map(sympy.latex, row)) for row in A)
    return '\\begin{pmatrix}\n' + l + '\n\\end{pmatrix}'


class Permutation:
    def __init__(self, args):
        if isinstance(args, list):
            self.p = {(i + 1): x for (i, x) in enumerate(args)}
        else:
            self.p = args

    def __getitem__(self, key):
        return self.p.get(key, key)

    # def __setitem__(self, key, value):
    #     self.p[key] = value

    def __mul__(self, other):
        new_p = {}
        for x in self.p.keys() | other.p.keys():
            new_p[x] = self[other[x]]
        return Permutation(new_p)

    def factor(self):
        xs = set(self.p.keys())
        cs = []
        while xs:
            cycle = [min(xs)]
            xs.remove(cycle[-1])
            while self[cycle[-1]] != cycle[0]:
                cycle.append(self[cycle[-1]])
                xs.remove(cycle[-1])
            cs.append(Cycle(cycle))
        return cs

    def __str__(self):
        vals = list(range(1, max(self.p.keys()) + 1))
        return to_latex([vals, [self[x] for x in vals]])

    def __pow__(self, p):
        if p == 0:
            return Permutation(dict())
        if p == 1:
            return self
        if p == -1:
            # vals = list(range(1, max(max(self.p.keys()), max(self.p.values())) + 1))
            return Permutation({y: x for (x, y) in self.p.items()})
        if p < 0:
            return (self ** -1) ** (-p)
        return (self ** (p // 2)) * (self ** (p - p // 2))


def str_factorization(fs):
    return '\n'.join(str(c) for c in fs)


class Cycle(Permutation):
    def __init__(self, args):
        self.cycle = args[:]
        args.append(args[0])
        super().__init__({x: y for (x, y) in zip(args, args[1:])})

    def __str__(self):
        return to_latex([self.cycle])


def main():
    # print(to_latex([[1, 1,], [2, 2]]))
    r = Permutation([4, 1, 8, 6, 7, 3, 5, 2])
    s = Permutation([5, 6, 2, 8, 1, 7, 4, 3])
    t = Permutation([5, 1, 3, 8, 4, 7, 6, 2])
    print(r)
    print(s)
    print(t)
    print('---')
    print(str_factorization(r.factor()))
    print('---')
    print(s ** -1)
    print('---')
    print(str_factorization((r * s ** -1).factor()))
    print('---')
    print((r * s ** -1) ** -173)
    print('---')
    print((r * s ** -1) ** -173 * t)


main()

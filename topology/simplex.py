"""
simplex.py (2026-10-19) - directed flag complex simplex counts (Reimann et al. 2017) + degree-preserving null.
Pure Python (pyflagser has no wheel for Python 3.14). A directed k-simplex = ordered (v0..vk) with an edge vi->vj for
every i<j (a transitive tournament); reciprocal pairs give two 1-simplices, as in flagser.
  count_simplices(out_bits)          -> [n_0, n_1, ...] (n_d = number of d-simplices)
  vertex_max_dim(out_bits)           -> highest simplex dimension each vertex belongs to
  degree_preserving_null(edges, ...) -> edge list with identical in- and out-degree of every node
  python simplex.py selftest         -> brute-force comparison on random graphs
"""
import sys, itertools, random


def to_bits(n, edges):
    out = [0] * n
    for a, b in edges:
        if a != b:
            out[a] |= 1 << b
    return out


def _walk(out, visit, max_dim):
    # each simplex has a unique vertex order, so a DFS that extends (v0..vk) by any w in out(v0)&...&out(vk) visits it once
    stack = [((v,), out[v]) for v in range(len(out))]
    while stack:
        simp, common = stack.pop()
        visit(simp)
        if len(simp) > max_dim:
            continue
        c = common
        while c:
            low = c & -c
            u = low.bit_length() - 1
            c ^= low
            stack.append((simp + (u,), common & out[u]))


def count_simplices(out, max_dim=60):
    counts = []

    def visit(s):
        d = len(s) - 1
        if d == len(counts):
            counts.append(0)
        counts[d] += 1
    _walk(out, visit, max_dim)
    return counts


def vertex_max_dim(out, max_dim=60):
    best = [0] * len(out)

    def visit(s):
        d = len(s) - 1
        for v in s:
            if d > best[v]:
                best[v] = d
    _walk(out, visit, max_dim)
    return best


def degree_preserving_null(edges, swaps_per_edge=10, seed=0):
    """Directed double-edge swaps a->b, c->d  =>  a->d, c->b (no self-loops, no duplicate edges)."""
    rng = random.Random(seed)
    E = list(edges)
    S = set(E)
    m = len(E)
    done, tries = 0, 0
    while done < swaps_per_edge * m and tries < 50 * swaps_per_edge * m:
        tries += 1
        i, j = rng.randrange(m), rng.randrange(m)
        (a, b), (c, d) = E[i], E[j]
        if a == d or c == b or a == c or b == d or (a, d) in S or (c, b) in S:
            continue
        S.difference_update([(a, b), (c, d)])
        S.update([(a, d), (c, b)])
        E[i], E[j] = (a, d), (c, b)
        done += 1
    return E, done


def brute_counts(n, edges):
    es = set(edges)
    counts = [n]
    for k in range(2, n + 1):
        c = sum(1 for p in itertools.permutations(range(n), k)
                if all((p[i], p[j]) in es for i in range(k) for j in range(i + 1, k)))
        if not c:
            break
        counts.append(c)
    return counts


def selftest():
    rng = random.Random(1)
    for t in range(200):
        n = rng.randint(2, 7)
        p = rng.random()
        edges = [(a, b) for a in range(n) for b in range(n) if a != b and rng.random() < p]
        got, exp = count_simplices(to_bits(n, edges)), brute_counts(n, edges)
        assert got == exp, (t, n, edges, got, exp)
    # null keeps every in/out degree and the edge count
    edges = [(a, b) for a in range(60) for b in range(60) if a != b and rng.random() < 0.1]
    E2, done = degree_preserving_null(edges, 10, seed=3)
    io = lambda E: (sorted((a, sum(1 for x, _ in E if x == a)) for a in range(60)), sorted((b, sum(1 for _, y in E if y == b)) for b in range(60)))
    assert io(edges) == io(E2) and len(set(E2)) == len(E2) == len(edges) and all(a != b for a, b in E2)
    # a transitive tournament on 6 nodes = one 5-simplex: counts C(6,k+1)
    tt = [(a, b) for a in range(6) for b in range(a + 1, 6)]
    assert count_simplices(to_bits(6, tt)) == [6, 15, 20, 15, 6, 1]
    assert vertex_max_dim(to_bits(6, tt)) == [5] * 6
    print("selftest ok: 200 random graphs match brute force; null preserves degrees (", done, "swaps ); tournament ok")


if __name__ == "__main__":
    selftest() if sys.argv[1:] == ["selftest"] else print(__doc__)

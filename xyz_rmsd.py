#!/usr/bin/env python3
"""xyz_rmsd.py - 比较两个 xyz 结构的 RMSD（先做刚体对齐）。

用法：
    python xyz_rmsd.py a.xyz b.xyz
    python xyz_rmsd.py a.xyz b.xyz --no-align     # 不做对齐，直接算

只依赖标准库。对齐用 Kabsch 算法，3x3 矩阵的 SVD 用 Jacobi 迭代手写实现，
避免 numpy 依赖。
"""
import argparse
import math
import sys


def read_xyz(path):
    with open(path, "r") as f:
        lines = [ln.rstrip("\n") for ln in f if ln.strip()]
    n = int(lines[0].split()[0])
    atoms = []
    for ln in lines[2:2 + n]:
        parts = ln.split()
        if len(parts) < 4:
            raise ValueError(f"{path}: 格式不对: {ln!r}")
        atoms.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
    if len(atoms) != n:
        raise ValueError(f"{path}: 声明 {n} 个原子，实际读到 {len(atoms)}")
    return atoms


def centroid(atoms):
    n = len(atoms)
    return (sum(a[1] for a in atoms) / n,
            sum(a[2] for a in atoms) / n,
            sum(a[3] for a in atoms) / n)


def centered(atoms, c):
    return [(a[0], a[1] - c[0], a[2] - c[1], a[3] - c[2]) for a in atoms]


def rmsd_direct(a, b):
    if len(a) != len(b):
        raise ValueError("原子数不同")
    s = 0.0
    for pa, pb in zip(a, b):
        dx = pa[1] - pb[1]
        dy = pa[2] - pb[2]
        dz = pa[3] - pb[3]
        s += dx * dx + dy * dy + dz * dz
    return math.sqrt(s / len(a))


def jacobi_eigen(m, iters=100):
    """对称 3x3 矩阵的特征分解，返回 (特征值, 特征向量按列)。"""
    a = [row[:] for row in m]
    v = [[1.0 if i == j else 0.0 for j in range(3)] for i in range(3)]
    for _ in range(iters):
        off = abs(a[0][1]) + abs(a[0][2]) + abs(a[1][2])
        if off < 1e-14:
            break
        for p, q in ((0, 1), (0, 2), (1, 2)):
            if abs(a[p][q]) < 1e-18:
                continue
            theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q])
            t = (1.0 if theta >= 0 else -1.0) / (abs(theta) + math.sqrt(theta * theta + 1))
            c = 1.0 / math.sqrt(t * t + 1)
            s = t * c
            for k in range(3):
                akp = a[k][p]
                akq = a[k][q]
                a[k][p] = c * akp - s * akq
                a[k][q] = s * akp + c * akq
            for k in range(3):
                apk = a[p][k]
                aqk = a[q][k]
                a[p][k] = c * apk - s * aqk
                a[q][k] = s * apk + c * aqk
            for k in range(3):
                vkp = v[k][p]
                vkq = v[k][q]
                v[k][p] = c * vkp - s * vkq
                v[k][q] = s * vkp + c * vkq
    vals = [a[i][i] for i in range(3)]
    vecs = [[v[i][j] for j in range(3)] for i in range(3)]
    return vals, vecs


def kabsch_rotation(a, b):
    """返回把 a 旋转到最贴近 b 的旋转矩阵 R。"""
    h = [[0.0] * 3 for _ in range(3)]
    for pa, pb in zip(a, b):
        va = (pa[1], pa[2], pa[3])
        vb = (pb[1], pb[2], pb[3])
        for i in range(3):
            for j in range(3):
                h[i][j] += va[i] * vb[j]
    # 简化：用 SVD 的替代方案——先求 H^T H 的特征分解构造旋转
    hth = [[sum(h[k][i] * h[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
    vals, vecs = jacobi_eigen(hth)
    # 对 H 做极分解 R = H (H^T H)^(-1/2)
    inv_sqrt = [[0.0] * 3 for _ in range(3)]
    for i in range(3):
        if vals[i] > 1e-12:
            for r in range(3):
                for c in range(3):
                    inv_sqrt[r][c] += vecs[r][i] * vecs[c][i] / math.sqrt(vals[i])
    r = [[sum(h[i][k] * inv_sqrt[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
    return r


def apply_rot(atoms, r):
    out = []
    for a in atoms:
        v = (a[1], a[2], a[3])
        out.append((a[0],
                    sum(r[0][k] * v[k] for k in range(3)),
                    sum(r[1][k] * v[k] for k in range(3)),
                    sum(r[2][k] * v[k] for k in range(3))))
    return out


def main():
    ap = argparse.ArgumentParser(description="比较两个 xyz 的 RMSD")
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--no-align", action="store_true", help="不做刚体对齐")
    args = ap.parse_args()

    a = read_xyz(args.a)
    b = read_xyz(args.b)

    if len(a) != len(b):
        print(f"原子数不同: {len(a)} vs {len(b)}", file=sys.stderr)
        return 1
    for pa, pb in zip(a, b):
        if pa[0] != pb[0]:
            print(f"原子顺序不一致: {pa[0]} vs {pb[0]}", file=sys.stderr)
            return 1

    ca, cb = centroid(a), centroid(b)
    ac = centered(a, ca)
    bc = centered(b, cb)

    raw = rmsd_direct(ac, bc)

    if args.no_align:
        print(f"未对齐 RMSD: {raw:.4f} A")
        return 0

    r = kabsch_rotation(ac, bc)
    aligned = apply_rot(ac, r)
    aligned_rmsd = rmsd_direct(aligned, bc)

    print(f"对齐前 RMSD: {raw:.4f} A")
    print(f"对齐后 RMSD: {aligned_rmsd:.4f} A")
    print(f"原子数: {len(a)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""gjf_gen.py - 从 xyz 文件生成 Gaussian 16 输入文件。

用法：
    python gjf_gen.py mol.xyz --job opt-freq
    python gjf_gen.py mol.xyz --job ts --method m062x --basis 6-311++g(d,p)
    python gjf_gen.py mol.xyz --job sp --nproc 12 --mem 8GB

支持的任务模板：
    opt        # opt
    opt-freq   # opt freq
    ts         # opt=(calcfc,ts,noeigen) freq
    sp         # 单点

生成的输入包含 %nprocshared / %mem / %chk / 坐标 / 末尾空行，
connectivity 需要另加（Gaussian 会用 guess 自己判断）。
"""
import argparse
import os
import sys

TEMPLATES = {
    "opt": "# opt {method}/{basis}",
    "opt-freq": "# opt freq {method}/{basis}",
    "ts": "# opt=(calcfc,ts,noeigen) freq {method}/{basis}",
    "sp": "# {method}/{basis}",
}


def read_xyz(path):
    with open(path, "r") as f:
        lines = [ln.rstrip("\n") for ln in f if ln.strip()]
    n = int(lines[0].split()[0])
    coords = [ln for ln in lines[2:2 + n]]
    if len(coords) != n:
        raise ValueError(f"{path}: 声明 {n} 个原子，实际 {len(coords)}")
    return coords


def main():
    ap = argparse.ArgumentParser(description="生成 Gaussian 16 输入")
    ap.add_argument("xyz", help="xyz 文件")
    ap.add_argument("--job", default="opt-freq", choices=sorted(TEMPLATES))
    ap.add_argument("--method", default="b3lyp")
    ap.add_argument("--basis", default="6-311+g(d,p)")
    ap.add_argument("--charge", type=int, default=0)
    ap.add_argument("--mult", type=int, default=1)
    ap.add_argument("--nproc", type=int, default=12)
    ap.add_argument("--mem", default="8GB")
    ap.add_argument("--out", help="输出 gjf 路径（默认同名 .gjf）")
    ap.add_argument("--title", default="Title Card Required")
    args = ap.parse_args()

    coords = read_xyz(args.xyz)
    base = os.path.splitext(os.path.abspath(args.xyz))[0]
    out_path = args.out or (base + ".gjf")
    chk = base + ".chk"

    route = TEMPLATES[args.job].format(method=args.method, basis=args.basis)

    with open(out_path, "w", newline="\n") as f:
        f.write(f"%nprocshared={args.nproc}\n")
        f.write(f"%mem={args.mem}\n")
        f.write(f"%chk={chk}\n")
        f.write(route + "\n")
        f.write("\n")
        f.write(args.title + "\n")
        f.write("\n")
        f.write(f"{args.charge} {args.mult}\n")
        for ln in coords:
            f.write(ln + "\n")
        f.write("\n")

    print(f"已生成: {out_path}")
    print(f"  任务: {args.job}  方法: {args.method}/{args.basis}")
    print(f"  原子数: {len(coords)}  电荷/多重度: {args.charge}/{args.mult}")
    print(f"  chk: {chk}")
    print("提示：如果原 xyz 带 connectivity，请手动补到坐标后面。")
    return 0


if __name__ == "__main__":
    sys.exit(main())

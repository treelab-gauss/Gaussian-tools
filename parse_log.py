#!/usr/bin/env python3
"""parse_log.py - 从 Gaussian 16 的 .log/.out 文件里抽取关键结果。

用法：
    python parse_log.py job.log
    python parse_log.py *.log
    python parse_log.py job.log --json

抽取内容：
    - 任务是否正常结束（Normal / Error termination）
    - 最终 SCF 能量（Hartree）
    - 原子数、电荷、自旋多重度
    - 虚频个数和虚频数值（cm-1）
    - 优化步数、最终力收敛指标
    - 墙钟时间

只依赖标准库。
"""
import argparse
import glob
import json
import os
import re
import sys

RE_NATOMS = re.compile(r"NAtoms=\s*(\d+)")
RE_CHARGE = re.compile(r"Charge\s*=\s*(-?\d+)\s+Multiplicity\s*=\s*(\d+)")
RE_SCF = re.compile(r"SCF Done:\s+E\(([^)]+)\)\s*=\s*(-?\d+\.\d+)")
RE_STEP = re.compile(r"Step number\s+(\d+)\s+out of a maximum of\s+(\d+)")
RE_MAXF = re.compile(r"Maximum Force\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s+(YES|NO)")
RE_RMSF = re.compile(r"RMS\s+Force\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s+(YES|NO)")
RE_TERM = re.compile(r"(Normal|Error) termination")
RE_TIME = re.compile(r"Job cpu time:\s*(\d+)\s*days\s*(\d+)\s*hours\s*(\d+)\s*minutes\s*(\d+\.\d+)\s*seconds")
RE_HARMONIC = re.compile(r"Harmonic frequencies \(cm\^-1\):")
RE_FREQ = re.compile(r"^\s+Frequencies --\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)")


def parse(path):
    info = {
        "file": os.path.basename(path),
        "termination": None,
        "natoms": None,
        "charge": None,
        "multiplicity": None,
        "scf_energy": None,
        "scf_method": None,
        "n_scf": 0,
        "opt_step": None,
        "opt_max_steps": None,
        "max_force": None,
        "max_force_ok": None,
        "rms_force": None,
        "rms_force_ok": None,
        "n_imag": None,
        "imag_freqs": [],
        "cpu_time": None,
    }
    try:
        with open(path, "r", errors="ignore") as f:
            lines = f.readlines()
    except OSError as e:
        info["error"] = str(e)
        return info

    in_freq = False
    for line in lines:
        m = RE_NATOMS.search(line)
        if m:
            info["natoms"] = int(m.group(1))

        m = RE_CHARGE.search(line)
        if m:
            info["charge"] = int(m.group(1))
            info["multiplicity"] = int(m.group(2))

        m = RE_SCF.search(line)
        if m:
            info["scf_method"] = m.group(1)
            info["scf_energy"] = float(m.group(2))
            info["n_scf"] += 1

        m = RE_STEP.search(line)
        if m:
            info["opt_step"] = int(m.group(1))
            info["opt_max_steps"] = int(m.group(2))

        m = RE_MAXF.search(line)
        if m:
            info["max_force"] = float(m.group(1))
            info["max_force_ok"] = m.group(3) == "YES"

        m = RE_RMSF.search(line)
        if m:
            info["rms_force"] = float(m.group(1))
            info["rms_force_ok"] = m.group(3) == "YES"

        m = RE_TERM.search(line)
        if m:
            info["termination"] = m.group(1)

        m = RE_TIME.search(line)
        if m:
            d, h, mi, s = (float(x) for x in m.groups())
            info["cpu_time"] = round(d * 86400 + h * 3600 + mi * 60 + s, 1)

        if RE_HARMONIC.search(line):
            in_freq = True
            continue
        if in_freq:
            m = RE_FREQ.search(line)
            if m:
                for g in m.groups():
                    v = float(g)
                    if v < 0:
                        info["imag_freqs"].append(v)

    if info["imag_freqs"]:
        info["n_imag"] = len(info["imag_freqs"])
    elif info["termination"] == "Normal":
        # 正常结束且没有负频率 -> 0 个虚频（极小点）
        info["n_imag"] = 0 if info["opt_step"] is not None else info["n_imag"]
    return info


def fmt(info):
    if info.get("error"):
        return f"{info['file']}: 读取失败 {info['error']}"

    lines = [f"文件: {info['file']}"]
    lines.append(f"  结束状态: {info['termination'] or '未找到'}")
    if info["natoms"]:
        lines.append(f"  原子数: {info['natoms']}  电荷 {info['charge']}  多重度 {info['multiplicity']}")
    if info["opt_step"]:
        lines.append(f"  优化步数: {info['opt_step']}/{info['opt_max_steps']}")
    if info["scf_energy"] is not None:
        lines.append(f"  最终能量: {info['scf_energy']:.8f} Ha  ({info['scf_method']}, {info['n_scf']} 次 SCF)")
    if info["max_force"] is not None:
        lines.append(f"  最终力: max {info['max_force']:.6f} ({'收敛' if info['max_force_ok'] else '未收敛'})"
                     f" / RMS {info['rms_force']:.6f} ({'收敛' if info['rms_force_ok'] else '未收敛'})")
    if info["n_imag"] is not None:
        if info["n_imag"] == 0:
            lines.append("  虚频: 0 个（极小点）")
        else:
            vals = ", ".join(f"{v:.1f}" for v in info["imag_freqs"])
            lines.append(f"  虚频: {info['n_imag']} 个（{vals}）")
    if info["cpu_time"]:
        h = info["cpu_time"] / 3600
        lines.append(f"  CPU 时间: {h:.2f} 小时")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description="解析 Gaussian 16 日志")
    ap.add_argument("files", nargs="+", help="日志文件，支持通配符")
    ap.add_argument("--json", action="store_true", help="输出 JSON")
    args = ap.parse_args()

    paths = []
    for pat in args.files:
        paths.extend(glob.glob(pat))
    if not paths:
        print("没有匹配的文件", file=sys.stderr)
        return 1

    results = [parse(p) for p in sorted(paths)]
    if args.json:
        print(json.dumps(results, ensure_ascii=False, indent=2))
    else:
        for i, r in enumerate(results):
            if i:
                print()
            print(fmt(r))
    return 0


if __name__ == "__main__":
    sys.exit(main())

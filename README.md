# Gaussian-tools

一组处理 Gaussian 16 计算的小工具，纯 Python 标准库，不需要 numpy。

写这套东西的动因：日常跑量化计算时，生成输入、看收敛、比结构这些事都是重复劳动，
用脚本串起来之后，一个任务的准备和检查从十几分钟压到几十秒。

## 脚本

### gjf_gen.py —— 从 xyz 生成 Gaussian 输入

```bash
python gjf_gen.py mol.xyz --job opt-freq
python gjf_gen.py mol.xyz --job ts --method m062x --basis 6-311++g(d,p)
python gjf_gen.py mol.xyz --job sp --nproc 16 --mem 16GB
```

支持四种任务模板：

| --job | 路由 |
|-------|------|
| opt | `# opt {method}/{basis}` |
| opt-freq | `# opt freq {method}/{basis}` |
| ts | `# opt=(calcfc,ts,noeigen) freq {method}/{basis}` |
| sp | `# {method}/{basis}` |

自动带入 `%nprocshared`、`%mem`、`%chk`，电荷和多重度可指定。

### parse_log.py —— 从 .log 抽结果

```bash
python parse_log.py job.log
python parse_log.py *.log
python parse_log.py job.log --json
```

一次抽取：结束状态（Normal / Error）、原子数、电荷、多重度、最终 SCF 能量、
优化步数、最终力是否收敛、虚频个数和数值、CPU 时间。

批量跑一堆 log 时特别省事，一眼看出哪个正常结束、哪个有虚频。

### xyz_rmsd.py —— 比较两个结构

```bash
python xyz_rmsd.py a.xyz b.xyz
python xyz_rmsd.py a.xyz b.xyz --no-align
```

先用 Kabsch 算法做刚体对齐（旋转 + 平移），再算 RMSD。用来判断
优化前后差多少、过渡态和产物差多少。

## 为什么不用现成的工具

GaussView 能看，但批量处理要一个个点。MultiWfn 很强，但那是分析波函数的，
不是做流水线的。这几个脚本是胶水层，配合现有流程用。

## 已知限制

1. `parse_log.py` 只处理 Gaussian 16 的英文输出格式。
2. `gjf_gen.py` 不生成 connectivity，复杂体系请手动补。
3. `xyz_rmsd.py` 要求两个文件原子顺序一致，原子符号会检查。

## 环境

Python 3.8+，无第三方依赖。

```bash
python -m py_compile *.py   # 语法自检
```

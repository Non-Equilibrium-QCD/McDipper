#!/usr/bin/env python3
# Slurm settings and module loads below are for the PC2 cluster; adapt them to yours.
import os
import subprocess
import sys

ETAMIN, DETA, NETA = -12.0, 0.2, 121
HALF = (NETA - 1) // 2          # index of y = 0

TEMPLATE = open("configs/table_template.yaml").read()

# Appended to the config only with --reuse.
REUSE_CONFIG = "    Reuse: True\n    ReuseTable: resume_gluon.dat\n"

# Inserted into the job script only with --reuse: keep the rows the killed
# job wrote and start from a clean Set0.
REUSE_SETUP = """
if [ -f tabs/IPSat/Set0/gluon_energy_density.dat ]; then
    mv tabs/IPSat/Set0/gluon_energy_density.dat resume_gluon.dat
fi
rm -rf tabs/IPSat
"""

JOB = """#!/bin/bash
#SBATCH --job-name=mcdip_table_eta{index:03d}
#SBATCH --output={repo}/log/%x_%j.out
#SBATCH --time=24:00:00
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={threads}
#SBATCH --mem=8G

cd {etadir}
echo "Start time: $(date)   eta = {eta}"
{reuse_setup}
module purge
module load compiler/GCC/10.3.0
module load numlib/GSL/2.7-GCC-10.3.0
module load lang/Python/3.9.5-GCCcore-10.3.0
module load devel/CMake/3.20.1-GCCcore-10.3.0
module load devel/Boost/1.76.0-GCC-10.3.0
module load math/Eigen/3.3.9-GCCcore-10.3.0
module load lib/libglvnd/1.3.3-GCCcore-10.3.0
module load data/HDF5/1.10.7-gompi-2021a

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export LD_LIBRARY_PATH=$HOME/LHAPDF/lib:$LD_LIBRARY_PATH
export PATH=$HOME/LHAPDF/bin:$PATH

./mcDipper -i configs/table.yaml

echo "End time: $(date)"
"""


def parse_args(argv):
    reuse = False
    threads = 64
    numbers = []
    for arg in argv:
        if arg == "--reuse":
            reuse = True
        elif arg.startswith("--threads="):
            threads = int(arg.split("=")[1])
        else:
            numbers.append(int(arg))
    # Only y >= 0 is computed; merge_table.py mirrors the rest.
    return reuse, threads, numbers or list(range(HALF, NETA))


def prepare_dir(repo, i):
    # One directory per rapidity: check_sets writes into the first
    # tabs/IPSat/SetN that does not exist, so jobs sharing a directory
    # would all pick the same set and overwrite each other.
    d = repo + "/table_build/eta_%03d" % i
    os.makedirs(d + "/configs", exist_ok=True)
    os.makedirs(d + "/tabs", exist_ok=True)

    # Copied, not linked: rebuilding while jobs are queued would otherwise
    # change the binary under them. Earlier runs left a symlink here, which cp
    # refuses to overwrite from its own target.
    binary = d + "/mcDipper"
    if os.path.islink(binary):
        os.remove(binary)
    subprocess.call(["cp", "-p", repo + "/mcDipper", binary])

    for link, target in [(d + "/ipsat", repo + "/ipsat"),
                         (d + "/tabs/IPSat-Dipole", repo + "/tabs/IPSat-Dipole")]:
        if not os.path.islink(link):
            os.symlink(target, link)
    return d


def write_config(d, i, eta, reuse):
    text = TEMPLATE.format(eta="%g" % eta, index=i)
    if reuse:
        text += REUSE_CONFIG
    else:
        # Drop a seed table left by an earlier run, so it cannot be picked
        # up later by mistake.
        stale = d + "/resume_gluon.dat"
        if os.path.exists(stale):
            os.remove(stale)
    open(d + "/configs/table.yaml", "w").write(text)


def write_job(d, repo, i, eta, reuse, threads):
    path = d + "/job.sh"
    open(path, "w").write(JOB.format(
        index=i, eta="%g" % eta, repo=repo, etadir=d, threads=threads,
        reuse_setup=REUSE_SETUP if reuse else ""))
    return path


repo = os.getcwd()
reuse, threads, indices = parse_args(sys.argv[1:])

for i in indices:
    eta = round(ETAMIN + i * DETA, 4)
    d = prepare_dir(repo, i)
    write_config(d, i, eta, reuse)
    subprocess.call(["sbatch", write_job(d, repo, i, eta, reuse, threads)])

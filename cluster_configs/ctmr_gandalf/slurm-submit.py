#!/usr/bin/env python3
import sys
import re
import argparse
import subprocess
from pathlib import Path

from snakemake.utils import read_job_properties

parser = argparse.ArgumentParser(add_help=False)
parser.add_argument(
    "--help", help="Display help message.", action="store_true")
parser.add_argument(
    "positional", action="append",
    nargs="?", metavar="POS",
    help="additional arguments not in slurm parser group to pass to sbatch")

# A subset of SLURM-specific arguments
slurm_parser = parser.add_argument_group("slurm-specific arguments")
slurm_parser.add_argument(
    "-a", "--array", help="job array index values")
slurm_parser.add_argument(
    "-A", "--account", help="charge job to specified account")
slurm_parser.add_argument(
    "--begin", help="defer job until HH:MM MM/DD/YY")
slurm_parser.add_argument(
    "-c", "--cpus-per-task", help="number of cpus required per task")
slurm_parser.add_argument(
    "-d", "--dependency",
    help="defer job until condition on jobid is satisfied")
slurm_parser.add_argument(
    "-D", "--workdir", help="set working directory for batch script")
slurm_parser.add_argument(
    "-e", "--error", help="file for batch script's standard error")
slurm_parser.add_argument(
    "-J", "--job-name", help="name of job")
slurm_parser.add_argument(
    "--mail-type", help="notify on state change: BEGIN, END, FAIL or ALL")
slurm_parser.add_argument(
    "--mail-user", help="who to send email notification for job state changes")
slurm_parser.add_argument(
    "-n", "--ntasks", help="number of tasks to run")
slurm_parser.add_argument(
    "-N", "--nodes", help="number of nodes on which to run (N = min[-max])")
slurm_parser.add_argument(
    "-o", "--output", help="file for batch script's standard output")
slurm_parser.add_argument(
    "-p", "--partition", help="partition requested")
slurm_parser.add_argument(
    "-q", "--qos", help="quality of service")
slurm_parser.add_argument(
    "-Q", "--quiet", help="quiet mode (suppress informational messages)")
slurm_parser.add_argument(
    "-t", "--time", help="time limit")
slurm_parser.add_argument(
    "--wrap", help="wrap command string in a sh script and submit")
slurm_parser.add_argument(
    "-C", "--constraint", help="specify a list of constraints")
slurm_parser.add_argument(
    "--mem", help="minimum amount of real memory")

args = parser.parse_args()

if args.help:
    parser.print_help()
    sys.exit(0)

jobscript = sys.argv[-1]
job_properties = read_job_properties(jobscript)

extras = ""
if args.positional:
    for m in args.positional:
        if m is not None:
            extras = extras + " " + m

arg_dict = dict(args.__dict__)


# Process resources
if "resources" in job_properties:
    resources = job_properties["resources"]
    if arg_dict["time"] is None:
        if "runtime" in resources:
            arg_dict["time"] = resources["runtime"]
        elif "walltime" in resources:
            arg_dict["time"] = resources["walltime"]
    if "mem" in resources and arg_dict["mem"] is None:
        arg_dict["mem"] = resources["mem"]

# Threads
if "threads" in job_properties:
    arg_dict["ntasks"] = job_properties["threads"]

opt_keys = {"array", "account", "begin", "cpus_per_task",
            "depedency", "workdir", "error", "job_name", "mail_type",
            "mail_user", "ntasks", "nodes", "output", "partition",
            "quiet", "time", "wrap", "constraint", "mem"}

# Set default partition
if arg_dict["partition"] is None:
    # partitions and SLURM - If not specified, the default behavior is to
    # allow the slurm controller to select the default partition as
    # designated by the system administrator.
    opt_keys.remove("partition")

# NOT APPLICABLE FOR CTMR UCP
#### Set default account
###if arg_dict["account"] is None:
###    raise Exception("Cannot submit Slurm jobs without account!")

# Ensure output folder for Slurm log files exist
# This is a bit hacky; it will try to create the folder
# for every Slurm submission...
if "output" in arg_dict:
    stdout_folder = Path(arg_dict["output"]).parent
    stdout_folder.mkdir(exist_ok=True, parents=True)
if "error" in arg_dict:
    stdout_folder = Path(arg_dict["error"]).parent
    stdout_folder.mkdir(exist_ok=True, parents=True)

opts = ""
for k, v in arg_dict.items():
    if k not in opt_keys:
        continue
    if v is not None:
        opts += " --{} \"{}\" ".format(k.replace("_", "-"), v)

if arg_dict["wrap"] is not None:
    cmd = "sbatch {opts}".format(opts=opts)
else:
    cmd = "sbatch {opts} {extras}".format(opts=opts, extras=extras)

try:
    res = subprocess.run(cmd, check=True, shell=True, stdout=subprocess.PIPE)
except subprocess.CalledProcessError as e:
    raise e

# Get jobid
res = res.stdout.decode()
try:
    m = re.search("Submitted batch job (\d+)", res)
    jobid = m.group(1)
    print(jobid)
except Exception as e:
    print(e)
    raise

from jinja2 import Environment, select_autoescape, FileSystemLoader
import pandas as pd
from pathlib import Path
import subprocess
import datetime
import os
import uuid
import argparse
import sys


def render(template_path, precice_config_params):
    base_path = Path(__file__).parent.absolute()

    env = Environment(
        loader=FileSystemLoader(base_path),
        autoescape=select_autoescape(['xml'])
    )

    precice_config_template = env.get_template(template_path)

    precice_config_name = base_path / "precice-config.xml"

    with open(precice_config_name, "w") as file:
        file.write(precice_config_template.render(precice_config_params))


def do_run(template_path, precice_config_params, participants):
    render(template_path, precice_config_params)
    print(f"{datetime.datetime.now()}: Start run with parameters {precice_config_params}")
    print("Running...")

    for pname, participant in participants.items():
        participant['logfile'] = f"stdout-{pname}.log"

    for participant in participants.values():
        with open(participant['root'] / participant['logfile'], "w") as outfile:
            cmd = participant["exec"] + participant["params"] + [f"{keyword}={value}" for keyword, value in participant['kwargs'].items()]
            p = subprocess.Popen(cmd,
                                 cwd=participant['root'],
                                 stdout=outfile)
            participant["proc"] = p

    for participant in participants.values():
        participant["proc"].wait()

    for participant in participants.values():
        if participant["proc"].returncode != 0:
            raise Exception(f'Experiment failed. See logs {[p["logfile"] for p in participants]}')

    print("Done.")
    print("Postprocessing...")
    time_window_size = precice_config_params['time_window_size']
    summary = {"time window size": time_window_size}

    for pname, participant in participants.items():
        summary[f"time step size {pname}"] = time_window_size

        t_end = precice_config_params['max_time']
        qoi = "Displacement0"  # quantity of interest


        if pname == "Solid":
            df_ref = pd.read_csv(f"watchpoint_{participant['case']}_ref", comment="#", delim_whitespace=True)
            try:
                qoi_ref_at_end = df_ref[df_ref["Time"]==t_end][qoi].to_list()[-1]
            except IndexError:
                qoi_ref_at_end = -1

            df = pd.read_csv(participant['root'] / f"precice-{pname}-watchpoint-Flap-Tip.log", comment="#", delim_whitespace=True)
            qoi_at_end = df[df["Time"]==t_end][qoi].to_list()[-1]
            summary[f"{qoi} {pname}"] = qoi_at_end
            summary[f"error {pname}"] = abs(qoi_at_end - qoi_ref_at_end)
        elif pname == "Fluid":
            pass  # watchpoint is empty for fluid-fake

    print("Done.")

    return summary


if __name__ == "__main__":
    n_supported_participants = 2

    parser = argparse.ArgumentParser(description="Solving perpendicular flap")
    parser.add_argument(
        "template_path",
        help="template for the preCICE configuration file",
        type=str)
    parser.add_argument(
        "-T",
        "--max-time",
        help="Max simulation time",
        type=float,
        default=5.0)
    parser.add_argument(
        "-dt",
        "--base-time-window-size",
        help="Base time window / time step size",
        type=float,
        default=0.001)
    parser.add_argument(
        "-w",
        "--time-window-refinements",
        help="Number of refinements by factor 2 for the time window size",
        type=int,
        default=1)

    args = parser.parse_args()

    df = pd.DataFrame()

    precice_config_params = {
        'time_window_size': None,  # will be defined later
        'max_time': args.max_time,
    }

    root_folder = Path(__file__).parent.absolute()

    participants = {
        "Fluid": {
            "case": "fluid-fake",
            "exec": ["./run.sh"],  # how to execute the participant, e.g. python3 script.py
            "params": [],  # list of positional arguments that will be used. Results in python3 script.py param1 ...
            "kwargs": {  # dict with keyword arguments that will be used. Results in python3 script.py param1 ... k1=v1 k2=v2 ...
            },
        },
        "Solid": {
            "case": "solid-fenics",
            "exec": ["./run.sh"],  # how to execute the participant, e.g. python3 script.py
            "params": [],  # list of positional arguments that will be used. Results in python3 script.py param1 ...
            "kwargs": {  # dict with keyword arguments that will be used. Results in python3 script.py param1 ... k1=v1 k2=v2 ...
            },
        },
    }

    run_id = uuid.uuid4()
    summary_file = Path("convergence-studies") / f"{run_id}.csv"
    watchpoint_folder = Path("convergence-studies") / str(run_id)
    watchpoint_folder.mkdir(parents=False, exist_ok=False)

    for dt in [args.base_time_window_size * 0.5**i for i in range(args.time_window_refinements)]:
        precice_config_params['time_window_size'] = dt

        for pname in participants.keys():
            participants[pname]["root"] = root_folder / participants[pname]["case"]


        summary = do_run(args.template_path, precice_config_params, participants)
        df = pd.concat([df, pd.DataFrame(summary, index=[0])], ignore_index=True)

        # store the watchpoint file
        (participants["Solid"]["root"] / "precice-Solid-watchpoint-Flap-Tip.log").rename(watchpoint_folder / f"watchpoint_{dt}")

        print(f"Write preliminary output to {summary_file}")
        df.to_csv(summary_file)

        term_size = os.get_terminal_size()
        print('-' * term_size.columns)
        print(df)
        print('-' * term_size.columns)

    df = df.set_index(['time window size'] + [f'time step size {p}' for p in participants.keys()])
    print(f"Write final output to {summary_file}")

    import git
    import precice

    repo = git.Repo(__file__, search_parent_directories=True)
    chash = str(repo.head.commit)[:7]
    if repo.is_dirty():
        chash += "-dirty"

    metadata = {
        "git repository": repo.remotes.origin.url,
        "git commit": chash,
        "precice.get_version_information()": precice.get_version_information(),
        "precice.__version__": precice.__version__,
        "run cmd": "python3 " + " ".join(sys.argv),
        "args": args,
        "precice_config_params": precice_config_params,
        "participants": participants,
    }

    summary_file.unlink()

    with open(summary_file, 'a') as f:
        for key, value in metadata.items():
            f.write(f"# {key}:{value}\n")
        df.to_csv(f)

    print('-' * term_size.columns)
    for key, value in metadata.items():
        print(f"{key}:{value}")
    print()
    print(df)
    print('-' * term_size.columns)

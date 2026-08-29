import hashlib
import subprocess
import threading
from .sources import resolve_tutorial_root, PRECICE_EXTERNAL_CACHE_DIR
from typing import List, Dict, Optional, Tuple
from jinja2 import Environment, FileSystemLoader
from dataclasses import dataclass, field
import shutil
from pathlib import Path
from paths import PRECICE_REL_OUTPUT_DIR, PRECICE_TOOLS_DIR, PRECICE_REL_REFERENCE_DIR, PRECICE_TESTS_DIR, PRECICE_TUTORIAL_DIR

from metadata_parser.metdata import Tutorial, CaseCombination, Case, ReferenceResult
from .SystemtestArguments import SystemtestArguments

from datetime import datetime, timedelta
import tarfile
import time

import unicodedata
import re
import logging
import os


GLOBAL_TIMEOUT = int(os.environ.get("PRECICE_SYSTEMTESTS_TIMEOUT", 180))
DEFAULT_BUILD_TIMEOUT = int(
    os.environ.get("PRECICE_SYSTEMTESTS_BUILD_TIMEOUT", 480))
DEFAULT_FIELDCOMPARE_RTOL = 3e-7
SHORT_TIMEOUT = 10

DIFF_RESULTS_DIR = "diff-results"
ITERATIONS_LOGS_DIR = "iterations-logs"
DIFF_VISUALIZER_LOG = "system-tests-compare-diff.log"
DIFF_VISUALIZER_TIMEOUT = int(
    os.environ.get("PRECICE_SYSTEMTESTS_DIFF_VISUALIZER_TIMEOUT", 900)
)

STAGE_LOG_FILES = {
    "build": "system-tests-build.log",
    "run": "system-tests-run.log",
    "compare": "system-tests-compare.log",
}

FAILURE_LOG_TAIL_LINES = 100


class _SystemtestLogSink:
    """Writes subprocess output incrementally to per-stage log files."""

    def __init__(self, system_test_dir: Path):
        self._system_test_dir = system_test_dir
        self._lock = threading.Lock()

    def begin_stage(self, stage: str) -> None:
        stage_path = self._system_test_dir / STAGE_LOG_FILES[stage]
        stage_path.write_text(f"=== {stage} ===\n", encoding="utf-8")

    def append_stdout(self, line: str, stage: str) -> None:
        with self._lock:
            stage_path = self._system_test_dir / STAGE_LOG_FILES[stage]
            with stage_path.open("a", encoding="utf-8") as log_file:
                log_file.write(line + "\n")

    def append_stderr(self, line: str, stage: str) -> None:
        with self._lock:
            stage_path = self._system_test_dir / STAGE_LOG_FILES[stage]
            with stage_path.open("a", encoding="utf-8") as log_file:
                log_file.write(f"[stderr] {line}\n")


def slugify(value, allow_unicode=False):
    """
    Taken from https://github.com/django/django/blob/master/django/utils/text.py
    Convert to ASCII if 'allow_unicode' is False. Convert spaces or repeated
    dashes to single dashes. Remove characters that aren't alphanumerics,
    underscores, or hyphens. Convert to lowercase. Also strip leading and
    trailing whitespace, dashes, and underscores.
    """
    value = str(value)
    if allow_unicode:
        value = unicodedata.normalize('NFKC', value)
    else:
        value = unicodedata.normalize('NFKD', value).encode(
            'ascii', 'ignore').decode('ascii')
    value = re.sub(r'[^\w\s-]', '', value.lower())
    return re.sub(r'[-\s]+', '-', value).strip('-_')


class Systemtest:
    pass


@dataclass
class DockerComposeResult:
    exit_code: int
    stdout_data: List[str]
    stderr_data: List[str]
    systemtest: Systemtest
    runtime: float  # in seconds


@dataclass
class FieldCompareResult:
    exit_code: int
    stdout_data: List[str]
    stderr_data: List[str]
    systemtest: Systemtest
    runtime: float  # in seconds


@dataclass
class SystemtestResult:
    success: bool
    stdout_data: List[str]
    stderr_data: List[str]
    systemtest: Systemtest
    build_time: float  # in seconds
    solver_time: float  # in seconds
    fieldcompare_time: float  # in seconds


def _success_status_symbol(success: bool) -> str:
    return "✅" if success else "❌"


def _read_log_tail(log_path: Path, max_lines: int = FAILURE_LOG_TAIL_LINES) -> str:
    lines = log_path.read_text(encoding="utf-8", errors="replace").splitlines()
    if not lines:
        return "(log file is empty)"
    return "\n".join(lines[-max_lines:])


def _append_failure_log_tails_to_summary(results: List[SystemtestResult]) -> None:
    summary_path = os.environ.get("GITHUB_STEP_SUMMARY")
    if not summary_path:
        return

    failed_results = [result for result in results if not result.success]
    if not failed_results:
        return

    with open(summary_path, "a", encoding="utf-8") as summary_file:
        print("\n## Failed test logs\n", file=summary_file)
        for result in failed_results:
            print(
                f"### {_success_status_symbol(False)} {result.systemtest}\n",
                file=summary_file,
            )
            run_dir = result.systemtest.get_system_test_dir()
            for log_name in STAGE_LOG_FILES.values():
                log_path = run_dir / log_name
                if not log_path.is_file():
                    continue
                tail = _read_log_tail(log_path)
                print("<details>", file=summary_file)
                print(f"<summary>{log_name} tail</summary>", file=summary_file)
                print("", file=summary_file)
                print("```text", file=summary_file)
                print(tail, file=summary_file)
                print("```", file=summary_file)
                print("</details>", file=summary_file)
                print("", file=summary_file)


def display_systemtestresults_as_table(results: List[SystemtestResult]):
    """
    Prints the result in a nice tabluated way to get an easy overview
    """
    print()

    def _get_length_of_name(results: List[SystemtestResult]) -> int:
        return max(len(str(result.systemtest)) for result in results)

    max_name_length = _get_length_of_name(results)

    header = f"| {'systemtest':<{max_name_length + 2}} "\
        f"| {'status':^7} "\
        f"| {'build':^11} "\
        f"| {'run':^11} "\
        f"| {'compare':^11} |"
    separator_plaintext = "+-" + "-" * (max_name_length + 2) + \
        "-+---------+-------------+-------------+-------------+"
    separator_markdown = "| --- | --- | --- | --- | --- |"

    print(separator_plaintext)
    print(header)
    print(separator_plaintext)

    if "GITHUB_STEP_SUMMARY" in os.environ:
        with open(os.environ["GITHUB_STEP_SUMMARY"], "a") as f:
            print(header, file=f)
            print(separator_markdown, file=f)

    for result in results:
        build_time = int(timedelta(seconds=result.build_time).total_seconds())
        build_time_m, build_time_s = divmod(build_time, 60)

        solver_time = int(timedelta(seconds=result.solver_time).total_seconds())
        solver_time_m, solver_time_s = divmod(solver_time, 60)

        fieldcompare_time = int(timedelta(seconds=result.fieldcompare_time).total_seconds())
        fieldcompare_time_m, fieldcompare_time_s = divmod(fieldcompare_time, 60)

        row = f"| {str(result.systemtest):<{max_name_length + 2}} "\
            f"| {_success_status_symbol(result.success):^5} "\
            f"|     {build_time_m:>2}m {build_time_s:02d}s "\
            f"|     {solver_time_m:>2}m {solver_time_s:02d}s "\
            f"|     {fieldcompare_time_m:>2}m {fieldcompare_time_s:02d}s |"
        print(row)
        print(separator_plaintext)
        if "GITHUB_STEP_SUMMARY" in os.environ:
            with open(os.environ["GITHUB_STEP_SUMMARY"], "a") as f:
                print(row, file=f)

    _append_failure_log_tails_to_summary(results)

    if "GITHUB_STEP_SUMMARY" in os.environ:
        with open(os.environ["GITHUB_STEP_SUMMARY"], "a") as f:
            print("\n\n", file=f)
            print(
                "In case a test fails, download the archive from the bottom of this page and inspect the per-stage logs (`system-tests-build.log`, `system-tests-run.log`, `system-tests-compare.log`, and `system-tests-compare-diff.log`). The stage runtimes might already give useful hints.",
                file=f)
            print(
                "See the [documentation](https://precice.org/dev-docs-system-tests.html#understanding-what-went-wrong).",
                file=f)


@dataclass
class Systemtest:
    """
    Represents a system test by specifing the cases and the corresponding Tutorial
    """

    tutorial: Tutorial
    arguments: SystemtestArguments
    case_combination: CaseCombination
    reference_result: ReferenceResult
    max_time: float | None = None
    max_time_windows: int | None = None
    timeout: int = GLOBAL_TIMEOUT
    tolerance: float = DEFAULT_FIELDCOMPARE_RTOL
    skip_compare: bool = False
    run_before: str | None = None
    run_after: str | None = None
    params_to_use: Dict[str, str] = field(init=False)
    build_arguments_by_component: Dict[str, Dict[str, str]] = field(init=False)
    env: Dict[str, str] = field(init=False)

    def __eq__(self, other) -> bool:
        if isinstance(other, Systemtest):
            return (
                self.tutorial == other.tutorial) and (
                self.arguments == other.arguments) and (
                self.case_combination == other.case_combination)
        return False

    def __hash__(self) -> int:
        return hash(f"{self.tutorial, self.arguments, self.case_combination}")

    def __post_init__(self):
        self.__init_args_to_use()
        self.env = {}
        self.build_timeout = self._resolve_build_timeout()

    def _resolve_build_timeout(self) -> int:
        """
        Wall-clock limit for the single ``docker compose build`` subprocess.

        Uses the maximum build_timeout of the distinct components in this test,
        so the step can run long enough for the slowest adapter. Components
        without build_timeout use DEFAULT_BUILD_TIMEOUT.
        """
        timeouts = []
        seen_components = set()
        for case in self.case_combination.cases:
            if case.component.name in seen_components:
                continue
            seen_components.add(case.component.name)
            if case.component.build_timeout is not None:
                timeouts.append(case.component.build_timeout)
            else:
                timeouts.append(DEFAULT_BUILD_TIMEOUT)
        return max(timeouts) if timeouts else DEFAULT_BUILD_TIMEOUT

    def __init_args_to_use(self):
        """
        Resolve build arguments per component from components.yaml defaults and
        global CLI overrides, then merge them for shared stages (prepare).

        Previously, this function was also checking if all required parameters were provided, and was raising exceptions for parameters not provided and not having a default value. This check made adding optional parameters with empty defaults (e.g., the TUTORIALS_PR) complicated, and it was removed.
        """
        provided_arguments = self.arguments.arguments
        explicit_cli_keys = frozenset(provided_arguments.keys())
        self.build_arguments_by_component = {}
        resolved_ref_cache: Dict[Tuple[str, str], str] = {}

        logging.debug(
            f"Substituting default build arguments and resolving git refs for {self}.")

        for case in self.case_combination.cases:
            component = case.component
            if component.name in self.build_arguments_by_component:
                continue

            logging.debug(f"Resolving build arguments for component {component.name}.")
            component_args: Dict[str, str] = {}
            for param in component.parameters:
                if param.key in explicit_cli_keys:
                    value = provided_arguments[param.key]
                elif (
                    param.key in provided_arguments
                    and len(provided_arguments[param.key]) == 40
                ):
                    value = provided_arguments[param.key]
                else:
                    logging.debug(
                        f"No argument provided for needed parameter {param.key} "
                        f"on component {component.name}. "
                        f"Substituting with {param.default}.")
                    value = param.default

                if param.key.endswith("_REF") and param.repository and value:
                    logging.debug(
                        f"The parameter {param.key} on component {component.name} "
                        f"points to the repository {param.repository}.")
                    # If a commit has already been resolved, it will be propagated to
                    # the next test in the test suite. To avoid resolving the same
                    # commit again, check if the value has the same length as the
                    # output of _resolve_branch_ref_to_commit.
                    cache_key = (str(param.repository), value)
                    if len(value) == 40:
                        logging.debug(
                            f"Git ref {value} is 40 characters long and probably already a commit.")
                    elif cache_key in resolved_ref_cache:
                        value = resolved_ref_cache[cache_key]
                    else:
                        value = self._resolve_branch_ref_to_commit(
                            param.repository, value)
                        resolved_ref_cache[cache_key] = value
                    if param.key in explicit_cli_keys:
                        provided_arguments[param.key] = value

                component_args[param.key] = value

            self.build_arguments_by_component[component.name] = component_args

        self.params_to_use = self._merged_build_arguments()

    def _primary_case(self) -> Case:
        return self.case_combination.cases[0]

    def _primary_platform(self) -> str | None:
        return self.build_arguments_by_component[
            self._primary_case().component.name
        ].get("PLATFORM")

    def _dockerfile_context_relative(self, platform: str) -> Path:
        return Path("..") / "tests" / "dockerfiles" / Path(platform)

    def _merged_build_arguments(self) -> Dict[str, str]:
        merged: Dict[str, str] = {}
        for case in self.case_combination.cases:
            component_args = self.build_arguments_by_component[case.component.name]
            for key, value in component_args.items():
                if key in merged and merged[key] != value:
                    logging.warning(
                        f"Build argument '{key}' differs between components "
                        f"({merged[key]!r} vs {value!r}); using {merged[key]!r} "
                        f"for shared build stages.")
                elif key not in merged:
                    merged[key] = value
        return merged

    def _validate_dockerfile_platform(self, platform: str, context: str = "") -> None:
        # Use an absolute path here only for validation that the requested
        # dockerfile context exists on the machine running the system tests.
        dockerfile_context = PRECICE_TESTS_DIR / "dockerfiles" / Path(platform)
        if not dockerfile_context.exists():
            suffix = f" for {context}" if context else ""
            raise ValueError(
                f"The path {dockerfile_context.resolve()} resulting from argument "
                f"PLATFORM={platform}{suffix} could not be found in the system")

    def __get_docker_services(self) -> Dict[str, str]:
        """
        Renders the service templates for each case using the parameters to use.

        Returns:
            A dictionary of rendered services per case name.
        """
        seen_platforms: set[str] = set()
        for case in self.case_combination.cases:
            platform = self.build_arguments_by_component[case.component.name].get(
                "PLATFORM")
            if platform and platform not in seen_platforms:
                self._validate_dockerfile_platform(
                    platform, f"component {case.component.name}")
                seen_platforms.add(platform)

        primary_platform = self._primary_platform()
        if not primary_platform:
            raise KeyError("Please specify a PLATFORM argument")
        self.dockerfile_context = (
            PRECICE_TESTS_DIR / "dockerfiles" / Path(primary_platform)
        )

        def render_service_template_per_case(case: Case) -> str:
            build_arguments = self.build_arguments_by_component[case.component.name]
            platform = build_arguments.get("PLATFORM")
            if not platform:
                raise KeyError(
                    f"Please specify a PLATFORM argument for component "
                    f"{case.component.name}")

            # Inside the individual system test directory (`self.system_test_dir`)
            # we copy a full `tests/` tree into the parent run directory
            # (see __copy_tools_and_tests). From the point of view of the system test
            # directory we therefore need to go one level up to reach the
            # shared `tests/` folder:
            #   <run_directory>/tests/dockerfiles/<PLATFORM>
            #   ^-------------^ parent of self.system_test_dir
            dockerfile_context_relative = self._dockerfile_context_relative(platform)

            render_dict = {
                # Use a relative path to the *parent* run directory so that
                # containers still see /runs/<tutorial_folder> like before,
                # while keeping the compose file independent of the CI
                # runner's absolute paths.
                'run_directory': "..",
                'tutorial_folder': self.tutorial_folder,
                'build_arguments': build_arguments,
                'params': build_arguments,
                'case_folder': case.path,
                'run': case.run_cmd,
                'dockerfile_context': dockerfile_context_relative,
            }
            jinja_env = Environment(loader=FileSystemLoader(PRECICE_TESTS_DIR))
            template = jinja_env.get_template(case.component.template)
            return template.render(render_dict)

        rendered_services = {}
        for case in self.case_combination.cases:
            rendered_services[case.name] = render_service_template_per_case(case)
        return rendered_services

    def __get_docker_compose_file(self):
        rendered_services = self.__get_docker_services()
        render_dict = {
            # See __get_docker_services: keep the docker-compose file
            # portable by referring to the parent run directory only.
            'run_directory': "..",
            'tutorial_folder': self.tutorial_folder,
            'tutorial': self.tutorial.path.name,
            'services': rendered_services,
            'build_arguments': self.params_to_use,
            # The dockerfile_context value inside the templates is only
            # used as a build context path and does not need to be
            # absolute – it will be resolved relative to the system test
            # directory.
            'dockerfile_context': self._dockerfile_context_relative(
                self._primary_platform()),
            'precice_output_folder': PRECICE_REL_OUTPUT_DIR,
        }
        jinja_env = Environment(loader=FileSystemLoader(PRECICE_TESTS_DIR))
        template = jinja_env.get_template("docker-compose.template.yaml")
        return template.render(render_dict)

    def __get_field_compare_compose_file(self):
        render_dict = {
            # Fieldcompare should also use only relative paths from inside
            # the system test directory so that the run directory can be
            # moved and re-executed elsewhere.
            'run_directory': "..",
            'tutorial_folder': self.tutorial_folder,
            'precice_output_folder': PRECICE_REL_OUTPUT_DIR,
            'reference_output_folder': PRECICE_REL_REFERENCE_DIR + "/" + self.reference_result.path.name.replace(".tar.gz", ""),
            'tolerance': self.tolerance,
            # The dockerfile_context value inside the templates is only
            # used as a build context path and does not need to be
            # absolute – it will be resolved relative to the system test
            # directory.
            'dockerfile_context': self._dockerfile_context_relative(
                self._primary_platform()),
        }
        jinja_env = Environment(loader=FileSystemLoader(PRECICE_TESTS_DIR))
        template = jinja_env.get_template(
            "docker-compose.field_compare.template.yaml")
        return template.render(render_dict)

    def _get_git_ref(self, repository: Path, abbrev_ref=False) -> Optional[str]:
        try:
            result = subprocess.run([
                "git",
                "-C", os.fspath(repository.resolve()),
                "rev-parse",
                "--abbrev-ref" if abbrev_ref else
                "HEAD"], stdout=subprocess.PIPE,
                stderr=subprocess.PIPE, text=True, check=True, timeout=60)
            current_ref = result.stdout.strip()
            return current_ref
        except Exception as e:
            raise RuntimeError(f"An error occurred while getting the current Git ref: {e}") from e

    def _fetch_pr(self, repository: Path, pr: str):
        try:
            result = subprocess.run([
                "git",
                "-C", os.fspath(repository.resolve()),
                "fetch",
                "origin",
                "pull/" + pr + "/head"
            ], check=True, timeout=60)

        except Exception as e:
            raise RuntimeError(f"git command returned code {result.returncode}")

    def _fetch_ref(self, repository: Path, ref: str):
        try:
            result = subprocess.run([
                "git",
                "-C", os.fspath(repository.resolve()),
                "fetch"
            ], check=True, timeout=60)
            if result.returncode != 0:
                raise RuntimeError(f"git command returned code {result.returncode}")

        except Exception as e:
            raise RuntimeError(
                f"An error occurred while fetching origin '{ref}':  {e}. Do the values in reference_versions.yaml point to (still) valid Git refs?")

    def _resolve_branch_ref_to_commit(self, repository: Path, ref: str) -> Optional[str]:
        try:
            git_ls_remote_output = subprocess.run([
                "git",
                "ls-remote",
                os.fspath(repository),
                ref,
            ], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True, timeout=60)

            # If an invalid ref is given, git ls-remote still returns success, but no list
            git_remote_refs = git_ls_remote_output.stdout.strip()
            if not git_remote_refs:
                raise ValueError(f"The git ref {ref} does not appear in the repository {repository}.")

            commit = git_remote_refs.split()[0]
            # The output assumes a URL of the form <repository>/commits/<commit>. Works for GitHub and Bitbucket.
            logging.debug(
                f"Resolved the git ref {ref} of the repository {repository} to {repository}/commits/{commit} .")
            return commit if commit else ref
        except Exception:
            logging.warning(
                f"Could not resolve git ref {ref} of the repository {repository} to a commit. Using the given git ref as-is.")
            return ref

    def _checkout_ref_in_subfolder(self, repository: Path, subfolder: Path, ref: str):
        try:
            result = subprocess.run([
                "git",
                "-C", os.fspath(repository.resolve()),
                "checkout", ref,
                "--", os.fspath(subfolder.resolve())
            ], check=True, timeout=60)
            if result.returncode != 0:
                raise RuntimeError(f"git command returned code {result.returncode}")

        except Exception as e:
            raise RuntimeError(f"An error occurred while checking out '{ref}' for folder '{repository}': {e}")

    def __copy_tutorial_into_directory(self, run_directory: Path):
        """
        Checks out the requested tutorial ref and copies the entire tutorial into a folder to prepare for running.
        """
        current_time_string = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        self.run_directory = run_directory
        current_ref = None
        ref_requested = None

        if self.tutorial.source.type == "local":
            pr_requested = self.params_to_use.get("TUTORIALS_PR")
            if pr_requested:
                logging.debug(f"Fetching the PR {pr_requested} HEAD reference")
                self._fetch_pr(PRECICE_TUTORIAL_DIR, pr_requested)
            current_ref = self._get_git_ref(PRECICE_TUTORIAL_DIR)
            ref_requested = self.params_to_use.get("TUTORIALS_REF")
            if ref_requested:
                logging.debug(f"Checking out tutorials {ref_requested} before copying")
                self._fetch_ref(PRECICE_TUTORIAL_DIR, ref_requested)
                self._checkout_ref_in_subfolder(
                    PRECICE_TUTORIAL_DIR, self.tutorial.path, ref_requested)

        self.tutorial_folder = slugify(
            f'{self.tutorial.path.name}_{self.case_combination.cases}_{current_time_string}')
        destination = run_directory / self.tutorial_folder
        # External sources are fetched and resolved once at parse time; reuse
        # that path here to avoid a redundant fetch (and duplicate log line).
        src = self.tutorial.resolved_root or resolve_tutorial_root(
            self.tutorial.path,
            self.tutorial.source,
            PRECICE_EXTERNAL_CACHE_DIR,
        )
        self.system_test_dir = destination
        shutil.copytree(src, destination)

        if self.tutorial.source.type == "local" and ref_requested:
            with open(destination / "tutorials_ref", 'w') as file:
                file.write(ref_requested)
            self._checkout_ref_in_subfolder(PRECICE_TUTORIAL_DIR, self.tutorial.path, current_ref)

    def __copy_tools_and_tests(self, run_directory: Path):
        src = PRECICE_TOOLS_DIR
        destination = run_directory / "tools"
        logging.debug(f"Copying tools from {src} to {destination}")
        try:
            shutil.copytree(src, destination)
        except FileExistsError as e:
            logging.debug(f"Tools directory has already been copied to the workspace - skipping.")
        except Exception as e:
            logging.warning(f"Something went wrong while copying the tools directory to the workspace: {e}")

        src = PRECICE_TESTS_DIR
        destination = run_directory / "tests"
        logging.debug(f"Copying tests from {src} to {destination}")
        try:
            shutil.copytree(src, destination)
        except FileExistsError as e:
            logging.debug(f"Tests directory has already been copied to the workspace - skipping.")
        except Exception as e:
            logging.warning(f"Something went wrong while copying the tests directory to the workspace: {e}")

    def __put_gitignore(self, run_directory: Path):
        # Create the .gitignore file with a single asterisk
        gitignore_file = run_directory / ".gitignore"
        with gitignore_file.open("w") as file:
            file.write("*")

    def __cleanup(self):
        shutil.rmtree(self.run_directory)

    def __get_uid_gid(self):
        try:
            uid = int(subprocess.check_output(["id", "-u"]).strip())
            gid = int(subprocess.check_output(["id", "-g"]).strip())
            return uid, gid
        except Exception as e:
            logging.error("Error getting group and user id: ", e)

    def __write_env_file(self):
        with open(self.system_test_dir / ".env", "w") as env_file:
            for key, value in self.env.items():
                env_file.write(f"{key}={value}\n")

    def __unpack_reference_results(self) -> Tuple[bool, str]:
        if not self.reference_result.path.exists():
            error_message = (
                f"Reference results archive was not found for {self}. "
                f"Expected file: {self.reference_result.path}. "
                "Please generate the reference results first or update tests.yaml accordingly.")
            logging.error(error_message)
            return False, error_message

        try:
            # Base directory where reference results should be extracted
            dest_dir = self.system_test_dir / PRECICE_REL_REFERENCE_DIR
            dest_dir.mkdir(parents=True, exist_ok=True)
            dest_dir_resolved = dest_dir.resolve()

            with tarfile.open(self.reference_result.path) as reference_results_tared:
                # Validate that each member will be extracted within dest_dir
                for member in reference_results_tared.getmembers():
                    member_path = dest_dir / member.name
                    member_path_resolved = member_path.resolve()
                    # Ensure the resolved member path is within the destination directory
                    if os.path.commonpath([str(dest_dir_resolved), str(
                            member_path_resolved)]) != str(dest_dir_resolved):
                        logging.error(
                            f"Unsafe path detected in reference results archive {self.reference_result.path} "
                            f"for {self}: {member.name}")
                        return False

                # All paths are safe; extract into the destination directory
                reference_results_tared.extractall(dest_dir)

            logging.debug(
                f"extracting {self.reference_result.path} into {dest_dir}")
            return True, ""
        except (tarfile.TarError, OSError) as e:
            error_message = (
                f"Could not unpack reference results archive {self.reference_result.path} for {self}: {e}")
            logging.error(error_message)
            return False, error_message

    def __init_run_logs(self) -> None:
        self._log_sink = _SystemtestLogSink(self.system_test_dir)

    def _run_docker_compose_subprocess(
        self,
        command: List[str],
        stage: str,
        timeout: int,
    ) -> Tuple[int, List[str], List[str]]:
        """
        Run a docker compose command, streaming stdout/stderr to log files as they arrive.
        """
        stdout_data: List[str] = []
        stderr_data: List[str] = []
        log_sink = getattr(self, "_log_sink", None)
        if log_sink is not None:
            log_sink.begin_stage(stage)
        logging.info(f"Docker compose {stage} for {self}")

        try:
            process = subprocess.Popen(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                bufsize=1,
                start_new_session=True,
                cwd=self.system_test_dir,
            )
        except Exception as e:
            logging.critical(f"Error starting docker compose {stage} command: {e}")
            return 1, stdout_data, stderr_data

        def read_stream(stream, is_stderr: bool) -> None:
            if stream is None:
                return
            for line in stream:
                line = line.rstrip("\n\r")
                if is_stderr:
                    stderr_data.append(line)
                    if log_sink is not None:
                        log_sink.append_stderr(line, stage)
                else:
                    stdout_data.append(line)
                    if log_sink is not None:
                        log_sink.append_stdout(line, stage)
            stream.close()

        stdout_thread = threading.Thread(
            target=read_stream, args=(process.stdout, False), daemon=True)
        stderr_thread = threading.Thread(
            target=read_stream, args=(process.stderr, True), daemon=True)
        stdout_thread.start()
        stderr_thread.start()

        try:
            exit_code = process.wait(timeout=timeout)
        except KeyboardInterrupt as k:
            process.kill()
            stdout_thread.join(timeout=SHORT_TIMEOUT)
            stderr_thread.join(timeout=SHORT_TIMEOUT)
            raise KeyboardInterrupt from k
        except subprocess.TimeoutExpired:
            logging.critical(
                f"Systemtest {self} timed out during docker compose {stage} "
                f"after {timeout}s. Killing the process.")
            process.kill()
            try:
                process.wait(timeout=SHORT_TIMEOUT)
            except subprocess.TimeoutExpired:
                pass
            exit_code = process.returncode if process.returncode is not None else 1
        except Exception as e:
            logging.critical(
                f"Systemtest {self} had serious issues during docker compose {stage}: {e}")
            process.kill()
            try:
                process.wait(timeout=SHORT_TIMEOUT)
            except subprocess.TimeoutExpired:
                pass
            exit_code = process.returncode if process.returncode is not None else 1

        stdout_thread.join(timeout=SHORT_TIMEOUT)
        stderr_thread.join(timeout=SHORT_TIMEOUT)
        if exit_code is None:
            exit_code = process.poll() or 1
        return exit_code, stdout_data, stderr_data

    def _cleanup_docker_networks(self):
        """
        Prunes the unused Docker networks, since there is an upper limit on the number of custom networks defined.
        """
        logging.debug(f"Deleting unused Docker networks...")
        stdout_data = []
        stderr_data = []
        try:
            # Execute docker-network-prune command
            process = subprocess.Popen(['docker',
                                        'network',
                                        'prune',
                                        '-f'],
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE,
                                       start_new_session=True,
                                       cwd=self.system_test_dir)
            try:
                stdout, stderr = process.communicate(timeout=self.timeout)
            except KeyboardInterrupt as k:
                process.kill()
                raise KeyboardInterrupt from k
        except Exception as e:
            logging.critical(
                f"Systemtest {self} could not prune the Docker networks. This might prevent tests from starting.")
            stdout_data.extend(stdout.decode().splitlines())
            stderr_data.extend(stderr.decode().splitlines())
            process.poll()

    def _run_field_compare(self):
        """
        Executes the field comparison step after unpacking reference results.

        Returns:
            A FieldCompareResult object containing the command outcome and logs.
        """
        logging.debug(f"Running fieldcompare for {self}")
        time_start = time.perf_counter()
        unpack_success, unpack_error_message = self.__unpack_reference_results()
        if not unpack_success:
            log_sink = getattr(self, "_log_sink", None)
            if log_sink is not None:
                log_sink.begin_stage("compare")
                log_sink.append_stderr(unpack_error_message, "compare")
            elapsed_time = time.perf_counter() - time_start
            return FieldCompareResult(1, [], [unpack_error_message], self, elapsed_time)
        docker_compose_content = self.__get_field_compare_compose_file()

        with open(self.system_test_dir / "docker-compose.field_compare.yaml", 'w') as file:
            file.write(docker_compose_content)
        exit_code, stdout_data, stderr_data = self._run_docker_compose_subprocess(
            [
                'docker',
                'compose',
                '--file',
                'docker-compose.field_compare.yaml',
                'up',
                '--exit-code-from',
                'field-compare',
            ],
            "compare",
            self.timeout,
        )
        elapsed_time = time.perf_counter() - time_start
        return FieldCompareResult(exit_code, stdout_data, stderr_data, self, elapsed_time)

    def __archive_fieldcompare_diffs(self) -> None:
        """
        Copy fieldcompare diff VTK files from precice-exports/ into diff-results/,
        preserving paths under precice-exports/ so nested outputs are not skipped
        and identical basenames in different folders do not overwrite each other.
        """
        exports_dir = self.system_test_dir / PRECICE_REL_OUTPUT_DIR
        if not exports_dir.is_dir():
            return
        suffixes = (".vtu", ".vtk", ".vtp")
        dest_root = self.system_test_dir / DIFF_RESULTS_DIR
        seen_resolved: set[Path] = set()
        archived_count = 0
        for path in exports_dir.rglob("*"):
            if not path.is_file():
                continue
            if path.suffix.lower() not in suffixes:
                continue
            if "diff" not in path.name.lower():
                continue
            resolved = path.resolve()
            if resolved in seen_resolved:
                continue
            try:
                rel = path.relative_to(exports_dir)
            except ValueError:
                continue
            seen_resolved.add(resolved)
            dest_path = dest_root / rel
            dest_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, dest_path)
            archived_count += 1
        if archived_count:
            logging.debug(
                "Archived %d fieldcompare diff file(s) to %s for %s",
                archived_count,
                dest_root,
                self,
            )

    def __get_diff_visualizer_compose_file(self) -> str:
        platform = self.params_to_use.get("PLATFORM")
        render_dict = {
            'dockerfile_context': (
                Path("..") / "tests" / "dockerfiles" / Path(platform)
            ),
            'build_arguments': self.params_to_use,
            'diff_results_folder': DIFF_RESULTS_DIR,
        }
        jinja_env = Environment(loader=FileSystemLoader(PRECICE_TESTS_DIR))
        template = jinja_env.get_template(
            "docker-compose.diff_visualizer.template.yaml")
        return template.render(render_dict)

    def __append_diff_visualizer_status(self, status: str, elapsed_s: float) -> None:
        log_path = self.system_test_dir / DIFF_VISUALIZER_LOG
        with log_path.open("a", encoding="utf-8") as log_file:
            log_file.write(f"\nstatus: {status}\nelapsed_s: {elapsed_s:.1f}\n")

    def __visualize_fieldcompare_diffs(self) -> None:
        """Best-effort rendering of archived fieldcompare diff VTK files via Docker."""
        diff_results_dir = self.system_test_dir / DIFF_RESULTS_DIR
        if not diff_results_dir.is_dir():
            return

        compose_path = self.system_test_dir / "docker-compose.diff_visualizer.yaml"
        log_path = self.system_test_dir / DIFF_VISUALIZER_LOG
        log_path.write_text("=== compare-diff ===\n", encoding="utf-8")
        log_lock = threading.Lock()
        time_start = time.perf_counter()

        try:
            compose_path.write_text(
                self.__get_diff_visualizer_compose_file(), encoding="utf-8")
        except OSError as error:
            elapsed_s = time.perf_counter() - time_start
            self.__append_diff_visualizer_status(f"error: {error}", elapsed_s)
            logging.warning(
                "Could not render fieldcompare diff visualizations for %s: %s",
                self,
                error,
            )
            return

        logging.info(
            "Rendering fieldcompare diff visualizations for %s "
            "(timeout %ss)",
            self,
            DIFF_VISUALIZER_TIMEOUT,
        )
        try:
            process = subprocess.Popen(
                [
                    "docker",
                    "compose",
                    "--file",
                    compose_path.name,
                    "up",
                    "--exit-code-from",
                    "diff-visualizer",
                    "--abort-on-container-exit",
                ],
                cwd=self.system_test_dir,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                bufsize=1,
            )
        except OSError as error:
            elapsed_s = time.perf_counter() - time_start
            self.__append_diff_visualizer_status(f"error: {error}", elapsed_s)
            logging.warning(
                "Could not render fieldcompare diff visualizations for %s: %s",
                self,
                error,
            )
            return

        def read_stream(stream, prefix: str) -> None:
            if stream is None:
                return
            for line in stream:
                line = line.rstrip("\n\r")
                with log_lock:
                    with log_path.open("a", encoding="utf-8") as log_file:
                        log_file.write(f"{prefix}{line}\n")
            stream.close()

        stdout_thread = threading.Thread(
            target=read_stream, args=(process.stdout, ""), daemon=True)
        stderr_thread = threading.Thread(
            target=read_stream, args=(process.stderr, "[stderr] "), daemon=True)
        stdout_thread.start()
        stderr_thread.start()

        timed_out = False
        try:
            exit_code = process.wait(timeout=DIFF_VISUALIZER_TIMEOUT)
        except subprocess.TimeoutExpired:
            timed_out = True
            process.kill()
            try:
                process.wait(timeout=SHORT_TIMEOUT)
            except subprocess.TimeoutExpired:
                pass
            exit_code = process.returncode if process.returncode is not None else 1

        stdout_thread.join(timeout=SHORT_TIMEOUT)
        stderr_thread.join(timeout=SHORT_TIMEOUT)
        elapsed_s = time.perf_counter() - time_start

        if timed_out:
            self.__append_diff_visualizer_status(
                f"timed out after {DIFF_VISUALIZER_TIMEOUT}s", elapsed_s
            )
            logging.warning(
                "Could not render fieldcompare diff visualizations for %s: "
                "timed out after %ss (visualizer ran %.1fs). "
                "See %s",
                self,
                DIFF_VISUALIZER_TIMEOUT,
                elapsed_s,
                DIFF_VISUALIZER_LOG,
            )
            return

        if exit_code != 0:
            self.__append_diff_visualizer_status(
                f"failed (exit {exit_code})", elapsed_s
            )
            logging.warning(
                "Rendering fieldcompare diff visualizations failed for %s "
                "after %.1fs (exit %s). See %s",
                self,
                elapsed_s,
                exit_code,
                DIFF_VISUALIZER_LOG,
            )
            return

        self.__append_diff_visualizer_status("ok", elapsed_s)
        logging.info(
            "Diff visualizations for %s took %.1fs",
            self,
            elapsed_s,
        )

    def __copy_rerun_system_test_script(self) -> None:
        """Copy tests/rerun-system-test.sh into the run directory for artifact replay."""
        rerun_src = PRECICE_TESTS_DIR / "rerun-system-test.sh"
        if not rerun_src.is_file():
            raise FileNotFoundError(
                f"Missing {rerun_src}. It is required for portable CI artifact replay.")
        rerun_dst = self.system_test_dir / "rerun-system-test.sh"
        shutil.copy2(rerun_src, rerun_dst)
        try:
            rerun_dst.chmod(rerun_dst.stat().st_mode | 0o111)
        except Exception:
            logging.debug(
                f"Could not mark {rerun_dst} as executable; continuing anyway.")

    @staticmethod
    def _sha256_file(path: Path) -> str:
        """Compute SHA-256 hex digest of a file."""
        h = hashlib.sha256()
        mv = memoryview(bytearray(128 * 1024))
        with open(path, 'rb', buffering=0) as f:
            while n := f.readinto(mv):
                h.update(mv[:n])
        return h.hexdigest()

    def _unpacked_reference_iterations_logs_dir(self) -> Path:
        """Iterations logs unpacked from the reference tar (fieldcompare runs first)."""
        stem = self.reference_result.path.name.replace(".tar.gz", "")
        return (
            self.system_test_dir
            / PRECICE_REL_REFERENCE_DIR
            / f"{stem}.{ITERATIONS_LOGS_DIR}"
        )

    def _collect_iterations_logs(
        self, system_test_dir: Path
    ) -> List[Tuple[str, Path]]:
        """
        Collect precice-*-iterations.log files from case dirs.
        Returns list of (relative_path, absolute_path) e.g. ("solid-fenics/precice-Solid-iterations.log", path).
        """
        collected = []
        for case in self.case_combination.cases:
            case_dir = system_test_dir / Path(case.path).name
            if not case_dir.exists():
                continue
            for log_file in case_dir.glob("precice-*-iterations.log"):
                if log_file.is_file():
                    rel = f"{Path(case.path).name}/{log_file.name}"
                    collected.append((rel, log_file))
        return collected

    def _reference_iterations_hashes(self) -> Optional[Dict[str, str]]:
        """
        Load expected iterations.log hashes from archived reference files.
        Returns None if no reference data is available.
        """
        ref_dir = self._unpacked_reference_iterations_logs_dir()
        if not ref_dir.is_dir():
            return None
        ref_hashes = {}
        for log_file in ref_dir.rglob("precice-*-iterations.log"):
            if log_file.is_file():
                rel = log_file.relative_to(ref_dir).as_posix()
                ref_hashes[rel] = self._sha256_file(log_file)
        return ref_hashes if ref_hashes else None

    def __archive_iterations_logs(self) -> None:
        """Copy precice-*-iterations.log from case dirs into iterations-logs/ for CI artifacts."""
        collected = self._collect_iterations_logs(self.system_test_dir)
        if not collected:
            return
        dest_dir = self.system_test_dir / ITERATIONS_LOGS_DIR
        dest_dir.mkdir(exist_ok=True)
        for rel, src in collected:
            dest_name = Path(rel).name
            if len(collected) > 1:
                prefix = Path(rel).parent.name + "_"
                dest_name = prefix + dest_name
            shutil.copy2(src, dest_dir / dest_name)
        logging.debug(
            "Archived %d iterations log(s) to %s for %s",
            len(collected),
            dest_dir,
            self,
        )

    def _append_compare_log(self, message: str, *, error: bool = False) -> None:
        log_sink = getattr(self, "_log_sink", None)
        if log_sink is None:
            return
        if error:
            log_sink.append_stderr(message, "compare")
        else:
            log_sink.append_stdout(message, "compare")

    def __compare_iterations_hashes(self) -> bool:
        """
        Compare current iterations.log hashes against reference data.
        Returns True if comparison passes (or is skipped). Returns False if hashes differ.
        """
        ref_hashes = self._reference_iterations_hashes()
        if ref_hashes is None:
            message = (
                f"Iterations.log hash check skipped (no reference data) for {self}"
            )
            logging.info(message)
            self._append_compare_log(message)
            return True
        collected = self._collect_iterations_logs(self.system_test_dir)
        current = {rel: self._sha256_file(p) for rel, p in collected}
        for rel, expected in ref_hashes.items():
            if rel not in current:
                message = (
                    f"Missing iterations log {rel} (expected from reference); "
                    f"{self} fails"
                )
                logging.critical(message)
                self._append_compare_log(message, error=True)
                return False
            if current[rel] != expected:
                message = (
                    f"Hash mismatch for {rel} (iterations.log regression); "
                    f"{self} fails"
                )
                logging.critical(message)
                self._append_compare_log(message, error=True)
                return False
        if len(current) != len(ref_hashes):
            extra = set(current) - set(ref_hashes)
            message = f"Unexpected iterations log(s) {extra}; {self} fails"
            logging.critical(message)
            self._append_compare_log(message, error=True)
            return False
        self._append_compare_log("=== Comparing the iterations.log files (checksums only) ===")
        for rel in sorted(ref_hashes):
            detail = f"  {rel}: sha256 ok"
            logging.debug(detail)
            self._append_compare_log(detail)
        message = (
            f"Iterations.log hash check passed for {self} ({len(ref_hashes)} file(s))"
        )
        logging.info(message)
        self._append_compare_log(message)
        return True

    def _build_docker(self):
        """
        Builds the docker image
        """
        logging.debug(f"Building docker image for {self}")
        logging.info(
            f"Using build timeout {self.build_timeout}s for {self}")
        time_start = time.perf_counter()
        docker_compose_content = self.__get_docker_compose_file()
        docker_compose_path = self.system_test_dir / "docker-compose.tutorial.yaml"
        with open(docker_compose_path, 'w') as file:
            file.write(docker_compose_content)

        field_compare_compose_path = self.system_test_dir / "docker-compose.field_compare.yaml"
        with open(field_compare_compose_path, 'w') as file:
            file.write(self.__get_field_compare_compose_file())

        self.__copy_rerun_system_test_script()

        exit_code, stdout_data, stderr_data = self._run_docker_compose_subprocess(
            [
                'docker',
                'compose',
                '--progress=plain',
                '--file',
                'docker-compose.tutorial.yaml',
                'build',
            ],
            "build",
            self.build_timeout,
        )
        elapsed_time = time.perf_counter() - time_start
        return DockerComposeResult(exit_code, stdout_data, stderr_data, self, elapsed_time)

    def _run_tutorial(self):
        """
        Runs precice couple

        Returns:
            A DockerComposeResult object containing the state.
        """
        logging.debug(f"Running tutorial {self}")
        time_start = time.perf_counter()
        exit_code, stdout_data, stderr_data = self._run_docker_compose_subprocess(
            [
                'docker',
                'compose',
                '--file',
                'docker-compose.tutorial.yaml',
                'up',
            ],
            "run",
            self.timeout,
        )
        elapsed_time = time.perf_counter() - time_start
        return DockerComposeResult(exit_code, stdout_data, stderr_data, self, elapsed_time)

    def __repr__(self):
        prefix = "External: " if getattr(self.tutorial.source, "type", "local") != "local" else ""
        return f"{prefix}{self.tutorial.name} {self.case_combination}"

    def __apply_max_time_override(self):
        """Overwrite <max-time> or <max-time-windows> value in precice-config.xml."""
        if self.max_time is None and self.max_time_windows is None:
            return
        config_path = self.system_test_dir / "precice-config.xml"
        text = config_path.read_text()
        new_text = text
        if self.max_time is not None:
            pattern = r'(<max-time\s+value=")[^"]*(")'
            new_text, count = re.subn(pattern, rf'\g<1>{self.max_time}\2', new_text)
            if count == 0:
                logging.warning(f"No <max-time> tag found in {config_path}")
            else:
                logging.info(f"Overwrote <max-time> to {self.max_time} in {config_path}")
        if self.max_time_windows is not None:
            pattern = r'(<max-time-windows\s+value=")[^"]*(")'
            new_text, count = re.subn(pattern, rf'\g<1>{self.max_time_windows}\2', new_text)
            if count == 0:
                logging.warning(f"No <max-time-windows> tag found in {config_path}")
            else:
                logging.info(f"Overwrote <max-time-windows> to {self.max_time_windows} in {config_path}")
        config_path.write_text(new_text)

    def _run_hook(self, stage: str, command: str | None) -> bool:
        """
        Run a shell command in the copied tutorial directory (e.g. run-before / run-after).
        """
        if not command:
            return True
        logging.info(f"Running {stage} for {self}: {command}")
        try:
            result = subprocess.run(
                command,
                shell=True,
                cwd=self.system_test_dir,
                capture_output=True,
                text=True,
                start_new_session=True,
            )
        except Exception as e:
            logging.critical(f"Failed to start {stage} for {self}: {e}")
            return False
        hook_output = (result.stdout or '') + (result.stderr or '')
        if hook_output.strip():
            logging.debug(f"{stage} output for {self}:\n{hook_output.rstrip()}")
        if result.returncode != 0:
            logging.critical(
                f"{stage} for {self} failed with exit code {result.returncode}: {command}")
            return False
        return True

    def __prepare_for_run(self, run_directory: Path):
        """
        Prepares the run_directory with folders and datastructures needed for every systemtest execution
        """
        self.__copy_tutorial_into_directory(run_directory)
        if not self._run_hook('run-before', self.run_before):
            raise RuntimeError(f"run-before hook failed for {self}")
        self.__apply_max_time_override()
        self.__copy_tools_and_tests(run_directory)
        self.__put_gitignore(run_directory)
        host_uid, host_gid = self.__get_uid_gid()
        self.params_to_use['PRECICE_UID'] = host_uid
        self.params_to_use['PRECICE_GID'] = host_gid
        for component_args in self.build_arguments_by_component.values():
            component_args['PRECICE_UID'] = host_uid
            component_args['PRECICE_GID'] = host_gid

    def run(self, run_directory: Path):
        """
        Runs the system test by generating the Docker Compose file, copying everything into a run folder, and executing docker-compose up.
        """
        try:
            self.__prepare_for_run(run_directory)
        except RuntimeError as e:
            logging.critical(str(e))
            return SystemtestResult(False, [], [str(e)], self, build_time=0, solver_time=0, fieldcompare_time=0)

        self.__init_run_logs()
        std_out: List[str] = []
        std_err: List[str] = []

        self._cleanup_docker_networks()
        docker_build_result = self._build_docker()
        std_out.extend(docker_build_result.stdout_data)
        std_err.extend(docker_build_result.stderr_data)
        if docker_build_result.exit_code != 0:
            logging.critical(f"Could not build the docker images, {self} failed")
            return SystemtestResult(
                False,
                std_out,
                std_err,
                self,
                build_time=docker_build_result.runtime,
                solver_time=0,
                fieldcompare_time=0)

        docker_run_result = self._run_tutorial()
        std_out.extend(docker_run_result.stdout_data)
        std_err.extend(docker_run_result.stderr_data)
        if docker_run_result.exit_code != 0:
            logging.critical(f"Could not run the tutorial, {self} failed")
            return SystemtestResult(
                False,
                std_out,
                std_err,
                self,
                build_time=docker_build_result.runtime,
                solver_time=docker_run_result.runtime,
                fieldcompare_time=0)

        if not self._run_hook('run-after', self.run_after):
            logging.critical(f"run-after hook failed for {self}")
            return SystemtestResult(
                False,
                std_out,
                std_err,
                self,
                build_time=docker_build_result.runtime,
                solver_time=docker_run_result.runtime,
                fieldcompare_time=0)

        if self.skip_compare:
            logging.info(f"Skipping fieldcompare for {self} (skip_compare=true)")
            fieldcompare_time = 0.0
        else:
            fieldcompare_result = self._run_field_compare()
            std_out.extend(fieldcompare_result.stdout_data)
            std_err.extend(fieldcompare_result.stderr_data)
            if fieldcompare_result.exit_code != 0:
                self.__archive_fieldcompare_diffs()
                self.__visualize_fieldcompare_diffs()
                logging.critical(f"Fieldcompare returned non zero exit code, therefore {self} failed")
                return SystemtestResult(
                    False,
                    std_out,
                    std_err,
                    self,
                    build_time=docker_build_result.runtime,
                    solver_time=docker_run_result.runtime,
                    fieldcompare_time=fieldcompare_result.runtime)
            fieldcompare_time = fieldcompare_result.runtime

        self.__archive_iterations_logs()
        if not self.__compare_iterations_hashes():
            logging.critical(
                f"Iterations.log hash comparison failed (regression), {self} failed"
            )
            return SystemtestResult(
                False,
                std_out,
                std_err,
                self,
                build_time=docker_build_result.runtime,
                solver_time=docker_run_result.runtime,
                fieldcompare_time=fieldcompare_time)

        # self.__cleanup()
        self._cleanup_docker_networks()
        return SystemtestResult(
            True,
            std_out,
            std_err,
            self,
            build_time=docker_build_result.runtime,
            solver_time=docker_run_result.runtime,
            fieldcompare_time=fieldcompare_time)

    def run_for_reference_results(self, run_directory: Path):
        """
        Runs the system test by generating the Docker Compose files to generate the reference results
        """
        try:
            self.__prepare_for_run(run_directory)
        except RuntimeError as e:
            logging.critical(str(e))
            return SystemtestResult(False, [], [str(e)], self, build_time=0, solver_time=0, fieldcompare_time=0)

        self.__init_run_logs()
        std_out: List[str] = []
        std_err: List[str] = []
        self._cleanup_docker_networks()
        docker_build_result = self._build_docker()
        std_out.extend(docker_build_result.stdout_data)
        std_err.extend(docker_build_result.stderr_data)
        if docker_build_result.exit_code != 0:
            logging.critical(f"Could not build the docker images, {self} failed")
            return SystemtestResult(
                False,
                std_out,
                std_err,
                self,
                build_time=docker_build_result.runtime,
                solver_time=0,
                fieldcompare_time=0)

        docker_run_result = self._run_tutorial()
        std_out.extend(docker_run_result.stdout_data)
        std_err.extend(docker_run_result.stderr_data)
        if docker_run_result.exit_code != 0:
            logging.critical(f"Could not run the tutorial, {self} failed")
            return SystemtestResult(
                False,
                std_out,
                std_err,
                self,
                build_time=docker_build_result.runtime,
                solver_time=docker_run_result.runtime,
                fieldcompare_time=0)

        if not self._run_hook('run-after', self.run_after):
            logging.critical(f"run-after hook failed for {self}")
            return SystemtestResult(
                False,
                std_out,
                std_err,
                self,
                build_time=docker_build_result.runtime,
                solver_time=docker_run_result.runtime,
                fieldcompare_time=0)

        self._cleanup_docker_networks()
        return SystemtestResult(
            True,
            std_out,
            std_err,
            self,
            build_time=docker_build_result.runtime,
            solver_time=docker_run_result.runtime,
            fieldcompare_time=0)

    def run_only_build(self, run_directory: Path):
        """
        Runs only the build commmand, for example to preheat the caches of the docker builder.
        """
        try:
            self.__prepare_for_run(run_directory)
        except RuntimeError as e:
            logging.critical(str(e))
            return SystemtestResult(False, [], [str(e)], self, build_time=0, solver_time=0, fieldcompare_time=0)

        self.__init_run_logs()
        std_out: List[str] = []
        std_err: List[str] = []
        docker_build_result = self._build_docker()
        std_out.extend(docker_build_result.stdout_data)
        std_err.extend(docker_build_result.stderr_data)
        if docker_build_result.exit_code != 0:
            logging.critical(f"Could not build the docker images, {self} failed")
            return SystemtestResult(
                False,
                std_out,
                std_err,
                self,
                build_time=docker_build_result.runtime,
                solver_time=0,
                fieldcompare_time=0)

        return SystemtestResult(
            True,
            std_out,
            std_err,
            self,
            build_time=docker_build_result.runtime,
            solver_time=0,
            fieldcompare_time=0)

    def get_system_test_dir(self) -> Path:
        return self.system_test_dir

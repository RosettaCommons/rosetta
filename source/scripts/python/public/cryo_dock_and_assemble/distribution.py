"""
Tools for distributing jobs
"""

import pathlib
import shutil
import os
import typing
import zlib
import subprocess
import tempfile
import glob


def trim_list_via_binary(l1: typing.List[typing.Any], l2: typing.List[bool]) -> typing.List[typing.Any]:
    if len(l1) != len(l2):
        raise RuntimeError("Trim list via binary list sizes were wrong!", len(l1), len(l2))
    ret_list = []
    for i, keep in enumerate(l2):
        if keep:
            ret_list.append(l1[i])
    return ret_list


def setup_temp_workdir(workdir: str) -> str:
    """
    Return orignal(current) directory, but switch workdir
    """
    og_dir = os.getcwd()
    if os.path.exists(workdir):
        shutil.rmtree(workdir, ignore_errors=True)
    pathlib.Path(workdir).mkdir(parents=True, exist_ok=False)
    os.chdir(workdir)
    return og_dir


class TmpDir:
    def __init__(self):
        pass

    def __enter__(self):
        self.og_dir = os.getcwd()
        self.T_D = tempfile.TemporaryDirectory()
        os.chdir(self.T_D.name)
        return self.T_D.name

    def __exit__(self, *args, **kwargs):
        os.chdir(self.og_dir)
        try:
            self.T_D.cleanup()
        except FileNotFoundError:
            pass


def run_subprocess(
    runcommand: typing.List[str],
    expected_outputs: typing.List[str],
    input_files: typing.Dict[str, bytes],
    task_id: str,
    max_tries: int = 2,
    save_broken_dir: str = "~/Databases/broken_runs/",
    save_logs: bool = False,
    dependencies: typing.Any = None,
    env: typing.Optional[typing.Dict[str, str]] = None,
) -> typing.Dict[str, bytes]:
    """
    Notes:
        pass in task_id because if we make task_id here in a multiprocessing context
        we can sometimes have the same seed.
    :param runcommand: what to run as a subprocess
    :param expected_outputs: Files that you want to collect once you're done with your computation
    :param input_files: Zlib compressed files that are required for computation (filename, bytes)
    :param task_id: A unique ID for this task
    :param max_tries: The maximum number of tries you want to perform before quitting/failing
    :param save_broken_dir: If set, will dump any commands that fail, and their directories to that location
    :param save_logs: If set we will save the logs to a file, and return those logs as 'out.log'
    """
    current_try = 0
    return_files = {}
    errors = []
    success = False
    with TmpDir():
        for input_filename, input_file_content in input_files.items():
            with open(input_filename, "wb") as fh:
                fh.write(zlib.decompress(input_file_content))
        while not success:
            current_try += 1
            if current_try > max_tries:
                break
            if save_logs:
                with open("out.log", "w") as output:
                    p = subprocess.call(runcommand, stdout=output, stderr=subprocess.STDOUT, env=env)
            else:
                p = subprocess.call(runcommand, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, env=env)
            if p:
                error = "Command: " + " ".join(runcommand) + f"\nreturned {p} in {os.getcwd()}"
                print(error)
                errors.append(error)
                continue

            globalert = False
            output_files = []
            glob_files = []
            for x in expected_outputs:
                if "*" in x:
                    current_glob = list(glob.glob(x))
                    if len(current_glob) == 0:
                        print("missing glob files", x, current_glob)
                        globalert = True
                        break
                    glob_files += current_glob
                else:
                    output_files.append(x)
            if globalert:
                error = f"missing glob files, continuing {expected_outputs}"
                print(error)
                errors.append(error)
                continue
            outputs_exist = [os.path.isfile(fn) for fn in output_files]
            if not all(outputs_exist):
                error = f"not all exist {list(zip(output_files, outputs_exist))}"
                print(error)
                errors.append(error)
                continue

            for output_file in output_files + glob_files:
                return_files[output_file] = zlib.compress(open(output_file, "rb").read())
            if save_logs:
                return_files["out.log"] = zlib.compress(open("out.log", "rb").read())
            success = True
            break
        if not success:
            if save_broken_dir:
                subprocess.call(["cp", "-r", os.getcwd(), os.path.expanduser(save_broken_dir)])
            error_content = (
                f"COMMAND: {' '.join(runcommand)} failed in the directory"
                f"\n{os.getcwd()} {list(input_files.keys())} {expected_outputs} tries: {current_try} : {max_tries}"
                f"\nerrors: {errors}\n"
                f"{os.path.join(os.path.expanduser(save_broken_dir), os.path.basename(os.getcwd()))}"
            )
            print(error_content)
            raise RuntimeError(error_content)
        # subprocess.call(["cp", "-r",  os.getcwd(), os.path.expanduser(save_broken_dir)])
        # print("Saved to", os.path.expanduser(save_broken_dir))
    return return_files

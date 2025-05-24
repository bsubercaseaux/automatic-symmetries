import time
from subprocess import TimeoutExpired, check_output, CalledProcessError, STDOUT
from plotter import plot_solution

def system_call(command, timeout=None):
    """
    params:
        command: list of strings, ex. ["ls", "-l"]
        timeout: number of seconds to wait for the command to complete.
    returns: output, return_code
    """
    try:
        cmd_output = check_output(command, stderr=STDOUT, timeout=timeout).decode()
        return_code = 0
    except CalledProcessError as e:
        cmd_output = e.output.decode()
        return_code = e.returncode
    except TimeoutExpired:
        cmd_output = f"Command timed out after {timeout} seconds"
        return_code = (
            -1
        )  # You can use any number that is not a valid return code for a success or normal failure
    return cmd_output, return_code


def timed_run_shell(commands, timeout=None):
    """
    params:
        command: list of strings, ex. ["ls", "-l"]
        timeout: number of seconds to wait for the command to complete.
    returns: output, return_code, elapsed_time in mseconds
    """
    start_time = time.perf_counter_ns()
    output, return_code = system_call(commands, timeout=timeout)
    elapsed_time = time.perf_counter_ns() - start_time
    return output, return_code, elapsed_time / 1e9


def test_realizability(input_file, timeout=None, realizer_path="localizer"):
    output, return_code, elapsed_time = timed_run_shell([realizer_path, input_file, "-i 10", "-r 30000"], timeout=timeout)
    if return_code == -1:
        return False
    else:
        return True
        
def test_and_output(input_file, timeout=None, points_output_file="out.txt", realizer_path="localizer", fix_file=None, sym_file=None):
    print("input_file = ", input_file, "points_output_file = ", points_output_file)
    cmds = [realizer_path, input_file, "-t", "4", "-i", "10", "-r", "30000", "-o", points_output_file]
    if fix_file is not None:
        cmds.extend(["-f", fix_file])
    if sym_file is not None:
        cmds.extend(["-c", sym_file])
    output, return_code, elapsed_time = timed_run_shell(cmds, timeout=timeout)
    if return_code == -1:
        return False
    else:
        print(f"output = {output}")
        plot_solution(points_output_file)
        return True

if __name__ == "__main__":
    import argparse
    argparser = argparse.ArgumentParser()
    argparser.add_argument("-f", "--file", type=str, required=True, help="Path to the input orientation file")
    argparser.add_argument("-t", "--timeout", type=int, default=5, help="Timeout in seconds")
    argparser.add_argument("-o", "--output", type=str, default="out.txt", help="Output file for points")
    args = argparser.parse_args()
    print(f"Realizable: {test_and_output(args.file, timeout=args.timeout, points_output_file=args.output)}")

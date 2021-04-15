import os
import numpy as np
import sys
import fileinput
import random
import subprocess
import time
import math
import ROOT
import psutil

path_to_build = "~/work/SHiP/build-SHiP_SBT_LScin"
path_to_source = "~/work/SHiP/SHiP_SBT_LScin"

def change_value(new_value, prefix, suffix, path_to_file):
    """ 
    Changes a line of text that contains a marker
    to a new one.

    Args:
        new_value (float): A value to be changed
        prefix (str): Text that marks the begining the desired line of code
        suffix (str): Text that won't change in this line
        path_to_file (str): Path to existing file that needs to be changed
    """

    new_text = prefix + str(new_value) + suffix

    x = fileinput.input(files=path_to_file, inplace=1)

    for line in x:
        if prefix in line:
            line = new_text
        sys.stdout.write(line)

    pass


def wait_while_running(proc_name, start):
    """
    Checks if there are screen sessions running with the name
    file_name. Displays for how long the script already runs

    Args:
        proc_name (str): names of process to be checked
        start (time.time()): time-stamp when program wath executed
    """
    still_running = proc_name in (p.name() for p in psutil.process_iter())

    while still_running:
        minutes = math.floor((time.time() - start)/60)
        seconds = int((time.time() - start)%60)
        print('Running for: {} m {} s'.format(minutes, seconds)) 
        time.sleep(1)
        still_running = proc_name in (p.name() for p in psutil.process_iter())

def compile(nThreads, path_to_build, path_to_source):
    os.system(f"(cd {path_to_build};cmake {path_to_source};make -j{nThreads})")

def run(nThreads, fileName):
    os.system(f"{path_to_build}/OpNovice -m {path_to_build}/run1.mac -t {nThreads}")


# compile(6, path_to_build, path_to_source)

change_value(150, "  SteelZ = ", "*mm;\n", "../src/OpNoviceDetectorConstruction.cc")

# os.system(f"cat {path_to_source}/src/OpNoviceDetectorConstruction.cc")

# length_box = 155
# length_WOM = 50

# start = time.time()

# while (length_box < 340):
#     for i in range(3):
#         print(f"Box length: {length_box} mm;", 
#               f"WOM length: {length_WOM} mm;")

#         change_value(length_box, "/WOMdir/boxWidthCmd ", " mm\n", path_to_build + "run1.mac")
#         change_value(f"test_{length_box}_{length_WOM}.root", "/analysis/setFileName ", "\n", path_to_build + "run1.mac")

#         run(6, f"test_{length_box}_{length_WOM}")
#         wait_while_running("OpNovice", start)

#         length_WOM += 50

#     length_WOM -= 120
#     length_box += 30
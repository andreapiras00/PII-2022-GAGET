import os
import sys
from os.path import exists, getsize
from scripts.graph_parser import parse_gfa2


def is_empty_file(fpath):
    return not fpath or not exists(fpath) or getsize(fpath) < 10

def get_input_name(fpath):
    file_name = fpath.rsplit('/', 1)[-1]
    name = file_name.rsplit('.', 1)[0]
    return name

def get_path_to_program(program, dirpath=None, min_version=None):
    """
    returns the path to an executable or None if it can't be found
    """
    def is_exe(fpath):
        if os.path.isfile(fpath) and os.access(fpath, os.X_OK):
            return True

    if dirpath:
        exe_file = os.path.join(dirpath, program)
        if is_exe(exe_file):
            return exe_file

    for path in os.environ["PATH"].split(os.pathsep):
        exe_file = os.path.join(path, program)
        if is_exe(exe_file):
            return exe_file
    return None

def parse_assembler_output(input_fpath, output_dirpath, min_edge_len):
	if not is_empty_file(input_fpath):
		if not input_fpath:
			sys.exit("ERROR! Failed parsing " + input_fpath + " file.")
		if input_fpath.endswith("gfa") or input_fpath.endswith("gfa2"):
			g = parse_gfa2(input_fpath)
	return g

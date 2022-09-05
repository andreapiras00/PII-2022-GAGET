import subprocess
from scripts.utils import is_empty_file, get_path_to_program, get_input_name
from os.path import join, exists

def run_graph_alignment(input_graph, input_ref, output_dir, n_threads):
	if not is_empty_file(input_graph):
		if not is_empty_file(input_ref):

			graphaligner_file = "GraphAligner_" + get_input_name(input_ref) + "_" + get_input_name(input_graph) + ".json" 
			graphaligner_out_fpath = join(output_dir, graphaligner_file)
			graphaligner_exec_path = get_path_to_program("GraphAligner")
		
			if not graphaligner_exec_path:
				print("GraphAligner is not found!")
				return None 

			if exists(graphaligner_out_fpath):
				print('Json file already present, skipping alignment process\n')
				return graphaligner_out_fpath

			parameter_preset = "vg"

			cmdline = [graphaligner_exec_path, "-g", input_graph, "-f", input_ref, "-a", graphaligner_out_fpath,
							"-t", str(n_threads), "-x", parameter_preset]#, "--all-alignments"]

			return_code = subprocess.call(cmdline, stdout=open("/dev/null", "w"))#, stderr=open("/dev/null", "w")) 

			if return_code != 0 or is_empty_file(graphaligner_out_fpath):
				print("Warning! GraphAligner failed aligning the assembly graph to the reference")

			return graphaligner_out_fpath

		else:
			print("ERROR: input graph file (.gfa) is not valid!!!")
			return None
	else:
		print("ERROR: input reference file (FASTA) is not valid!!!")
		return None
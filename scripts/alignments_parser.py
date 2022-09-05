import json
import sys
import re
from scripts.utils import is_empty_file
from os.path import exists
from Bio import SeqIO
from collections import defaultdict
from gfapy.sequence import rc
from skbio.alignment import StripedSmithWaterman


def parse_json(alignments_fpath):
	print("Parsing " + alignments_fpath + "...")
	if exists(alignments_fpath) and not is_empty_file(alignments_fpath):
		alignments = [json.loads(line) for line in open(alignments_fpath,'r')]
		return alignments
	else:
		sys.exit("ERROR while opening alignment results file!")


def parse_reference(ref_fpath):
	print("Parsing " + ref_fpath + "...")
	if exists(ref_fpath) and not is_empty_file(ref_fpath):
		reference = {}
		fasta_seqs = SeqIO.parse(open(ref_fpath),'fasta')
		for seq in fasta_seqs:
			seq.id = str(seq.id).split()[0]
			reference[seq.id] = {'seq':str(seq.seq).upper(), 'len':len(str(seq.seq)),'pos_mapped':defaultdict(list)}
		print("Finish parsing\n")
		return reference
	else:
		sys.exit("ERROR while opening reference file!")


def get_value_alignment(alignment, key):
	try:
		return alignment[key]
	except:
		if key in ["score", "query_position"]:
			return 0
		else:
			return None


def get_value_path(path, key):
	try:
		return path[key]
	except:
		if key in ["is_circular"]:
			return False
		else:
			return None


def get_value_mapping(mapping, key):
	try:
		return mapping[key]
	except:
		if key in ["rank"]:
			return 0
		else:
			return None


def get_value_position(position, key):
	try:
		return position[key]
	except:
		if key in ["offset"]:
			return 0
		elif key in ["is_reverse"]:
			return False
		else:
			return None


def get_value_edit(edit, key):
	try:
		return edit[key]
	except:
		if key in ["sequence"]:
			return ""
		elif key in ["to_length", "from_length"]:
			return 0
		else:
			return None

def align_info2(input_data, graph, ref, n_align):
	i_al = input_data[0]
	al = input_data[1]

	print("Parsing alignment {}/{}".format(i_al, n_align))
	ref_name = get_value_alignment(al,"name").split()[0]
	info = ([], [], i_al, ref_name)
	ref_seq = ref[ref_name]["seq"]
	i_ref = int(get_value_alignment(al,"query_position"))
	path = get_value_alignment(al,"path")


	for mapping in get_value_path(path,"mapping"):
		node = get_value_position(mapping["position"],"node_id")
		offset = int(get_value_position(mapping["position"],"offset"))
		reverse = get_value_position(mapping["position"],"is_reverse")
		if reverse:
			node_seq = rc(graph.nodes[node]["seq"])
		else:
			node_seq = graph.nodes[node]["seq"]
		i_node = offset
		
		for edit in get_value_mapping(mapping,"edit"):
			score = 0
			edit_seq = get_value_edit(edit,"sequence").upper()
			if get_value_edit(edit,"from_length") == get_value_edit(edit,"to_length"): #and get_value_edit(edit,"from_length")!=0:
				if edit_seq == "":
					edit_seq = ref_seq[i_ref:i_ref+int(get_value_edit(edit,"to_length"))] 
				if edit_seq != ref_seq[i_ref:i_ref+int(get_value_edit(edit,"to_length"))]:
					sys.exit("ERROR: edit_sequence is not equal to ref_sequence:\ne_seq = " + edit_seq +"\nr_seq = " + ref_seq[i_ref:i_ref+int(get_value_edit(edit,"to_length"))])
				sm_ref = i_ref
				em_ref = i_ref
				#i_node = offset
				sm_node = i_node
				em_node = i_node
				mapped = False
				count_match = 0
				count_mismatch = 0
				CIGAR = ""	
					
				i_edit = 0

				for i_edit in range(0,int(get_value_edit(edit,"to_length"))):
						# print("node len = {}".format(graph.nodes[node]["size"]))
						# print(f"node index = {i_node}")
						# print(f"edit len = {len(edit_seq)}")
						# print(f"edit index = {i_edit}")

						if node_seq[i_node] != edit_seq[i_edit]:
							if mapped:
								mapped = False
								count_mismatch = 1
								CIGAR = CIGAR + str(count_match) + '='
								score += count_match
								count_match = 0
							else:
								count_mismatch += 1
						else: 
							if not mapped:
								mapped = True 
								count_match = 1
								if count_mismatch != 0:
									CIGAR = CIGAR + str(count_mismatch) + 'X' 
									score -= count_mismatch
									count_mismatch = 0
							else:
								count_match += 1
						i_ref += 1
						i_node += 1
				
				em_node = i_node - 1
				em_ref = i_ref - 1
				if mapped:
					CIGAR = CIGAR + str(count_match) + '='
					score += count_match
					count_match = 0
				else:
					CIGAR = CIGAR + str(count_mismatch) + 'X'
					score -= count_mismatch
					count_mismatch = 0
				info[1].append((sm_ref, em_ref, CIGAR, score, node))
				if reverse:
					info[0].append((node, (graph.nodes[node]["size"])-em_node-1, (graph.nodes[node]["size"]-sm_node-1), CIGAR, score, ref_name, sm_ref, em_ref, reverse))
				else:
					info[0].append((node, sm_node, em_node, CIGAR, score, ref_name, sm_ref, em_ref, reverse))

			elif get_value_edit(edit,"to_length") != 0:
				if edit_seq != ref_seq[i_ref:i_ref+int(get_value_edit(edit,"to_length"))]:
					sys.exit("ERROR: edit_sequence is not equal to ref_sequence:\ne_seq = " + edit_seq +"\nr_seq = " + ref_seq[i_ref:i_ref+int(get_value_edit(edit,"to_length"))])
				if get_value_edit(edit,"from_length")!=0:
					query = StripedSmithWaterman(node_seq[offset:])
					loc_al = query(edit_seq)
					overlap_operations = re.split('(\d+)', loc_al.cigar.strip())
					i_ref += loc_al.target_begin
					sm_ref = i_ref
					em_ref = i_ref
					i_node += loc_al.query_begin
					sm_node = i_node
					em_node = i_node
					mapped = False
					i_edit = loc_al.target_begin
					count_match = 0
					count_mismatch = 0
					CIGAR = ""

					for i in range(0, len(overlap_operations) - 1):
						if not overlap_operations[i]:
							continue

						if overlap_operations[i+1] == 'M':
							for i_edit in range(i_edit, i_edit+int(overlap_operations[i])):
								if node_seq[i_node] != edit_seq[i_edit]:
									if mapped:
										mapped = False
										count_mismatch = 1
										CIGAR = CIGAR + str(count_match) + '='
										score += count_match
										count_match = 0
									else:
										count_mismatch += 1
								else: 
									if not mapped:
										mapped = True 
										count_match = 1
										if count_mismatch != 0:
											CIGAR = CIGAR + str(count_mismatch) + 'X' 
											score -= count_mismatch
											count_mismatch = 0
									else:
										count_match += 1
								i_ref += 1
								i_node += 1

							if mapped:
								CIGAR = CIGAR + str(count_match) + '='
								score += count_match
								count_match = 0
							else:
								CIGAR = CIGAR + str(count_mismatch) + 'X'
								score -= count_mismatch
								count_mismatch = 0

						elif overlap_operations[i+1] == 'I':
							mapped = False
							if overlap_operations[i-1] == 'M':
								i_edit += 1
							i_node += int(overlap_operations[i])
							CIGAR = CIGAR + overlap_operations[i] +'I' 
							score -= int(overlap_operations[i])
							
						elif overlap_operations[i+1] == 'D':
							mapped = False
							if overlap_operations[i-1] == 'M':
								i_edit += 1
							i_ref += int(overlap_operations[i])
							i_edit += int(overlap_operations[i])
							CIGAR = CIGAR + overlap_operations[i] +'D' 
							score -= int(overlap_operations[i])	
					
					i_ref += int(get_value_edit(edit,"to_length")) - i_edit
					if overlap_operations[-1] == 'M':
						i_ref -= 1
						i_node -= 1
					em_node = i_node #-1
					em_ref = i_ref-1

					if i_edit != loc_al.target_end_optimal or i_node != offset+loc_al.query_end:
						print("Qualcosa non va...")
						print(str(i_edit) + "!=" + str(loc_al.target_end_optimal))
						print(str(i_node) + "!=" + str(offset+loc_al.query_end))
				else:
					sm_ref = i_ref
					em_ref = i_ref
					sm_node = i_node
					em_node = i_node
					i_ref +=  get_value_edit(edit,"to_length")
					em_ref = i_ref-1
					CIGAR = str(get_value_edit(edit,"to_length")) + 'D'
					score -= get_value_edit(edit,"to_length")

				
				info[1].append((sm_ref, em_ref, CIGAR, score, node))
				if reverse:
					info[0].append((node, (graph.nodes[node]["size"])-em_node-1, (graph.nodes[node]["size"]-sm_node-1), CIGAR, score, ref_name, sm_ref, em_ref, reverse))
				else:
					info[0].append((node, sm_node, em_node, CIGAR, score, ref_name, sm_ref, em_ref, reverse))

			elif get_value_edit(edit,"to_length") == 0:
				sm_ref = i_ref
				em_ref = i_ref
				sm_node = i_node
				em_node = i_node
				i_node += get_value_edit(edit,"from_length")
				em_node = i_node-1
				CIGAR = str(get_value_edit(edit,"from_length")) + 'I'
				score -= get_value_edit(edit,"from_length")

				info[1].append((sm_ref, em_ref, CIGAR, score, node))
				if reverse:
					info[0].append((node, (graph.nodes[node]["size"])-em_node-1, (graph.nodes[node]["size"]-sm_node-1), CIGAR, score, ref_name, sm_ref, em_ref, reverse))
				else:
					info[0].append((node, sm_node, em_node, CIGAR, score, ref_name, sm_ref, em_ref, reverse))

			else:
				print("WARNING: edit's info not parsed!")

	return info

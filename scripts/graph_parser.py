import re
from collections import defaultdict

import gfapy
import networkx as nx


def parse_gfa2(gfa_fpath): 

    g = nx.DiGraph()

    #print("Parsing " + gfa_fpath + "...")
    gfa = gfapy.Gfa.from_file(gfa_fpath, vlevel = 0)
    
    for segment in gfa.segments:
        g.add_node(segment.name, seq=segment.sequence.upper(), size=len(segment.sequence), unique_size=len(segment.sequence), pos_mapped=defaultdict(list))
    for edge in gfa.edges:
        if not "GFAPY_virtual_line" in str(edge).split()[-1]:
            g.add_edge(edge.from_segment.name, edge.to_segment.name, from_ori=edge.from_orient, to_ori=edge.to_orient, overlap=str(edge.overlap))

    for node in list(g.nodes()):
        min_overlap = None
        for in_edge in g.in_edges(node):
            overlap = 0
            overlap_operations = re.split('(\d+)', g.edges[in_edge]["overlap"].strip())
            for i in range(0, len(overlap_operations) - 1):
                if not overlap_operations[i]:
                    continue
                if overlap_operations[i+1] == 'M' or overlap_operations[i+1] == 'I': #not sure, maybe should check for D instead
                    overlap += int(overlap_operations[i])
            if not min_overlap:
                min_overlap = overlap
            elif overlap < min_overlap:
                min_overlap = overlap
        if not min_overlap:
            min_overlap = 0
        g.nodes()[node]["unique_size"] = g.nodes()[node]["size"] - min_overlap

    return g

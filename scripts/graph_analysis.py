from collections import defaultdict
import copy

class Alignment(object):
    ''' Alignment class is used to convert the alignments from the aligner tool, in order to make access and
    modifications to them easier'''
    __slots__ = ("path", "path_id", "s_ref", "e_ref", "s_node", "node1", "e_node", "node2", "score", "ref")

    def __init__(self, path, path_id, s_ref, e_ref, s_node, node1, e_node, node2, score, ref):
        self.path, self.path_id, self.s_ref, self.e_ref, self.s_node, self.node1, self.e_node, self.node2, self.score, self.ref = path, path_id, s_ref, e_ref, s_node, node1, e_node, node2, score, ref
    
    def __str__(self):
        return ' '.join(str(x) for x in ['ref:', self.ref, 'start:', self.s_ref, 'end:', self.e_ref, '|', 'start node:', self.node1, 'position:', self.s_node, '|', 'end node:', self.node2, 'position:', self.e_node])

    def clone(self):
        # clone method creates a deep copy of the current alignment.
        return Alignment(copy.deepcopy(self.path), self.path_id, self.s_ref, self.e_ref, self.s_node, self.node1, self.e_node, self.node2, self.score, self.ref)
    
    def start(self):
        # start method returns the alignment's leftmost position on the reference.
        return min(self.s_ref, self.e_ref)

    def end(self):
        # start method returns the alignment's rightmost position on the reference.
        return max(self.s_ref, self.e_ref)

    def remove_last_node(self, g):
        # remove_last_node method removes the last node of the alignment's path and updates all required info, returning
        # the score of the removed node.
        removed = self.path.pop(-1)
        self.e_ref = self.path[-1][1]
        self.score -= removed[3]
        self.node2 = self.path[-1][4]
        self.e_node = g.nodes[self.node2]['pos_mapped'][self.path_id][-1][1]
        return removed[3]  

    def remove_first_node(self, g):
        # remove_first_node method removes the first node of the alignment's path and updates all required info, returning
        # the score of the removed node.
        removed = self.path.pop(0)
        self.s_ref = self.path[0][0]
        self.score -= removed[3]
        self.node1 = self.path[0][4]
        self.s_node = g.nodes[self.node1]['pos_mapped'][self.path_id][0][0]
        return removed[3] 

    def distance_on_ref(self, other):
        # distance_on_ref method returns the distance in bp between the alignments and another one.
        return other.start() - self.end()    


def score(set, curr_align, g):
    # score method removes the overlapping between the last alignments of the set and the current alignment,
    # adds the align to the set, where one of the two has been reduced in case of an overlap and calculates the set's
    # new score
    prev_align = set[0][-1]
    prev_score = set[1]
    al_score = curr_align.score
    new_score = prev_score + al_score
    # new align and set don't overlap, return set with new align
    if curr_align.start() >= prev_align.end():
        set[0].append(curr_align)
        # remove gap distance from score
        new_score -= prev_align.distance_on_ref(curr_align)
        return (set[0], new_score)  
    else: 
        # reduce last align of set until no overlap with new align
        prev_align_copy = prev_align.clone()
        set[0][-1] = prev_align_copy
        reduced_set_score = reduced_align_score = new_score
        while len(set[0][-1].path) > 1 and curr_align.start() <= prev_align_copy.end():
            reduced_set_score -= set[0][-1].remove_last_node(g)
        # remove gap distance from score    
        reduced_set_score -= set[0][-1].distance_on_ref(curr_align)
        # if last align of set ends after last node of new align cannot reduce new align, return reduced set
        if prev_align.end() > curr_align.path[-1][0]:
            set[0].append(curr_align)
            return (set[0], reduced_set_score)
        else:    
            # reduce new align until no overlap with set 
            set[0][-1] = prev_align    
            curr_align_copy = curr_align.clone()
            while len(curr_align_copy.path) > 0 and curr_align_copy.start() < prev_align.end():
                reduced_align_score -= curr_align_copy.remove_first_node(g)
            # remove gap distance from score  
            reduced_set_score -= prev_align.distance_on_ref(curr_align_copy)
            # return the set with the best score (reduced align or reduced set)
            if reduced_set_score > reduced_align_score:
                set[0][-1] = prev_align_copy
                set[0].append(curr_align)
                new_score = reduced_set_score
            else:  
                set[0].append(curr_align_copy)
                new_score = reduced_align_score  

            return (set[0], new_score)   

def score_set(alignment, sets, g):
    # score_set method calculates the score for each candidate set with the current alignment and returns the
    # one with the biggest score
    candidates = []
    for al_set in sets:
        if not al_set[0]:
            candidates.append(([alignment], alignment.score))    
        # HEURISTIC: skip sets on which the last alignment's first node ends after the curr align's start position.
        elif alignment.start() <= al_set[0][-1].path[0][1]:
            continue
        else:
            candidates.append(score(al_set, alignment, g))           
    # return argmaxB∈BestSets Score(B ∪ ai )   
    max_alignment = candidates[0]
    for candidate in candidates:
        if candidate[1] > max_alignment[1]:
            max_alignment = candidate     
    return max_alignment

def best_alignments_selection(g, reference):  
    # best_alignments_selection method returns the best set of non-overlapping alignments that map the reference. 
    # BestSets ← {(EmptyAlignment, 0)}
    BEST = defaultdict(lambda: defaultdict())
    for ref in reference:
        pos_mapped = reference[ref]['pos_mapped']
        # initialized with empty alignment and score 0
        best_ref_sets = [([], 0)]
        alignments = [] 
        for path in pos_mapped:
            # score: total bases mapped
            score = sum([pos[3] for pos in pos_mapped[path]])
            start = pos_mapped[path][0]
            end = pos_mapped[path][-1]
            start_ref = int(start[0])
            end_ref = int(end[1])
            start_node_id = start[4]
            end_node_id = end[4]
            start_node = int(g.nodes[start_node_id]['pos_mapped'][path][0][0])
            end_node = int(g.nodes[end_node_id]['pos_mapped'][path][-1][1])
            align = Alignment(pos_mapped[path], path, start_ref, end_ref, start_node, start_node_id, end_node, end_node_id, score, ref)
            alignments.append(align)
        # order alignments based on last position mapped    
        sorted_alignments = sorted(alignments, key=lambda x : x.end())  

        # finding best alignments based on their position on the reference
        for al in sorted_alignments:
            # besti ← argmaxB∈BestSets Score(B ∪ ai )
            best = score_set(al, best_ref_sets, g)
            # BestSets ← BestSets ∪ (besti, Score(besti ))
            best_ref_sets.append(best)
        # return argmaxB∈BestSets Score(B)
        max_set = best_ref_sets[0]
        for walk in best_ref_sets:
            if walk[1] > max_set[1]:
                max_set = walk                        
        BEST[ref]["reference"] = max_set  
    return BEST
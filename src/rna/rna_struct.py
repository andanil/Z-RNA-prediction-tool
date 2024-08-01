import RNA

from ..utils import reverse_compl
from .utils import dotbracket_to_bp
from .rna_graph import RNAGraph
from .bulge import Bulge


class RNAStructure():
    def __init__(self, seq):
        self._seq = seq.replace("T", "U")
        self._reversed_seq = reverse_compl(seq).replace("T", "U")
        self._forward_struct = self._calculate_struct(self._seq)
        self._reversed_struct = self._calculate_struct(self._reversed_seq)
        self.forward = True
        self._graph = None

    def __len__(self):
        return len(self._seq)

    @property
    def seq(self):
        return self._seq if self.forward else self._reversed_seq

    def get_struct(self, forward=True):
        return self._forward_struct if forward else self._reversed_struct

    def reset_graph_type(self, forward=True):
        self.forward = forward
        self._graph = None

    def get_stems(self):
        return self._get_graph_instance(self.forward).stem_iterator()

    def get_stems_union(self):
        stems = [s for s in self._get_graph_instance(self.forward).stem_iterator()]
        stems_union = []
        for i in range(1, len(stems)):
            left = self.get_stem_coord(stems[i-1])
            right = self.get_stem_coord(stems[i])

            if left[1] == right[0] - 1 or left[2] == right[3] + 1:
                n_unpaired = max(right[0] - left[1] - 1, left[2] - right[3] - 1)
                if right[0] - left[1] - 1 > left[2] - right[3] - 1:
                    unpaired = [left[1] + 1, right[0] - 1]
                    opposite = [right[3], left[2]]
                else:
                    unpaired = [right[3] + 1, left[2] - 1]
                    opposite = [left[1], right[0]]
                if n_unpaired <= 10:
                    stems_union += [[left[0], right[1], right[2], left[3], [Bulge(unpaired, opposite)]]]
        i = 1
        changed = True
        while changed:
            changed = False
            i = 1
            while i < len(stems_union):
                if (stems_union[i-1][1] > stems_union[i][0] and stems_union[i-1][1] < stems_union[i][1]) or (stems_union[i-1][2] < stems_union[i][3] and stems_union[i-1][2] > stems_union[i][2]):
                    stems_union[i-1] = [stems_union[i-1][0], stems_union[i][1], stems_union[i][2], stems_union[i-1][3], [*stems_union[i-1][4], *stems_union[i][4]]]
                    del stems_union[i]
                    changed = True
                i += 1
        return stems_union

    def get_stem_bp(self, s):
        return self._get_graph_instance(self.forward).stem_bp_iterator(s)

    def get_stems_coord(self):
        return [self.get_stem_coord(s) for s in self._get_graph_instance(self.forward).stem_iterator()]

    def get_stem_coord(self, s):
        return self._get_graph_instance(self.forward).stem_range(s)

    def get_pairing_partner(self, nt):
        return self._get_graph_instance(self.forward).pairing_partner(nt)

    def get_elem(self, inx):
        return self._get_graph_instance(self.forward).get_node_from_residue_num(inx)

    def get_bp(self, forward=True):
        return dotbracket_to_bp(self.get_struct(forward))

    def _calculate_struct(self, seq):
        (ss, _) = RNA.fold(seq)
        return ss

    def _get_graph_instance(self, forward):
        if not self._graph:
            self._graph = RNAGraph.from_dotbracket(self.get_struct(forward))
        return self._graph

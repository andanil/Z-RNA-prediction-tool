from collections import defaultdict
import itertools as it
import sys

from .utils import pairtable_to_tuples, dotbracket_to_pairtable, any_difference_of_one


class RNAGraph():
    def __init__(self, tuples, seq_length):
        self.defines = {}
        self.edges = defaultdict(set)
        self.weights = {}
        self.seq_length = seq_length
        self._name_counter = 0
        self._from_tuples(tuples)

    @classmethod
    def from_dotbracket(cls, dotbracket):
        """
        Create a graph from a dotbracket string.
        """
        tuples = pairtable_to_tuples(dotbracket_to_pairtable(dotbracket))
        rna_g = cls(tuples, len(dotbracket))
        return _cleaned_graph(rna_g)

    def connections(self, bulge):
        def sort_key(x):
            if self.defines[x] and self.defines[x][0] == 0:
                return -1
            return self.define_a(x)[0]
        connections = list(self.edges[bulge])
        connections.sort(key=sort_key)
        return connections

    def stem_length(self, key):
        """
        Get the length of a particular element. If it's a stem, it's equal to
        the number of paired bases. If it's an interior loop, it's equal to the
        number of unpaired bases on the strand with less unpaired bases. If
        it's a multiloop, then it's the number of unpaired bases.
        """
        d = self.defines[key]
        if key[0] == 's' or key[0] == 'y':
            return (d[1] - d[0]) + 1
        elif key[0] == 'f':
            return self._get_bulge_dimensions(key)[0]
        elif key[0] == 't':
            return self._get_bulge_dimensions(key)[1]
        elif key[0] == 'h':
            return self._get_bulge_dimensions(key)[0]
        else:
            return min(self._get_bulge_dimensions(key))

    def stem_bp_iterator(self, stem):
        """
        Iterate over all the base pairs in the stem.
        """
        assert stem[0] == "s"
        d = self.defines[stem]
        for i in range(self.stem_length(stem)):
            yield (d[0] + i, d[3] - i)

    def stem_range(self, stem):
        assert stem[0] == "s"
        return self.defines[stem]

    def stem_iterator(self):
        """
        Iterator over all of the stems in the structure.
        """
        for d in self.defines.keys():
            assert d[0] in "ftsmih", "stem_iterator should only be called after relabelling of nodes during GraphConstruction"
            if d[0] == 's':
                yield d

    def get_node_from_residue_num(self, base_num):
        """
        Iterate over the defines and see which one encompasses this base.
        """
        for key in self.defines:
            for drange in self._define_range_iterator(key):
                if drange[0] <= base_num <= drange[1]:
                    return key
        raise LookupError("Base number {} not found in the defines {}.".format(base_num, self.defines))

    def pairing_partner(self, nucleotide_number):
        """
        Return the base pairing partner of the nucleotide at position
        nucleotide_number. If this nucleotide is unpaired, return None.
        :param nucleotide_number: The position of the query nucleotide in the
                                  sequence or a RESID instance.
        :return: The number of the nucleotide base paired with the one at
                 position nucleotide_number.
        """
        for d in self.stem_iterator():
            for (r1, r2) in self.stem_bp_iterator(d):
                if r1 == nucleotide_number:
                    return r2
                elif r2 == nucleotide_number:
                    return r1
        return None

    def length_one_stem_basepairs(self):
        """
        Return a list of basepairs that correspond to length-1 stems.
        """
        stems_to_dissolve = [s for s in self.stem_iterator() if self.stem_length(s) == 1]
        bps_to_dissolve = []
        for s in stems_to_dissolve:
            bps_to_dissolve.extend(self.stem_bp_iterator(s))
        return bps_to_dissolve

    def to_pair_tuples(self, remove_basepairs=None):
        """
        Create a list of tuples corresponding to all of the base pairs in the
        structure. Unpaired bases will be shown as being paired with a
        nucleotide numbered 0.
        i.e. [(1,5),(2,4),(3,0),(4,2),(5,1)]
        :param remove_basepairs: A list of 2-tuples containing
                                 basepairs that should be removed
        """
        table = []

        for d in self.defines:
            for b in self._define_residue_num_iterator(d):
                p = self.pairing_partner(b)
                if p is None:
                    p = -1
                table.append((b, p))

        if remove_basepairs:
            nt = []
            for p in table:
                to_add = p
                for s in remove_basepairs:
                    if sorted(p) == sorted(s):
                        to_add = (p[0], -1)
                        break
                nt += [to_add]
            table = nt

        return table

    def define_a(self, elem):
        if elem[0] == "i":
            conns = self.connections(elem)
            s1 = self.defines[conns[0]]
            s2 = self.defines[conns[1]]
            return [s1[1], s2[0], s2[3], s1[2]]
        else:
            if self.defines[elem] == []:
                return self._define_a_zerolength(elem)
            return self._define_a_nonzero(elem)

    def _define_residue_num_iterator(self, node):
        """
        Iterate over the residue numbers that belong to this node.
        :param node: The name of the node
        """
        visited = set()
        for r in self._define_range_iterator(node):
            for i in range(r[0], r[1] + 1):
                if i not in visited:
                    visited.add(i)
                    yield i

    def _define_range_iterator(self, node):
        """
        Return the ranges of the nucleotides in the define.
        In other words, if a define contains the following: [1,2,7,8]
        The ranges will be [1,2] and [7,8].
        :return: A list of two-element lists
        """
        define = self.defines[node]
        if define:
            yield [define[0], define[1]]
            if len(define) > 2:
                yield [define[2], define[3]]

    def _merge_vertices(self, vertices):
        """
        This is done when two of the outgoing strands of a stem
        go to different bulges
        It is assumed that the two ends are on the same sides because
        at least one vertex has a weight of 2, implying that it accounts
        for all of the edges going out of one side of the stem
        :param vertices: A list of vertex names to combine into one.
        """
        new_vertex = self._get_vertex()
        self.weights[new_vertex] = 0

        connections = set()

        for v in vertices:
            for item in self.edges[v]:
                connections.add(item)

            if v[0] == 's':
                self.defines[new_vertex] = self.defines.get(new_vertex, []) + [self.defines[v][0],
                                                            self.defines[v][2]] + [self.defines[v][1], self.defines[v][3]]
            else:
                self.defines[new_vertex] = self.defines.get(new_vertex, []) + self.defines[v]

            self.weights[new_vertex] += 1
            self._remove_vertex(v)
            self._reduce_defines()

        for connection in connections:
            self.edges[new_vertex].add(connection)
            self.edges[connection].add(new_vertex)

        return new_vertex

    def _reduce_defines(self):
        """
        Make defines like this:
        define x0 2 124 124 3 4 125 127 5 5
        Into this:
        define x0 2 3 5 124 127
        That is, consolidate contiguous bulge region defines.
        """
        for key in self.defines.keys():
            if key[0] != 's':
                assert (len(self.defines[key]) % 2 == 0)
                new_j = 0

                while new_j < len(self.defines[key]):

                    j = new_j
                    new_j += j + 2

                    (f1, t1) = (int(self.defines[key][j]), int(self.defines[key][j + 1]))

                    # remove bulges of length 0
                    if f1 == -1 and t1 == -2:
                        del self.defines[key][j]
                        del self.defines[key][j]

                        new_j = 0
                        continue

                    # merge contiguous bulge regions
                    for k in range(j + 2, len(self.defines[key]), 2):
                        if key[0] == 'y':
                            # we can have stems with defines like: [1,2,3,4]
                            # which would imply a non-existant loop at its end
                            continue

                        (f2, t2) = (int(self.defines[key][k]), int(self.defines[key][k + 1]))

                        if t2 + 1 != f1 and t1 + 1 != f2:
                            continue

                        if t2 + 1 == f1:
                            self.defines[key][j] = str(f2)
                            self.defines[key][j + 1] = str(t1)
                        elif t1 + 1 == f2:
                            self.defines[key][j] = str(f1)
                            self.defines[key][j + 1] = str(t2)

                        del self.defines[key][k]
                        del self.defines[key][k]

                        new_j = 0

                        break

    def _remove_vertex(self, v):
        """
        Delete a node after merging it with another
        :param v: The name of the node
        """
        # delete all edges to this node
        for key in self.edges[v]:
            self.edges[key].remove(v)

        for edge in self.edges:
            if v in self.edges[edge]:
                self.edges[edge].remove(v)

        # delete all edges from this node
        del self.edges[v]
        del self.defines[v]

    def _compare_stems(self, b):
        return (self.defines[b][0], 0)

    def _compare_bulges(self, b):
        return self.define_a(b)

    def _compare_hairpins(self, b):
        connections = self.connections(b)
        return (self.defines[connections[0]][1], sys.maxsize)

    def _get_vertex(self, name=None):
        """
        Return a new unique vertex name starting with `x`.
        At this stage stems and bulges are not distinguished
        """
        if name is None:
            name = "x{}".format(self._name_counter)
            self._name_counter += 1
        return name

    def _get_bulge_dimensions(self, bulge):
        """
        Return the dimensions of the bulge.
        If it is single stranded it will be (x, -1) for h,t,f or (x, 1000) for m.
        Otherwise it will be (x, y).
        :param bulge: The name of the bulge.
        :return: A pair containing its dimensions
        """
        if bulge[0] == "s":
            raise ValueError("Stems are not allowed in get_bulge_dimensions")
        bd = self.defines[bulge]
        c = self.connections(bulge)

        if bulge[0] == 'i':
            s1 = self.defines[c[0]]
            s2 = self.defines[c[1]]
            return get_define_len([s1[1] + 1, s2[0] - 1]), get_define_len([s2[3] + 1, s1[2] - 1])
        else:
            if bd:
                dim0 = get_define_len([bd[0], bd[1]])
            else:
                dim0 = 0
            if bulge[0] == "m":
                dim1 = 1000
            else:
                dim1 = -1
            return dim0, dim1

    def _from_tuples(self, tuples):
        """
        Create a graph from a list of pair tuples.
        Unpaired nucleotides have a pairing partner of -1
        """
        stems = []
        bulges = []

        tuples.sort()
        tuples = iter(tuples)
        (prev_from, prev_to) = next(tuples)

        start_from = prev_from
        start_to = prev_to
        last_paired = prev_from

        for from_bp, to_bp in tuples:
            if abs(to_bp - prev_to) == 1 and prev_to != -1:  # adjacent basepairs on 3' strand
                # stem
                if (((prev_to - prev_from > 0 and to_bp - from_bp > 0) or
                     (prev_to - prev_from < 0 and to_bp - from_bp < 0)) and
                        (to_bp - prev_to) == -(from_bp - prev_from)):
                    (prev_from, prev_to) = (from_bp, to_bp)
                    last_paired = from_bp
                    continue

            if to_bp == -1 and prev_to == -1:
                # bulge
                (prev_from, prev_to) = (from_bp, to_bp)
                continue
            else:
                if prev_to != -1:
                    new_stem = tuple(sorted([tuple(sorted([start_from, start_to])),
                                             tuple(sorted([prev_from, prev_to]))]))
                    if new_stem not in stems:
                        stems += [new_stem]

                    last_paired = from_bp
                    start_from = from_bp
                    start_to = to_bp
                else:
                    new_bulge = ((last_paired, prev_from))
                    bulges += [new_bulge]

                    start_from = from_bp
                    start_to = to_bp

            prev_from = from_bp
            prev_to = to_bp

        # Take care of the last element
        if prev_to != -1:
            new_stem = tuple(sorted([tuple(sorted([start_from, start_to])),
                                     tuple(sorted([prev_from, prev_to]))]))
            if new_stem not in stems:
                stems += [new_stem]
        if prev_to == -1:
            new_bulge = ((last_paired, prev_from))
            bulges += [new_bulge]

        self._from_stems_and_bulges(stems, bulges)

    def _from_stems_and_bulges(self, stems, bulges):
        for i in range(len(stems)):
            ss1 = stems[i][0][0]
            ss2 = stems[i][0][1]
            se1 = stems[i][1][0]
            se2 = stems[i][1][1]
            self.defines['y%d' % (i)] = [min(ss1, se1), max(ss1, se1),
                                         min(ss2, se2), max(ss2, se2)]
            self.weights['y%d' % (i)] = 1

        for i, bulge in enumerate(bulges):
            self.defines['b%d' % (i)] = sorted([bulge[0], bulge[1]])
            self.weights['b%d' % (i)] = 1

        self._create_bulge_graph(stems, bulges)
        self._create_stem_graph(stems, len(bulges))
        self._collapse()
        self._sort_defines()
        self._relabel_nodes()
        self._remove_degenerate_nodes()

    def _create_bulge_graph(self, stems, bulges):
        """
        Find out which stems connect to which bulges
        Stems and bulges which share a nucleotide are considered connected.
        :param stems: A list of tuples of tuples of the form [((s1, e1), (s2, e2))]
                      where s1 and e1 are the nucleotides at one end of the stem
                      and s2 and e2 are the nucleotides at the other.
        :param bulges: A list of tuples of the form [(s, e)] where s and e are the
                       numbers of the nucleotides at the start and end of the bulge.
        """
        for i in range(len(stems)):
            stem = stems[i]
            for j, bulge in enumerate(bulges):
                if any_difference_of_one(stem, bulge):
                    self.edges['y{}'.format(i)].add('b{}'.format(j))
                    self.edges['b{}'.format(j)].add('y{}'.format(i))

    def _create_stem_graph(self, stems, bulge_counter):
        """
        Determine which stems are connected to each other. A stem can be connected to
        another stem when there is an interior loop with an unpaired nucleotide on
        one side. In this case, a bulge will be created on the other side, but it
        will only consist of the two paired bases around where the unpaired base
        would be if it existed.
        The defines for these bulges will be printed as well as the connection
        strings for the stems they are connected to.
        :param stems: A list of tuples of tuples of the form [((s1, e1), (s2, e2))]
                      where s1 and e1 are the nucleotides at one end of the stem
                      and s2 and e2 are the nucleotides at the other.
        :param bulge_counter: The number of bulges that have been encountered so far.
        :returns: None
        """
        for i, j in it.combinations(range(len(stems)), 2):
            for k1, k2, l1, l2 in it.product(range(2), repeat=4):
                s1 = stems[i][k1][l1]
                s2 = stems[j][k2][l2]
                if k1 == 1 and stems[i][0][l1] == stems[i][1][l1]:
                    continue
                if k2 == 1 and stems[j][0][l2] == stems[j][1][l2]:
                    continue
                if abs(s1 - s2) == 1:
                    bn = 'b{}'.format(bulge_counter)
                    self.defines[bn] = []
                    self.weights[bn] = 1

                    self.edges['y{}'.format(i)].add(bn)
                    self.edges[bn].add('y{}'.format(i))

                    self.edges['y{}'.format(j)].add(bn)
                    self.edges[bn].add('y{}'.format(j))

                    bulge_counter += 1

        for d in list(self.defines.keys()):
            if d[0] != 'y':
                continue

            (s1, e1, s2, e2) = self.defines[d]
            if abs(s2 - e1) == 1:
                bn = 'b{}'.format(bulge_counter)

                self.defines[bn] = []
                self.weights[bn] = 1

                self.edges[bn].add(d)
                self.edges[d].add(bn)

                bulge_counter += 1
        return

    def _sort_defines(self):
        """
        Sort the defines of interior loops and stems so that the 5' region
        is always first.
        """
        for k in self.defines.keys():
            d = self.defines[k]

            if len(d) == 4:
                if d[0] > d[2]:
                    new_d = [d[2], d[3], d[0], d[1]]
                    self.defines[k] = new_d
    
    def _collapse(self):
        """
        If any vertices form a loop, then they are either a bulge region or
        a fork region. The bulge (interior loop) regions will be condensed
        into one node.
        """
        new_vertex = True
        while new_vertex:
            new_vertex = False
            bulges = [k for k in self.defines if k[0] != 'y']

            for (b1, b2) in it.combinations(bulges, r=2):
                if b1 in self.edges and b2 in self.edges and self.edges[b1] == self.edges[b2] and len(self.edges[b1]) == 2:
                    connections = self.connections(b1)

                    all_connections = [sorted((self._get_sides_plus(connections[0], b1)[0],
                                               self._get_sides_plus(connections[0], b2)[0])),
                                       sorted((self._get_sides_plus(connections[1], b1)[0],
                                               self._get_sides_plus(connections[1], b2)[0]))]

                    if all_connections == [[1, 2], [0, 3]]:
                        self._merge_vertices([b1, b2])
                        new_vertex = True

    def _relabel_nodes(self):
        """
        Change the labels of the nodes to be more indicative of their nature.
        s: stem
        h: hairpin
        i: interior loop
        m: multiloop
        f: five-prime unpaired
        t: three-prime unpaired
        """
        stems = []
        hairpins = []
        interior_loops = []
        multiloops = []
        fiveprimes = []
        threeprimes = []

        for d in self.defines.keys():
            if d[0] == 'y' or d[0] == 's':
                stems += [d]
                continue

            if len(self.defines[d]) == 0 and len(self.edges[d]) == 1:
                hairpins += [d]
                continue

            if len(self.defines[d]) == 0 and len(self.edges[d]) == 2:
                multiloops += [d]
                continue

            if len(self.edges[d]) <= 1 and self.defines[d][0] == 0:
                fiveprimes += [d]
                continue

            if len(self.edges[d]) == 1 and self.defines[d][1] == self.seq_length - 1:
                threeprimes += [d]
                continue

            if (len(self.edges[d]) == 1 and self.defines[d][0] != 0 and self.defines[d][1] != self.seq_length - 1):
                hairpins += [d]
                continue

            if d[0] == 'm' or (d[0] != 'i' and len(self.edges[d]) == 2 and self.weights[d] == 1 and 
                               self.defines[d][0] != 0 and self.defines[d][1] != self.seq_length - 1):
                multiloops += [d]
                continue

            if d[0] == 'i' or self.weights[d] == 2:
                interior_loops += [d]

        stems.sort(key=self._compare_stems)
        hairpins.sort(key=self._compare_hairpins)
        multiloops.sort(key=self._compare_bulges)
        interior_loops.sort(key=self._compare_stems)

        if fiveprimes:
            d, = fiveprimes
            self._relabel_node(d, 'f0')
        if threeprimes:
            d, = threeprimes
            self._relabel_node(d, 't0')
        for i, d in enumerate(stems):
            self._relabel_node(d, 's%d' % (i))
        for i, d in enumerate(interior_loops):
            self._relabel_node(d, 'i%d' % (i))
        for i, d in enumerate(multiloops):
            self._relabel_node(d, 'm%d' % (i))
        for i, d in enumerate(hairpins):
            self._relabel_node(d, 'h%d' % (i))

    def _remove_degenerate_nodes(self):
        """
        Remove all hairpins that have no length.
        """
        to_remove = []
        for d in self.defines:
            if d[0] == 'h' and len(self.defines[d]) == 0:
                to_remove += [d]

        for r in to_remove:
            self._remove_vertex(r)

    def _define_a_nonzero(self, elem):
        """
        Get a define including the adjacent nucleotides.
        """
        define = self.defines[elem]
        new_def = []
        for i in range(0, len(define), 2):
            new_def.append(max(define[i] - 1, 1))
            new_def.append(min(define[i + 1] + 1, self.seq_length))
        return new_def
    
    def _define_a_zerolength(self, elem):
        """
        Return the define with adjacent nucleotides for a zero-length element.
        Hereby we define that in cases of ambiuigity, the alphabetically first
        zero-length element comes at the lowest nucleotide position etc.
        :param elem: An element, e.g. "m0
        """
        if self.defines[elem] != []:
            raise ValueError("{} does not have zero length".format(elem))
        edges = self.edges[elem]
        if len(edges) == 1:  # Hairpin
            stem, = edges
            define = self.defines[stem]
            if define[2] == define[1] + 1:
                return [define[1], define[2]]
            raise GraphIntegrityError("Very strange zero-length hairpin {} "
                                      "(not?) connected to {}".format(elem, stem))
        if len(edges) == 2:
            stem1, stem2 = edges
            # See if this is the only element connecting the two stems.
            connections = self.edges[stem1] & self.edges[stem2]
            zl_connections = []
            for conn in connections:
                if self.defines[conn] == []:
                    zl_connections.append(conn)
            assert elem in zl_connections
            # We define the 0-length connections to be sorted alphabetically by position
            zl_connections.sort()
            zl_coordinates = self._zerolen_defines_a_between(stem1, stem2)
            if len(zl_connections) != len(zl_coordinates):
                raise GraphIntegrityError("Expecting stems {} and {} to have {} zero-length connections " 
                                          "at nucleotide positions {}, however, found {} elements: {}"
                                          .format(stem1, stem2, len(zl_coordinates), zl_coordinates, len(zl_connections), zl_connections))
            zl_coordinates = list(zl_coordinates)
            zl_coordinates.sort()
            i = zl_connections.index(elem)
            return list(zl_coordinates[i])
        raise GraphIntegrityError("Very strange zero length bulge {} with more than 2 adjacent "
                                  "elements: {}.".format(elem, edges))

    def _zerolen_defines_a_between(self, stem1, stem2):
        zl_coordinates = set()
        for k, l in it.product(range(4), repeat=2):
            if abs(self.defines[stem1][k] - self.defines[stem2][l]) == 1:
                d = [self.defines[stem1][k], self.defines[stem2][l]]
                d.sort()
                zl_coordinates.add(tuple(d))
        return zl_coordinates

    def _get_sides_plus(self, s1, bulge):
        """
        Get the side of s1 that is next to b.
        s1e -> s1b -> b
        :param s1: The stem.
        :param b: The bulge.
        :return: A tuple indicating the corner of the stem that connects
                 to the bulge as well as the corner of the bulge that connects
                 to the stem.
                 These sides are equivalent to the indices of the define.
        """
        if bulge not in self.edges[s1]:
            raise ValueError("_get_sides_plus expects stem to be connected to bulge!")

        s1d = self.defines[s1]
        bd = self.defines[bulge] 

        if not bd:
            bd = self._define_a_zerolength(bulge)
            bd[0] += 1
            bd[1] -= 1

        # before the stem on the 5' strand
        if s1d[0] - bd[1] == 1:
            return (0, 1)
        # after the stem on the 5' strand
        if bd[0] - s1d[1] == 1:
            return (1, 0)
        # before the stem on the 3' strand
        if s1d[2] - bd[1] == 1:
            return (2, 1)
        # after the stem on the 3' strand
        if bd[0] - s1d[3] == 1:
            return (3, 0)

        raise GraphIntegrityError("Faulty bulge {}:{} connected to {}:{}".format(bulge, bd, s1, s1d))

    def _relabel_node(self, old_name, new_name):
        """
        Change the name of a node.
        param old_name: The previous name of the node
        param new_name: The new name of the node
        """
        # replace the define name
        define = self.defines[old_name]

        del self.defines[old_name]
        self.defines[new_name] = define

        # replace the index into the edges array
        edge = self.edges[old_name]
        del self.edges[old_name]
        self.edges[new_name] = edge

        # replace the name of any edge that pointed to old_name
        for k in self.edges.keys():
            new_edges = set()
            for e in self.edges[k]:
                if e == old_name:
                    new_edges.add(new_name)
                else:
                    new_edges.add(e)
            self.edges[k] = new_edges


def _cleaned_graph(rna_g):
    """
    Return a new RNA_Graph with a cleaned secondary structure
    or a reference to the original RNA_Graph.
    This function performs the cleaning operation: 
    all stems with only a single basepair are removed from the structure.

    :returns: A modified copy of or a reference to the original RNA_Graph
    """
    bps_to_remove = []       
    bps_to_remove.extend(rna_g.length_one_stem_basepairs())
    if bps_to_remove:
        new_tuples = rna_g.to_pair_tuples(bps_to_remove)
        return type(rna_g)(new_tuples, rna_g.seq_length)
    else:
        return rna_g

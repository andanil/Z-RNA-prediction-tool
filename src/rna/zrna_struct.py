from .rna_struct import RNAStructure


class ZRNAStructure(RNAStructure):
    def __init__(self, seq, snp_pos):
        super(ZRNAStructure, self).__init__(seq)
        self.snp_pos = snp_pos
        self.zh_score = None

    def __gt__(self, other):
        return self.zh_score > other.zh_score

    def __lt__(self, other):
        self.zh_score < other.zh_score

    def __eq__(self, other):
        self.zh_score == other.zh_score

    def __str__(self):
        return f'ZH-score: {self.zh_score}\nLength: {self.len()}\nSeq: {self._seq}\nStructure (forward): {self._forward_struct}\nStructure (reversed): {self._reversed_struct}'

    def set_zhscore(self, zh_score):
        self.zh_score = zh_score

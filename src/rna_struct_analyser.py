import numpy as np


from .utils import get_seq_around_pos, reverse_compl, replace
from .rna.zrna_struct import ZRNAStructure
from .score_services.zhunt_score_service import ZHuntScoreService


def get_structs_around_snp(genome, chrom, pos, windows):
    struct_list = []
    for window in windows:
        struct_list += [get_zrna_struct(get_seq_around_pos(genome, chrom, pos, window), window[0])]
    struct_list.sort(key=lambda x: x.zh_score, reverse=True)
    return struct_list


def get_structs_around_edited_snp(genome, chrom, pos, modif, windows):
    struct_list = []
    for window in windows:
        seq = get_seq_around_pos(genome, chrom, pos, window)
        seq = replace(seq, window[0], modif)
        struct_list += [get_zrna_struct(seq, window[0])]
    struct_list.sort(key=lambda x: x.zh_score, reverse=True)
    return struct_list


def get_edited_struct(genome, chrom, center, modif_pos, modif, window):
    seq = get_seq_around_pos(genome, chrom, center, window)
    seq = replace(seq, modif_pos, modif)
    return get_zrna_struct(seq, window[0])


def get_zrna_struct(seq, snp_pos):
    struct = ZRNAStructure(seq, snp_pos)
    score_service = ZHuntScoreService()
    struct.set_zhscore(calculate_zscore(seq, struct, score_service))
    return struct


def calculate_zscore(seq, struct, score_service):
    bp = struct.get_bp()
    z_score = score_service.calculate_scores(seq)[np.array(bp).flatten()].sum() if bp else 0
    if score_service.with_reversed:
        bp = struct.get_bp(forward=False)
        z_score += score_service.calculate_scores(reverse_compl(seq))[np.array(bp).flatten()].sum() if bp else 0
    return z_score

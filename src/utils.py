def get_seq_around_pos(genome, chrm, pos, window):
    return genome.fetch(chrm, (pos - window[0] - 1), (pos + window[1] - 1)).upper()  # -1 because start and end denote 0-based


def get_seq(genome, chrm, start, end):
    return genome.fetch(chrm, start, end).upper()


def reverse_compl(seq):
    nn = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A', 'N': 'N'}
    return "".join(nn[n] for n in reversed(seq))

def from_chrom_to_seq_coord_single(x, pos_chrom, pos_seq):
    return pos_seq - pos_chrom + x


def from_chrom_to_seq_coord(region_chrom, pos_chrom, pos_seq, length):
    return (max(0, pos_seq - pos_chrom + region_chrom[0]), min(pos_seq - pos_chrom + region_chrom[1], length - 1))


def from_chrom_to_seq_coord_padding_based(region_chrom, padding):
    return (padding, padding + (region_chrom[1] - region_chrom[0]))


def replace(s, position, character):
    return (s[:position] + character + s[position+1:]).upper()


def reverse_region(region, length):
    reversed_ind = [i for i in range(length - 1, -1, -1)]
    return (reversed_ind[region[1]], reversed_ind[region[0]])


def reverse_pos(pos, length):
    reversed_ind = [i for i in range(length - 1, -1, -1)]
    return reversed_ind[pos]

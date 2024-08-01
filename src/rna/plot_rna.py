import RNA
import numpy as np
import math

import matplotlib.pyplot as plt
import matplotlib.colors as mc
from matplotlib.backends.backend_pdf import PdfPages
import colorsys

from .utils import vec_distance, magnitude
from ..utils import reverse_region, reverse_pos


def save_image(filename):
    p = PdfPages(filename)
    for n in plt.get_fignums():
        plt.figure(n).savefig(p, format='pdf')
    p.close()
    plt.close('all')


def plot_rna(rna_struct, title='', snp_pos=None, zdna_coord=None, forward=True):
    RNA.cvar.rna_plot_type = 1

    bp_to_color = {'G': 'red',
                   'A': 'yellow',
                   'U': 'blue',
                   'T': 'blue',
                   'C': 'green'}
    lighten = 0.5

    fig_size = 10 * (len(rna_struct) / 100 + 1)
    _, ax = plt.subplots(figsize=(fig_size, fig_size))

    coords = []
    vrna_coords = RNA.get_xy_coordinates(rna_struct.get_struct(forward))
    for i in range(len(rna_struct)):
        coords.append((vrna_coords.get(i).X, vrna_coords.get(i).Y))
    coords = np.array(coords)

    bkwargs = {"color": "black", "zorder": 0, "linewidth": 1}
    ax.plot(coords[:, 0], coords[:, 1], **bkwargs)

    basepairs = []
    rna_struct.reset_graph_type(forward)
    for s in rna_struct.get_stems():
        for p1, p2 in rna_struct.get_stem_bp(s):
            basepairs.append([coords[p1], coords[p2]])

    if basepairs:
        basepairs = np.array(basepairs)
        bpkwargs = {"color": "black", "zorder": 0, "linewidth": 3}
        ax.plot(basepairs[:, :, 0].T, basepairs[:, :, 1].T, **bpkwargs)

    if not forward:
        zdna_coord = reverse_region(zdna_coord, len(rna_struct))

    if snp_pos:
        snp_pos = snp_pos if forward else reverse_pos(snp_pos, len(rna_struct))

    for i, coord in enumerate(coords):
        r = 8 if snp_pos and i == snp_pos else 5
        c = bp_to_color[rna_struct.seq[i]]
        h, l, s = colorsys.rgb_to_hls(*mc.to_rgb(c))
        l += (1 - l) * min(1, lighten)
        c = colorsys.hls_to_rgb(h, l, s)
        if zdna_coord:
            ed_color = 'red' if i >= zdna_coord[0] and i <= zdna_coord[1] else None
        else:
            ed_color = None
        circle = plt.Circle((coord[0], coord[1]), facecolor=c, radius=r,
                            edgecolor=ed_color, lw=4)
        ax.add_artist(circle)

        size = 18 if snp_pos and i == snp_pos else 10
        fw = 'bold' if snp_pos and i == snp_pos else 'black'
        ax.annotate(rna_struct.seq[i], xy=coord, ha="center", va="center",
                    fontweight=fw, fontsize=size)

    all_coords = list(coords)
    for nt in range(10, len(rna_struct), 10):
        annot_pos = _find_annot_pos_on_circle(nt, all_coords, rna_struct)
        if annot_pos is not None:
            ax.annotate(
                str(nt), xy=coords[nt-1], xytext=annot_pos,
                arrowprops={"width": 1, "headwidth": 1, "color": "gray"},
                ha="center", va="center", zorder=0, color='gray',
                fontweight='black')
            all_coords.append(annot_pos)

    datalim = ((min(list(coords[:, 0]) + [ax.get_xlim()[0]]),
                min(list(coords[:, 1]) + [ax.get_ylim()[0]])),
               (max(list(coords[:, 0]) + [ax.get_xlim()[1]]),
                max(list(coords[:, 1]) + [ax.get_ylim()[1]])))

    ax.set_aspect('equal', 'datalim')
    ax.update_datalim(datalim)
    ax.autoscale_view()
    ax.set_axis_off()

    ax.title.set_text(title)
    ax.title.set_size(20 * 0.5 * (len(rna_struct) / 100 + 1))

    return (ax, coords)


def _find_annot_pos_on_circle(nt, coords, rna_struct):
    for i in range(5):
        for sign in [-1, 1]:
            a = np.pi / 4 * i * sign
            if rna_struct.get_elem(nt)[0] == "s":
                bp = rna_struct.get_pairing_partner(nt)
                anchor = coords[bp-1]
            else:
                anchor = np.mean([coords[nt - 2], coords[nt]], axis=0)
            vec = coords[nt - 1] - anchor
            mag = magnitude(vec)
            if mag == 0:
                continue
            vec = vec / mag
            rotated_vec = np.array(
                [vec[0] * math.cos(a) - vec[1] * math.sin(a),
                 vec[0] * math.sin(a) + vec[1] * math.cos(a)])
            annot_pos = coords[nt - 1] + rotated_vec * 18
            if _clashfree_annot_pos(annot_pos, coords):
                return annot_pos
    return None


def _clashfree_annot_pos(pos, coords):
    for c in coords:
        dist = vec_distance(c, pos)
        if dist < 14:
            return False
    return True

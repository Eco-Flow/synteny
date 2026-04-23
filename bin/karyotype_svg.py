#!/usr/bin/env python3
"""Drop-in replacement for `python -m jcvi.graphics.karyotype`.
Removes the nseqids < 5 guard so chromosome labels are shown even when a
species has fewer than 5 chromosomes (JCVI's default behaviour hides them).
Accepts all the same arguments as the original command, including --format svg.
"""
import sys
import jcvi.graphics.karyotype as _k


def _draw(self, chrstyle="auto", keep_chrlabels=False, plot_label=True,
          plot_circles=True, pad=0.03, vpad=0.09):
    if self.empty:
        return
    y = self.y
    color = self.color
    ax = self.ax
    xstart = self.xstart
    gap = self.gap
    va = self.va
    nseqids = len(self.seqids)
    tr = self.tr
    for i, sid in enumerate(self.seqids):
        size = self.sizes[sid]
        rsize = self.ratio * size
        xend = xstart + rsize
        hc = _k.HorizontalChromosome(
            ax, xstart, xend, y,
            height=self.height, lw=self.lw, fc=color, style=chrstyle,
        )
        hc.set_transform(tr)
        si = sid if keep_chrlabels else _k.make_circle_name(sid, self.rev)
        xx = (xstart + xend) / 2
        xstart = xend + gap
        step = 2 if nseqids <= 40 else 10
        if nseqids >= 2 * _k.MaxSeqids and (i + 1) % step != 0:
            continue
        # nseqids < 5 guard removed — always label chromosomes
        hpad = -pad if va == "bottom" else pad
        if plot_circles:
            _k.TextCircle(ax, xx, y + hpad, si, fc="w", color=color, size=10, transform=tr)
    label = _k.markup(self.label)
    c = color if color != "gainsboro" else "k"
    if plot_label:
        if self.label_va == "top":
            x, y = self.x, self.y + vpad
        elif self.label_va == "bottom":
            x, y = self.x, self.y - vpad
        else:
            x, y = self.xstart - vpad / 2, self.y
        ax.text(x, y, label, ha="center", va="center", color=c, transform=tr)


_k.Track.draw = _draw

from jcvi.graphics.karyotype import main
main(sys.argv[1:])

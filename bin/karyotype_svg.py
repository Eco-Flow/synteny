#!/usr/bin/env python3
"""Drop-in replacement for `python -m jcvi.graphics.karyotype`.

Changes from jcvi default:
  - Always draws chromosome labels, regardless of nseqids count (removes both
    the nseqids < 5 and the nseqids >= 2*MaxSeqids skip guards).
  - For dense tracks (nseqids > MaxSeqids), labels are staggered on alternating
    rows and scaled slightly smaller so neighbouring labels don't overlap.
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

    # For dense tracks, scale circles down and stagger labels on two rows so
    # adjacent chromosome labels sit at different heights instead of overlapping.
    dense = nseqids > _k.MaxSeqids
    if dense:
        scale = max(0.5, _k.MaxSeqids / nseqids)
        label_size = max(5, round(10 * scale))
        label_radius = max(0.01, 0.02 * scale)
    else:
        label_size = 10
        label_radius = 0.02

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

        # Stagger: alternate between 1× and 1.8× pad so adjacent labels sit on
        # different rows, halving the effective horizontal label density.
        base_hpad = -pad if va == "bottom" else pad
        hpad = base_hpad * (1.8 if (dense and i % 2 == 1) else 1.0)

        if plot_circles:
            _k.TextCircle(ax, xx, y + hpad, si, fc="w", color=color,
                          size=label_size, radius=label_radius, transform=tr)

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

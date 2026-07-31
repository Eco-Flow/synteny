# Local-only image for the SYNGRAPH process (modules/local/algo/syngraph.nf).
# Not published to quay.io yet -- build and tag it yourself before running
# with -profile docker:
#   docker build --platform linux/amd64 -t synteny-syngraph:local containers/syngraph/
#
# Syngraph (https://github.com/A-J-F-Mackintosh/syngraph) targets Python 3.8
# with fairly old pinned dependencies (networkx==2.4, numpy==1.20.3,
# graph-tool, ete3, pygraphviz...). Their yaml/syngraph.yaml is a full
# `conda env export` snapshot pinned to exact build hashes from ~2022, which
# no longer solves against current channel indices -- this instead follows
# the looser install command documented in their README, which lets conda
# pick currently-available builds satisfying the same version constraints.
# Only linux-64 (amd64) builds of graph-tool etc. are available, hence the
# explicit --platform above (needed when building on Apple Silicon).
# Deliberately hardcoded (not $TARGETPLATFORM): graph-tool and friends are
# only reliably available as linux-64 conda builds, so this always targets
# amd64 regardless of the host/build platform (Docker will emulate on
# arm64 hosts such as Apple Silicon). Ignore the linter's "should not use a
# constant" hint here -- that's exactly what's wanted.
FROM --platform=linux/amd64 continuumio/miniconda3:24.7.1-0

RUN git clone --depth 1 https://github.com/A-J-F-Mackintosh/syngraph.git /opt/syngraph

RUN conda create -y -n syngraph -c etetoolkit -c conda-forge -c bioconda -c defaults \
        python=3.8 \
        pandas \
        docopt \
        ete3 \
        pygraphviz \
        matplotlib \
        tqdm \
        "networkx=2.4" \
        "numpy=1.20.3" \
        "pillow<10" \
        more-itertools \
    && conda clean -afy

ENV PATH="/opt/conda/envs/syngraph/bin:/opt/syngraph:${PATH}"

# The repo's main entrypoint script is `syngraph` at the repo root (invoked as
# `syngraph <module> ...` in modules/local/algo/syngraph.nf) -- confirm this
# still matches the vendored commit and adjust PATH/permissions if it moves.
RUN chmod +x /opt/syngraph/syngraph

WORKDIR /data

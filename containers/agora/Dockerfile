# Local-only image for the AGORA_PREP and AGORA processes
# (modules/local/algo/agora_prep.nf, modules/local/algo/agora.nf).
# Not published to quay.io yet -- build and tag it yourself before running
# with -profile docker:
#   docker build -t synteny-agora:local containers/agora/
#
# AGORA (https://github.com/DyogenIBENS/Agora) ships its own multi-stage
# Dockerfile pinned to python:3.7.7-slim + requirements.txt (just `psutil`).
# This mirrors that base image/dependency rather than re-deriving it, using a
# single stage since a slim final image isn't a concern for a local-only build.
FROM python:3.7.7-slim

# python:3.7.7-slim is based on Debian buster, which is EOL and has been
# moved to archive.debian.org -- repoint apt there before installing git.
RUN sed -i \
        -e 's|deb.debian.org/debian |archive.debian.org/debian |g' \
        -e 's|security.debian.org/debian-security|archive.debian.org/debian-security|g' \
        -e '/buster-updates/d' \
        /etc/apt/sources.list \
    && apt-get -o Acquire::Check-Valid-Until=false update \
    && apt-get install -y --no-install-recommends git procps \
    && rm -rf /var/lib/apt/lists/*

RUN git clone --depth 1 https://github.com/DyogenIBENS/Agora.git /opt/agora \
    && pip install --no-cache-dir -r /opt/agora/requirements.txt

RUN chmod +x /opt/agora/src/*.py
ENV PATH="/opt/agora/src:${PATH}"

# Confirm agora-generic.py's exact positional CLI signature against the
# vendored commit's doc/HowTo.md (see the comment in
# modules/local/algo/agora.nf) before relying on it in production.
WORKDIR /data

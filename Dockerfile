# FROM mambaorg/micromamba:git-4073267-bionic
FROM mambaorg/micromamba@sha256:59fe75f753f5f67d61331575a5ddb342b1882680ffd2085dcf9d66b63369fc1d

ENV PYTHONPATH="/data/src"
ARG QUARTO_VERSION=1.3.361
ARG QUARTO_RELEASE=https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.deb
USER root
RUN apt-get update && \
  apt-get install -y curl && \
  curl -fsSL ${QUARTO_RELEASE} -o quarto.deb && \
  dpkg -i quarto.deb
USER mambauser

COPY --chown=$MAMBA_USER:$MAMBA_USER env.yml /tmp/env-full.yml

# use PYTHONPATH not pip for local package
RUN head -n -3 /tmp/env-full.yml > /tmp/env.yml
# pretty prompt
RUN echo "PS1='🐳 \[\033[1;36m\]\h \[\033[1;34m\]\W\[\033[0;35m\] \[\033[1;93m\]# \[\033[0m\]'" > /home/mambauser/.bashrc

RUN micromamba install -y -n base -f /tmp/env.yml && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1  # (otherwise python will not be found)
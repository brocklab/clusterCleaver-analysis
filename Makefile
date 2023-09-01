SHELL = bash
MAMBA := $(if \
	$(shell command -v micromamba),micromamba,\
	$(if $(shell command -v mamba),mamba,conda))

locked.yml: env.yml
	$(MAMBA) env create -p ./env -f $<
	$(MAMBA) env -p ./env export | sed -r '/(prefix|name):.*/d' > $@
	echo -e '- pip:\n  - -e .' >> $@

env: locked.yml ## generate conda enviroment
	$(MAMBA) env create -p ./env -f $<

dag.svg:
	snakemake --dag --forceall | dot -Tsvg > dag.svg

-include .task.cfg.mk

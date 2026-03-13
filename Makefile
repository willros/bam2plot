PIXI ?= $(if $(wildcard ./.pixi-home/bin/pixi),./.pixi-home/bin/pixi,pixi)
PIXI_CACHE_DIR ?= $(CURDIR)/.pixi-cache
PIXI_MANIFEST ?= $(CURDIR)/pixi.toml

PYTHON_FILES = setup.py bam2plot/__init__.py bam2plot/main.py tests/*.py

.PHONY: clean build upload format format-check test

clean:
	rm -rf dist/ build/ *.egg-info
build:
	python -m build

upload:
	twine upload dist/*
format:
	PIXI_PROJECT_MANIFEST=$(PIXI_MANIFEST) PIXI_CACHE_DIR=$(PIXI_CACHE_DIR) $(PIXI) run format
format-check:
	PIXI_PROJECT_MANIFEST=$(PIXI_MANIFEST) PIXI_CACHE_DIR=$(PIXI_CACHE_DIR) $(PIXI) run format-check
test:
	PIXI_PROJECT_MANIFEST=$(PIXI_MANIFEST) PIXI_CACHE_DIR=$(PIXI_CACHE_DIR) $(PIXI) run test

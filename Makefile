PIXI ?= $(if $(wildcard ./.pixi-home/bin/pixi),./.pixi-home/bin/pixi,pixi)
PIXI_CACHE_DIR ?= $(CURDIR)/.pixi-cache
PIXI_MANIFEST ?= $(CURDIR)/pixi.toml

.PHONY: clean build upload test

clean:
	rm -rf dist/ build/ *.egg-info
build:
	python -m build

upload:
	twine upload dist/*
test:
	PIXI_PROJECT_MANIFEST=$(PIXI_MANIFEST) PIXI_CACHE_DIR=$(PIXI_CACHE_DIR) $(PIXI) run test

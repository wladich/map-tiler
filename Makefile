define help
Available targets:
	help
	black
	check
	venv - create or update venv for development
	clean - clean caches
endef
export help

TOOL_PREFIX=./.venv/bin/

help:
	@echo "$$help"

black:
	$(TOOL_PREFIX)black --diff --check -q .

check: black
	@echo All checks passed.

.PHONY: venv
venv:
	uv sync

clean:
	rm -r ./.venv ./__pycache__

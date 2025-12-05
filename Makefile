# =============================================================================
# intronIC Development Makefile
# =============================================================================
#
# This Makefile provides common development commands for intronIC.
# Run `make help` to see all available targets.
#
# Prerequisites:
#   - pixi (recommended): https://pixi.sh/
#   - uv (optional): https://docs.astral.sh/uv/
#
# Quick start:
#   make install    # Set up development environment
#   make test       # Run tests
#   make help       # Show all commands
#
# =============================================================================

# Default shell
SHELL := /bin/bash

# Detect available tools
PIXI := $(shell command -v pixi 2> /dev/null)
UV := $(shell command -v uv 2> /dev/null)

# Default runner (prefer pixi if available)
ifdef PIXI
    RUNNER := pixi run
    PYTHON := pixi run python
    PYTEST := pixi run pytest
else ifdef UV
    RUNNER := 
    PYTHON := python
    PYTEST := pytest
else
    RUNNER := 
    PYTHON := python
    PYTEST := pytest
endif

# =============================================================================
# HELP
# =============================================================================

.DEFAULT_GOAL := help

.PHONY: help
help:  ## Show this help message
	@echo ""
	@echo "intronIC Development Commands"
	@echo "=============================="
	@echo ""
	@echo "Usage: make <target>"
	@echo ""
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  %-18s %s\n", $$1, $$2}'
	@echo ""

# =============================================================================
# SETUP & INSTALLATION
# =============================================================================

.PHONY: install
install:  ## Install all dependencies (uses pixi if available, else uv)
ifdef PIXI
	@echo "Installing with pixi..."
	pixi install
else ifdef UV
	@echo "Installing with uv..."
	uv sync
else
	@echo "Installing with pip..."
	pip install -e ".[dev]"
endif
	@echo "✓ Installation complete"

.PHONY: install-pixi
install-pixi:  ## Install dependencies using pixi
	@echo "Installing with pixi..."
	pixi install
	@echo "✓ Pixi installation complete"

.PHONY: install-uv
install-uv:  ## Install dependencies using uv
	@echo "Installing with uv..."
	uv sync
	@echo "✓ uv installation complete"

.PHONY: install-pip
install-pip:  ## Install dependencies using pip (editable mode)
	@echo "Installing with pip..."
	pip install -e ".[dev]"
	@echo "✓ pip installation complete"

# =============================================================================
# DEPENDENCY MANAGEMENT
# =============================================================================

.PHONY: lock
lock:  ## Update all lock files (run after changing pyproject.toml or pixi.toml)
	@echo "Updating lock files..."
ifdef UV
	@echo "  → Running uv lock..."
	uv lock
else
	@echo "  → uv not found, skipping uv.lock"
endif
ifdef PIXI
	@echo "  → Running pixi install..."
	pixi install
else
	@echo "  → pixi not found, skipping pixi.lock"
endif
	@echo "✓ Lock files updated"
	@echo ""
	@echo "Remember to commit: uv.lock pixi.lock"

.PHONY: lock-uv
lock-uv:  ## Update only uv.lock
	@echo "Updating uv.lock..."
	uv lock
	@echo "✓ uv.lock updated"

.PHONY: lock-pixi
lock-pixi:  ## Update only pixi.lock
	@echo "Updating pixi.lock..."
	pixi install
	@echo "✓ pixi.lock updated"

.PHONY: upgrade
upgrade:  ## Upgrade all dependencies to latest compatible versions
	@echo "Upgrading dependencies..."
ifdef UV
	uv lock --upgrade
endif
ifdef PIXI
	pixi update
endif
	@echo "✓ Dependencies upgraded"

# =============================================================================
# CODE QUALITY
# =============================================================================

.PHONY: format
format:  ## Format code with black and isort
	@echo "Formatting code..."
	$(RUNNER) black src/ tests/
	$(RUNNER) isort src/ tests/
	@echo "✓ Code formatted"

.PHONY: lint
lint:  ## Run all linters (ruff)
	@echo "Running linters..."
	$(RUNNER) ruff check src/ tests/
	@echo "✓ Linting complete"

.PHONY: lint-fix
lint-fix:  ## Run linters and auto-fix issues
	@echo "Running linters with auto-fix..."
	$(RUNNER) ruff check --fix src/ tests/
	@echo "✓ Linting complete"

.PHONY: typecheck
typecheck:  ## Run type checking with mypy
	@echo "Running type checker..."
	$(RUNNER) mypy src/
	@echo "✓ Type checking complete"

.PHONY: check
check: lint typecheck  ## Run all code quality checks (lint + typecheck)
	@echo "✓ All checks passed"

# =============================================================================
# TESTING
# =============================================================================

.PHONY: test
test:  ## Run full test suite
	@echo "Running tests..."
	$(PYTEST) tests/ -v
	@echo "✓ Tests complete"

.PHONY: test-quick
test-quick:  ## Run quick tests only (skip slow integration tests)
	@echo "Running quick tests..."
	$(PYTEST) tests/unit/ -v --tb=short
	@echo "✓ Quick tests complete"

.PHONY: test-unit
test-unit:  ## Run unit tests only
	@echo "Running unit tests..."
	$(PYTEST) tests/unit/ -v
	@echo "✓ Unit tests complete"

.PHONY: test-integration
test-integration:  ## Run integration tests only
	@echo "Running integration tests..."
	$(PYTEST) tests/integration/ -v
	@echo "✓ Integration tests complete"

.PHONY: test-cov
test-cov:  ## Run tests with coverage report
	@echo "Running tests with coverage..."
	$(PYTEST) tests/ -v --cov=src/intronIC --cov-report=term-missing --cov-report=html
	@echo "✓ Coverage report generated in htmlcov/"

.PHONY: test-small
test-small:  ## Run small test on Chr19 data (pixi only)
ifdef PIXI
	@echo "Running small test..."
	pixi run test-small
else
	@echo "test-small requires pixi"
endif

# =============================================================================
# BUILD & RELEASE
# =============================================================================

.PHONY: build
build:  ## Build distribution packages (wheel and sdist)
	@echo "Building distribution packages..."
	$(PYTHON) -m build
	@echo "✓ Build complete - see dist/"

.PHONY: build-check
build-check:  ## Check built packages with twine
	@echo "Checking built packages..."
	$(RUNNER) twine check dist/*
	@echo "✓ Package check complete"

.PHONY: clean
clean:  ## Remove build artifacts and caches
	@echo "Cleaning build artifacts..."
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info/
	rm -rf src/*.egg-info/
	rm -rf .pytest_cache/
	rm -rf .mypy_cache/
	rm -rf .ruff_cache/
	rm -rf htmlcov/
	rm -rf .coverage
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete 2>/dev/null || true
	@echo "✓ Clean complete"

.PHONY: clean-all
clean-all: clean  ## Remove all artifacts including environments
	@echo "Removing environments..."
	rm -rf .venv/
	rm -rf .pixi/
	@echo "✓ Full clean complete"

# =============================================================================
# DEVELOPMENT HELPERS
# =============================================================================

.PHONY: run
run:  ## Run intronIC with arguments (use: make run ARGS="-h")
	$(RUNNER) intronIC $(ARGS)

.PHONY: shell
shell:  ## Open a Python shell in the development environment
ifdef PIXI
	pixi shell
else
	$(PYTHON)
endif

.PHONY: version
version:  ## Show current version
	@$(PYTHON) -c "import intronIC; print(intronIC.__version__)"

.PHONY: info
info:  ## Show development environment info
	@echo ""
	@echo "Development Environment"
	@echo "========================"
	@echo "Python: $$($(PYTHON) --version 2>&1)"
	@echo "Runner: $(RUNNER)"
ifdef PIXI
	@echo "Pixi: $$(pixi --version)"
endif
ifdef UV
	@echo "uv: $$(uv --version)"
endif
	@echo ""
	@echo "Package Info"
	@$(PYTHON) -c "import intronIC; print('intronIC:', intronIC.__version__)" 2>/dev/null || echo "intronIC: not installed"

# =============================================================================
# WORKFLOW SHORTCUTS
# =============================================================================

.PHONY: pre-commit
pre-commit: format lint test-quick  ## Run all pre-commit checks
	@echo "✓ Pre-commit checks passed"

.PHONY: ci
ci: install lint test  ## Run full CI pipeline
	@echo "✓ CI pipeline complete"

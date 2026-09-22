PYTHON ?= python3
SAGE ?= sage
export PYTHONPATH := src

.PHONY: check backend small standard acceptance scaling figures demo forge sage-check clean-results

check:
	$(PYTHON) scripts/check_environment.py

backend:
	sh scripts/build_backend.sh

small: backend
	$(PYTHON) scripts/run_small_validation.py

standard: backend
	$(PYTHON) scripts/run_standard_128.py --seed 0

acceptance: backend
	$(PYTHON) scripts/run_acceptance.py --source kat

figures:
	$(PYTHON) scripts/visualize.py

scaling: backend
	$(PYTHON) scripts/measure_structure_scaling.py

demo: backend
	$(PYTHON) scripts/make_attack_demo.py

forge: backend
	$(PYTHON) scripts/forge_server.py

sage-check: backend
	$(SAGE) -python scripts/verify_with_sage.py

clean-results:
	rm -rf results/acceptance figures .mplconfig

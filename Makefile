.PHONY: all
all: install

#* Installation
.PHONY: install
install:
	poetry lock -n && poetry export --without-hashes > requirements.txt
	poetry install -n
	#* -poetry run mypy --install-types --non-interactive ./

#* Formatters
.PHONY: codestyle
codestyle:
	#* poetry run pyupgrade --exit-zero-even-if-changed --py39-plus **/*.py
	poetry run isort --settings-path pyproject.toml ./
	poetry run black --config pyproject.toml ./
	poetry run pylint pinr tests --rcfile=.pylint.ini --extension-pkg-whitelist='pydantic'
	poetry run pydocstyle pinr tests

#* Python sehll
.PHONY: shell
shell:
	poetry run ipython
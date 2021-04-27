# rxn_rebuild

[![Anaconda-Server Badge](https://anaconda.org/brsynth/rxn_rebuild/badges/latest_release_date.svg)](https://anaconda.org/brsynth/rxn_rebuild) [![Anaconda-Server Badge](https://anaconda.org/brsynth/rxn_rebuild/badges/version.svg)](https://anaconda.org/brsynth/rxn_rebuild)

Rebuild full reaction from reaction rule by adding co-factors removed to generate the rule.

## Prerequisites

* Python 3

## Install

### Prerequisite

The conda package manager is required. Fresh instructions on how to install conda are [available online](https://docs.conda.io/projects/conda/en/latest/user-guide/install/).

### conda environment

In case a new conda environment `<my_env>` need to be set up, first start with:
```shell
conda create -n my_env python=3
```

### conda package

Install in the `<my_env>` conda environment:
```shell
conda install -c brsynth -c conda-forge -n <my_env> rxn_rebuild 
```

## Run

### rxn_rebuild process
**From CLI**
```sh
python -m rxn_rebuild <rxn_rule_id> <transfo> [<ori_rxn_id>]
```
**From Python code**
```python
from rr_cache import rrCache
from rxn_rebuild import rebuild_rxn, build_args_parser

parser = build_args_parser()
args   = parser.parse_args()

compl_rxn_smi = rebuild_rxn(
    cache = rrCache(db='file', attrs=None),
    rxn_rule_id = args.rxn_rule_id,
    trans_smi = args.trans_smi,
    ori_rxn_id = args.ori_rxn_id
)
```

## Tests
Test can be run with the following commands:

### Natively
```bash
cd tests
pytest -v
```

# CI/CD
For further tests and development tools, a CI toolkit is provided in `ci` folder (see [ci/README.md](ci/README.md)).


<!-- ### How to cite rxn_rebuild? -->

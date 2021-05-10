# rxn_rebuild

[![Anaconda-Server Badge](https://anaconda.org/brsynth/rxn_rebuild/badges/latest_release_date.svg)](https://anaconda.org/brsynth/rxn_rebuild) [![Anaconda-Server Badge](https://anaconda.org/brsynth/rxn_rebuild/badges/version.svg)](https://anaconda.org/brsynth/rxn_rebuild)

Rebuild full reaction from a reaction rule ID and a chemical transformation by adding co-factors removed when the rule has been generated. The algorithm proceeds to the following:
- seeks for the reaction(s) rule(s) of the given reaction rule ID, in the cache
- for each reaction rule
    - seeks for the template reaction of the given reaction rule ID, in the cache
    - find compounds to add by doing the difference between the template reaction and the reaction rule
  
## Input
- reaction rule ID
- chemical transfomation as SMILES (xxx.xxx>>xxx.xxx) or IDs (CMPD_ID_1 + CMPD_ID_2 = CMPD_ID_3)
- (Optional) template reaction ID

## Output
- Prints out the full completed transformation
- Returns a dictionary with the following keys:
  - 'full_transfo': completed transformation (SMILES or IDs)
  - 'added_cmpds': sub-dictionary of compounds to add to the transformation with:
    - keys: IDs
    - values: infos (SMILES, InChI, InChIKey, formula, name)


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

completed_transfos = rebuild_rxn(
    cache = rrCache(db='file', attrs=['rr_reactions', 'template_reactions','cid_strc']),
    rxn_rule_id = args.rxn_rule_id,
    transfo = args.trans_smi,
    tmpl_rxn_id = args.ori_rxn_id
)
```
If `cache` is not provided, it ill be automatically loaded within `rebuild_rxn` function but it could be much slower if called inside a loop.

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

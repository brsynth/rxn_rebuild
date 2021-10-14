# rxn_rebuild
Rebuild full reaction from reaction rule
| Name | Downloads | Version | Platforms |
| --- | --- | --- | --- |
| [![Conda Recipe](https://img.shields.io/badge/recipe-rxn_rebuild-green.svg)](https://anaconda.org/conda-forge/rxn_rebuild) | [![Conda Downloads](https://img.shields.io/conda/dn/conda-forge/rxn_rebuild.svg)](https://anaconda.org/conda-forge/rxn_rebuild) | [![Conda Version](https://img.shields.io/conda/vn/conda-forge/rxn_rebuild.svg)](https://anaconda.org/conda-forge/rxn_rebuild) | [![Conda Platforms](https://img.shields.io/conda/pn/conda-forge/rxn_rebuild.svg)](https://anaconda.org/conda-forge/rxn_rebuild) |

## Description
*rxn_rebuild* provides a cache for RetroRules and MetaNetX compounds and reactions


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


## Install

The conda package manager is required. Fresh instructions on how to install conda are [available online](https://docs.conda.io/projects/conda/en/latest/user-guide/install/).

```shell
conda install -c conda-forge rxn_rebuild 
```

## Run

### rxn_rebuild process
**From CLI**
```sh
python -m rxn_rebuild <rxn_rule_id> <transfo> [<ori_rxn_id>]
```
**From Python code**
```python
from rxn_rebuild import rebuild_rxn, build_args_parser

parser = build_args_parser()
args   = parser.parse_args()

completed_transfos = rebuild_rxn(
    rxn_rule_id=args.rxn_rule_id,
    transfo=args.trans_smi,
    tmpl_rxn_id=args.ori_rxn_id
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

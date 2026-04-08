from argparse import ArgumentParser

DEFAULTS = {
    "cspace": "mnx3.1",
}


def add_arguments(parser: ArgumentParser) -> ArgumentParser:

    parser.add_argument("rxn_rule_id", type=str, help="Reaction rule identifier")
    parser.add_argument(
        "transfo",
        type=str,
        help="Transformation to complete with flatten compounds ids (e.g. MNXM181 + MNXM4 = MNXM1 + MNXM1 + MNXM1144) or in SMILES string (e.g. [H]OC(=O)C([H])=C([H])C([H])=C([H])C(=O)O[H]>>[H]Oc1c([H])c([H])c([H])c([H])c1O[H].O=O.O=O)",
    )
    parser.add_argument(
        "--tmpl_rxn_id", type=str, help="Template (original) reaction identifier"
    )
    parser.add_argument(
        "--to-ignore",
        type=str,
        help="Name of the file containing the list of compounds to ignore (default: None)",
        default=None,
    )
    parser.add_argument(
        "--chemical-space",
        dest="cspace",
        default=DEFAULTS["cspace"],
        type=str,
        help="Chemical space to use (e.g. mnx3.1, mnx4.4...). Determines which configuration files and folders to use both the cache and the input cache (default: %(default)s).",
    )

    return parser

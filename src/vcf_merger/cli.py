"""Command-line interface for vcf-merger."""

from __future__ import annotations

import argparse
import logging
import sys
from typing import Optional

from vcf_merger import __version__
from vcf_merger.exceptions import VcfMergerError
from vcf_merger.harmonizer import harmonize_vcfs
from vcf_merger.inspect import inspect_vcf_json
from vcf_merger.models import AnalysisMode, EnsembleStrategy
from vcf_merger.normalization import normalize_vcf


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="vcf-merger",
        description=(
            "Technical VCF normalization, evidence aggregation, and harmonization. "
            "Does not perform ACMG/AMP clinical classification."
        ),
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Increase logging verbosity",
    )

    sub = parser.add_subparsers(dest="command", required=True)

    p_inspect = sub.add_parser("inspect", help="Inspect a VCF/BCF and print metadata JSON")
    p_inspect.add_argument("vcf", help="Input VCF/VCF.GZ/BCF")
    p_inspect.add_argument(
        "--caller",
        default=None,
        help="Explicit caller name (overrides detection)",
    )

    p_norm = sub.add_parser("normalize", help="Reference-aware normalization via bcftools norm")
    p_norm.add_argument("--reference", "-r", required=True, help="Reference FASTA")
    p_norm.add_argument("--input", "-i", required=True, help="Input VCF")
    p_norm.add_argument("--output", "-o", required=True, help="Output VCF/VCF.GZ")
    p_norm.add_argument(
        "--multiallelic",
        default="-any",
        help="bcftools norm -m value (default: -any)",
    )

    p_merge = sub.add_parser("merge", help="Harmonize per-caller VCFs")
    p_merge.add_argument(
        "--mode",
        choices=[m.value for m in AnalysisMode],
        default=AnalysisMode.GERMLINE.value,
    )
    p_merge.add_argument("--reference", "-r", default=None, help="Reference FASTA (required unless --no-normalize)")
    p_merge.add_argument(
        "--input",
        "-i",
        action="append",
        dest="inputs",
        required=True,
        help="Input VCF (repeatable)",
    )
    p_merge.add_argument(
        "--caller",
        action="append",
        dest="callers",
        default=None,
        help="Explicit caller for each --input (repeatable, same order)",
    )
    p_merge.add_argument(
        "--strategy",
        choices=[s.value for s in EnsembleStrategy],
        default=EnsembleStrategy.UNION.value,
    )
    p_merge.add_argument(
        "--consensus-n",
        type=int,
        default=2,
        help="Minimum distinct callers for consensus strategy",
    )
    p_merge.add_argument(
        "--include-caller",
        action="append",
        default=None,
        help="Only include these callers (repeatable)",
    )
    p_merge.add_argument(
        "--exclude-caller",
        action="append",
        default=None,
        help="Exclude these callers (repeatable)",
    )
    p_merge.add_argument("--output", "-o", required=True, help="Output harmonized VCF/VCF.GZ")
    p_merge.add_argument("--tumor-sample", default=None, help="Tumor sample name (somatic)")
    p_merge.add_argument("--normal-sample", default=None, help="Normal sample name (somatic)")
    p_merge.add_argument(
        "--no-normalize",
        action="store_true",
        help="Skip bcftools normalization (not recommended)",
    )
    p_merge.add_argument(
        "--allow-gvcf",
        action="store_true",
        help="Do not reject gVCF inputs (dangerous; unsupported)",
    )

    return parser


def main(argv: Optional[list[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    level = logging.WARNING
    if args.verbose == 1:
        level = logging.INFO
    elif args.verbose >= 2:
        level = logging.DEBUG
    logging.basicConfig(level=level, format="%(levelname)s: %(message)s")

    try:
        if args.command == "inspect":
            print(inspect_vcf_json(args.vcf, caller=args.caller))
            return 0

        if args.command == "normalize":
            out = normalize_vcf(
                args.input,
                args.output,
                reference=args.reference,
                multiallelic=args.multiallelic,
            )
            print(out)
            return 0

        if args.command == "merge":
            result = harmonize_vcfs(
                args.inputs,
                args.output,
                mode=args.mode,
                strategy=args.strategy,
                reference=args.reference,
                normalize=not args.no_normalize,
                consensus_n=args.consensus_n,
                include_callers=args.include_caller,
                exclude_callers=args.exclude_caller,
                callers=args.callers,
                tumor_sample=args.tumor_sample,
                normal_sample=args.normal_sample,
                allow_gvcf=args.allow_gvcf,
                command_line=["vcf-merger", *sys.argv[1:]],
            )
            print(
                f"Wrote {result['paths']['vcf']} "
                f"({result['variant_count']} variants, mode={result['mode']}, "
                f"strategy={result['strategy']})"
            )
            print(f"Evidence: {result['paths']['evidence']}")
            print(f"Provenance: {result['paths']['provenance']}")
            return 0

        parser.error(f"Unknown command: {args.command}")
        return 2
    except VcfMergerError as exc:
        logging.error("%s", exc)
        return 1
    except Exception as exc:  # noqa: BLE001
        logging.exception("Unexpected error: %s", exc)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

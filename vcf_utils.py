# vcf_utils.py
#
# Shared core logic for VCF merging.
# Used by both Merge_All_VCFs.py and Merge_all_VCF_Groups.py.

import os
import io

import pandas as pd

# -----------------------------------------------------------------
# Constants
# -----------------------------------------------------------------

# Columns that uniquely identify a variant across callers.
# Using REF + ALT ensures multi-allelic sites at the same position
# are never incorrectly collapsed into a single row.
MERGE_KEY = ['#CHROM', 'POS', 'REF', 'ALT']

# Caller priority (index 0 = highest priority).
# When the same variant is found by multiple callers, QUAL / FILTER /
# FORMAT / sample data are taken from the highest-priority caller.
PRIORITY_ORDER = ['.MU.vcf', '.FB.vcf', '.HC.vcf', '.ST.vcf', '.DV.vcf']

# VCF INFO header line for the CALLERS tag we inject into every output row.
_CALLERS_INFO_HDR = (
    '##INFO=<ID=CALLERS,Number=.,Type=String,'
    'Description="Comma-separated list of variant callers that identified this variant">'
)

# Maps filename tags to short caller names used in the CALLERS INFO field.
_CALLER_TAGS = {
    '.MU.': 'MU',   # GATK Mutect2
    '.FB.': 'FB',   # FreeBayes
    '.HC.': 'HC',   # GATK HaplotypeCaller
    '.ST.': 'ST',   # Strelka
    '.DV.': 'DV',   # DeepVariant
}

# Numeric sort key for standard chromosomes (handles both 'chr1' and '1' style).
_CHROM_ORDER = {str(i): i for i in range(1, 23)}
_CHROM_ORDER.update({'X': 23, 'Y': 24, 'MT': 25, 'M': 25})


# -----------------------------------------------------------------
# Internal helpers
# -----------------------------------------------------------------

def _get_caller_name(vcf_path: str) -> str:
    """Return a short caller label (e.g. 'FB') extracted from the filename."""
    basename = os.path.basename(vcf_path)
    for tag, name in _CALLER_TAGS.items():
        if tag in basename:
            return name
    return 'UNKNOWN'


def _chrom_sort_key(chrom_series: pd.Series) -> pd.Series:
    """Convert a CHROM column to a numeric sort key."""
    return (
        chrom_series
        .str.replace('chr', '', regex=False)
        .map(lambda x: _CHROM_ORDER.get(x, 99))
    )


def _merge_meta_headers(all_meta_headers: list[str]) -> list[str]:
    """Deduplicate ## header lines (preserving first occurrence) and inject CALLERS tag."""
    seen: set[str] = set()
    merged: list[str] = []
    for line in all_meta_headers:
        if line not in seen:
            seen.add(line)
            merged.append(line)
    if _CALLERS_INFO_HDR not in seen:
        merged.append(_CALLERS_INFO_HDR)
    return merged


# -----------------------------------------------------------------
# Public API
# -----------------------------------------------------------------

def read_vcf(file: str) -> tuple[list[str], pd.DataFrame]:
    """Read a VCF file.

    Returns
    -------
    meta_headers : list[str]
        All ``##`` meta-information lines (without trailing newline).
    df : pd.DataFrame
        Variant data rows, including the ``#CHROM`` header row as column names.
    """
    meta_headers: list[str] = []
    data_lines: list[str] = []

    with open(file, 'r', encoding='iso-8859-1') as f:
        for line in f:
            if line.startswith('##'):
                meta_headers.append(line.rstrip('\n'))
            else:
                data_lines.append(line)

    df = pd.read_csv(
        io.StringIO(''.join(data_lines)),
        sep='\t',
        dtype=str,        # keep all values as strings to avoid type coercion
        low_memory=False,
    )
    return meta_headers, df


def sort_vcf_files(vcf_files: list[str]) -> list[str]:
    """Sort VCF files by caller priority (MU > FB > HC > ST > DV)."""
    vcf_files.sort(
        key=lambda x: next(
            (i for i, ext in enumerate(PRIORITY_ORDER) if ext in x),
            len(PRIORITY_ORDER),   # unknown callers go last
        )
    )
    return vcf_files


def merge_vcfs(vcf_files: list[str], output_file: str) -> None:
    """Merge multiple per-caller VCF files into a single output VCF.

    Strategy
    --------
    For each unique variant (CHROM + POS + REF + ALT):

    1. Collect every row for that variant from every caller.
    2. Use the **complete row** (QUAL, FILTER, INFO, FORMAT, sample column)
       from the highest-priority caller that called the variant.
       This prevents mixing FORMAT tags between callers, which would produce
       an internally inconsistent FORMAT/sample pair.
    3. Append ``CALLERS=<names>`` to the INFO field so downstream tools can
       see which callers supported the variant.

    All ``##`` meta-information header lines are collected, deduplicated, and
    written at the top of the output file, producing a valid VCF.

    Parameters
    ----------
    vcf_files : list[str]
        Paths to input VCF files.  Will be sorted by caller priority internally.
    output_file : str
        Path for the merged output VCF.
    """
    vcf_files = sort_vcf_files(vcf_files)

    all_meta_headers: list[str] = []
    tagged_frames: list[pd.DataFrame] = []

    for priority_idx, vcf in enumerate(vcf_files):
        caller = _get_caller_name(vcf)
        meta, df = read_vcf(vcf)
        all_meta_headers.extend(meta)

        df = df.copy()
        # Remove intra-file duplicates before merging
        df = df.drop_duplicates(subset=MERGE_KEY)
        df['_caller'] = caller
        df['_priority'] = priority_idx
        tagged_frames.append(df)

    all_variants = pd.concat(tagged_frames, ignore_index=True)

    # For each unique variant select the best caller's row and record all callers.
    result_rows: list[pd.Series] = []
    for _key, group in all_variants.groupby(MERGE_KEY, sort=False):
        group_sorted = group.sort_values('_priority')

        callers_str = ','.join(group_sorted['_caller'].tolist())
        best_row = group_sorted.iloc[0].copy()

        # Append CALLERS tag to INFO so provenance is traceable.
        info_val = best_row.get('INFO', '.')
        if pd.isna(info_val) or info_val == '.':
            best_row['INFO'] = f'CALLERS={callers_str}'
        else:
            best_row['INFO'] = f'{info_val};CALLERS={callers_str}'

        result_rows.append(best_row)

    result_df = (
        pd.DataFrame(result_rows)
        .drop(columns=['_caller', '_priority'], errors='ignore')
    )

    # Sort output variants by chromosome then position.
    result_df['_pos_int'] = pd.to_numeric(result_df['POS'], errors='coerce')
    result_df['_chrom_sort'] = _chrom_sort_key(result_df['#CHROM'])
    result_df = (
        result_df
        .sort_values(['_chrom_sort', '_pos_int'])
        .drop(columns=['_pos_int', '_chrom_sort'], errors='ignore')
    )

    # Write valid VCF: ## headers first, then column header + data rows.
    merged_headers = _merge_meta_headers(all_meta_headers)
    with open(output_file, 'w') as out:
        for h in merged_headers:
            out.write(h + '\n')
        result_df.to_csv(out, sep='\t', index=False)

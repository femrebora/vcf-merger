"""Typed domain models for VCF harmonization."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Optional, Protocol, runtime_checkable


class AnalysisMode(str, Enum):
    GERMLINE = "germline"
    SOMATIC = "somatic"


class EnsembleStrategy(str, Enum):
    UNION = "union"
    PASS_UNION = "pass-union"
    CONSENSUS = "consensus"
    CALLER_SPECIFIC = "caller-specific"


class DetectionSource(str, Enum):
    CLI = "cli"
    HEADER = "header"
    FILENAME = "filename"
    UNKNOWN = "unknown"


class FileKind(str, Enum):
    VCF = "vcf"
    VCF_GZ = "vcf.gz"
    GVCF = "gvcf"
    GVCF_GZ = "gvcf.gz"
    BCF = "bcf"
    UNKNOWN = "unknown"


@runtime_checkable
class VrsIdentifier(Protocol):
    """Optional future GA4GH VRS identifier abstraction (not required offline)."""

    def identify(self, assembly: str, contig: str, pos: int, ref: str, alt: str) -> str:
        ...


@dataclass(frozen=True, slots=True)
class CanonicalVariant:
    """Normalized small-variant identity."""

    assembly: str
    contig: str
    pos: int
    ref: str
    alt: str
    vrs_id: Optional[str] = None

    @property
    def key(self) -> tuple[str, str, int, str, str]:
        return (self.assembly, self.contig, self.pos, self.ref, self.alt)

    def __str__(self) -> str:
        return f"{self.assembly}:{self.contig}:{self.pos}:{self.ref}>{self.alt}"


@dataclass(slots=True)
class SampleEvidence:
    """Per-sample technical evidence from one caller record."""

    sample: str
    genotype: Optional[str] = None
    phased: bool = False
    ploidy: Optional[int] = None
    depth: Optional[int] = None
    ref_depth: Optional[int] = None
    alt_depth: Optional[int] = None
    allele_fraction: Optional[float] = None
    genotype_quality: Optional[float] = None
    raw_format: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class CallerEvidence:
    """Technical evidence from one caller for one canonical allele."""

    caller: str
    caller_version: Optional[str]
    source_file: str
    filter_status: str
    passed: bool
    qual: Optional[float]
    samples: list[SampleEvidence] = field(default_factory=list)
    caller_specific_metrics: dict[str, Any] = field(default_factory=dict)
    raw_info: dict[str, Any] = field(default_factory=dict)
    original_chrom: str = ""
    original_pos: int = 0
    original_ref: str = ""
    original_alt: str = ""
    original_record_identifier: str = ""
    detection_source: DetectionSource = DetectionSource.UNKNOWN
    # Somatic-oriented optional fields
    tumor_allele_fraction: Optional[float] = None
    normal_allele_fraction: Optional[float] = None


@dataclass(slots=True)
class HarmonizedVariant:
    """Canonical allele plus aggregated caller evidence."""

    canonical: CanonicalVariant
    evidence: list[CallerEvidence]
    gt_conflict: bool = False
    representative_gt: Optional[str] = None
    unsupported: bool = False
    unsupported_reason: Optional[str] = None

    @property
    def callers_all(self) -> list[str]:
        return sorted({e.caller for e in self.evidence})

    @property
    def callers_pass(self) -> list[str]:
        return sorted({e.caller for e in self.evidence if e.passed})

    @property
    def caller_count(self) -> int:
        return len(self.callers_all)

    @property
    def pass_caller_count(self) -> int:
        return len(self.callers_pass)


@dataclass(slots=True)
class ContigInfo:
    name: str
    length: Optional[int] = None


@dataclass(slots=True)
class InputVCFMetadata:
    path: str
    file_kind: FileKind
    fileformat: Optional[str] = None
    assembly: Optional[str] = None
    reference: Optional[str] = None
    contigs: list[ContigInfo] = field(default_factory=list)
    samples: list[str] = field(default_factory=list)
    caller: Optional[str] = None
    caller_version: Optional[str] = None
    caller_detection_source: DetectionSource = DetectionSource.UNKNOWN
    is_gvcf: bool = False
    gvcf_reasons: list[str] = field(default_factory=list)
    record_count_estimate: Optional[int] = None
    header_info_ids: list[str] = field(default_factory=list)
    header_format_ids: list[str] = field(default_factory=list)
    header_filter_ids: list[str] = field(default_factory=list)


@dataclass(slots=True)
class ProvenanceRecord:
    tool_name: str
    tool_version: str
    timestamp: str
    command_line: list[str]
    analysis_mode: str
    ensemble_strategy: str
    reference_path: Optional[str]
    reference_checksum: Optional[str]
    normalization: dict[str, Any]
    inputs: list[dict[str, Any]]
    samples: list[str]
    notes: list[str] = field(default_factory=list)


@dataclass(slots=True)
class OriginalAllele:
    chrom: str
    pos: int
    ref: str
    alt: str
    source_caller: str
    source_file: str

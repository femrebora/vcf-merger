# Migration from v0.1 scripts

## Breaking changes

1. **Priority merge is no longer the default scientific model.**  
   Evidence from all callers is retained. Default strategy is `union`.

2. **Normalization + reference are required by default** for `merge` (`--reference`).  
   Use `--no-normalize` only for debugging.

3. **gVCFs are rejected** for ensemble merging.

4. **Batch “≥3 callers” gate removed** from the legacy batch wrapper.

5. **INFO tag** `CALLERS` is replaced by `VM_CALLERS` / `VM_PASS_CALLERS` (and related `VM_*` fields).

6. **Caller detection** prefers CLI → header → filename. Names such as `sample.FB.norm.vcf` no longer break detection.

7. **Sample sets must match** across inputs; silent pandas column alignment is gone.

## Legacy wrappers

- `Merge_All_VCFs.py` — edit paths, then run; emits `DeprecationWarning`.
- `Merge_all_VCF_Groups.py` — `python Merge_all_VCF_Groups.py IN_DIR OUT_DIR`.
- `vcf_utils.merge_vcfs(...)` — shim calling `harmonize_vcfs` with `strategy=union`.

Prefer:

```bash
vcf-merger merge --mode germline --reference REF.fa -i a.vcf -i b.vcf -o out.vcf.gz
```

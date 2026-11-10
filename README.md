# Mutational signatures from scratch

Find mutational processes in cancer genomes with nothing but numpy. Everything
here follows the five-part post "5 cancer genomics projects you can do at home
with Claude" (@thebiosciencedesk).

The pipeline:

1. Load 21 breast cancer whole genomes as a 96-channel mutation catalogue (183,916 substitutions).
2. Count the APOBEC-associated channels before fitting anything.
3. Factorise the 96 × 21 matrix into K signatures with Lee & Seung multiplicative updates (~10 lines).
4. Match the learned signatures against COSMIC v3.4 by cosine similarity.
5. Test how stable the fit is across random starts.

## Setup

```bash
pip install -r requirements.txt
python run_analysis.py
python tests/test_slides.py   # or: python -m pytest -q
```

Python 3.10+ and numpy are required. scipy is only needed to regenerate the data.

## Data

- `data/breast21_sbs96.txt`: 96-channel counts for 21 breast cancer genomes. Converted
  from `21_breast_WGS_substitutions.mat`, shipped with the SigProfilerAssignment and
  SigProfilerExtractor pip packages. Regenerate with `python prepare_data.py <mat> data/breast21_sbs96.txt`.
- `data/cosmic_v3.4_sbs_grch37.txt`: COSMIC v3.4 SBS reference signatures (GRCh37), 86 signatures,
  from the SigProfilerAssignment package.

## Files

- `signatures.py`: `tsv`, `cos`, `nmf`, `pair`, and the APOBEC channel count (close to the slide code).
- `run_analysis.py`: the five steps, with printed output.
- `prepare_data.py`: converts the MATLAB matrix to text.
- `tests/test_slides.py`: checks the numbers from the slides.

## Results

| Step | Slide | This code |
|---|---|---|
| Total mutations | 183,916 | 183,916 |
| Largest tumour | PD4120a, 70,690 (47.6× quietest) | same |
| TCA/TCT share of C>T/C>G | 60.2% vs 12.5% expected (4.8×) | same |
| K=3 reconstruction cosine | 0.99968 | 0.99968 |
| K=2 to K=6 fit | 0.99955 to 0.99989 | same |
| Signatures found | SBS2, SBS13 (APOBEC), SBS3 (HR repair / BRCA) | SBS2, SBS13, SBS3 |

## Known differences

- **Signature cosines.** The slide reports SBS13 at 0.847 and SBS3 at 0.900 (runner-up SBS40c at 0.845).
  This code gives SBS13 at 0.913 and SBS3 at 0.899 (runner-up SBS40c at 0.850). The same three signatures are
  recovered, but the exact cosine values do not match the slide.
- **Repeatability.** The slide reports 0.982 (K=3) and 0.860 / 0.748 / 0.918 (K=4 to K=6). This code gives
  0.959, 0.829, 0.806 and 0.872 for the same K values, using seeds 0 to 7. The fit values match exactly, which
  suggests the restarts converge to different local optima. The seeding scheme in the post may differ from this code.

These differences are reported as found, not adjusted to match the slide.

## Caveats (from the post)

- 21 tumours from one study, not a representative cohort.
- No opportunity correction: the genome does not contain the 96 contexts in equal numbers.
- A signature is a pattern, not a mechanism. Calling one "APOBEC" is an inference from cosine similarity.
- Flat signatures (such as SBS3 vs SBS40c) are hard to tell apart. This is an open problem in the field.

## Credit

Analysis outline and code from @thebiosciencedesk. Breast cancer genomes and the matrix are distributed with
SigProfiler. COSMIC signatures are from the COSMIC project (Sanger Institute).

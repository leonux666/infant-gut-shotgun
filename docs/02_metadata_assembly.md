# 02 — Metadata assembly

How to rebuild, from scratch, the joined metadata table that drives sample selection for this
project. Run on the Amarel login node; every step here is network or light I/O, no compute.

**Dataset**: Shao et al. 2019, *Nature* 574:117–121 — Baby Biome Study infant gut shotgun
metagenomes. ENA study accessions `ERP115334` and `ERP024601`.

**End state**: `data/metadata_joined.tsv` — one row per metagenomic run, carrying delivery mode,
sampling time point, 20 further clinical covariates, and FASTQ download paths.

Verified end to end on 2026-08-14.

---

## Prerequisites

```bash
conda activate nfcore
```

`pandas` is already present as an indirect dependency of `nf-core/tools`. `openpyxl` is not, and
`pandas.read_excel()` needs it as the parsing engine:

```bash
mamba install -y -c conda-forge openpyxl
```

`-y` skips the confirmation prompt. Not pinned: this is a one-off reader, not a pipeline
component. Installing packages on the login node is acceptable (I/O bound, not compute).

Confirm the install did not disturb the environment. `mamba` relinks a large number of packages
when it resolves, including `nf-core` itself, so a successful transaction is not by itself
evidence that the environment still works:

```bash
python -c "import pandas, openpyxl; print('pandas', pandas.__version__, '| openpyxl', openpyxl.__version__)" \
  && nextflow -version | head -4 \
  && nf-core --version
```

Expect `pandas 3.0.5 | openpyxl 3.1.5`, Nextflow `26.04.6`, nf-core `4.1.0`.

---

## Step 1 — Pull the ENA run report

The ENA Portal API `filereport` endpoint returns a flat table for a given accession. `result=`
selects the record level; `read_run` is the sequencing-run level, which is where the FASTQ files
hang. `fields=` selects columns, **returned in the order requested** — this matters, because every
`awk`/`cut` below indexes by position.

```bash
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=ERP115334&result=read_run&fields=run_accession,sample_accession,library_strategy,library_source,library_layout,instrument_model,read_count,base_count,fastq_bytes,fastq_ftp&format=tsv" -o ~/ena_ERP115334.tsv && wc -l ~/ena_ERP115334.tsv
```

Expect `2388` lines: 1 header + 2387 runs.

The `&& wc -l` is the verification half. `curl -s -o` writes silently and exits 0 on an empty
response, so counting lines checks the target state rather than the absence of an error.

---

## Step 2 — Establish what the study actually contains

Do not assume the study is homogeneous. Count the actual values of `library_source` (column 4) and
`library_layout` (column 5):

```bash
cut -f4,5 ~/ena_ERP115334.tsv | sort | uniq -c
```

```
 708 GENOMIC      PAIRED
   1 library_source library_layout   <- header row, counted as data
1679 METAGENOMIC  PAIRED
```

`sort` before `uniq -c` is mandatory: `uniq` only collapses *adjacent* duplicates.

**`ERP115334` is mixed.** 1679 faecal metagenomes plus 708 bacterial isolate genomes. The 1679
matches the sample count reported in the paper exactly — one run per sample, no resequencing.
Everything downstream must filter on `library_source == "METAGENOMIC"`.

All runs are paired-end.

The second accession is out of scope, and this is worth confirming rather than assuming:

```bash
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=ERP024601&result=read_run&fields=run_accession,sample_accession,library_strategy,library_source,library_layout,instrument_model,read_count,base_count,fastq_bytes,fastq_ftp&format=tsv" -o ~/ena_ERP024601.tsv && cut -f4,5 ~/ena_ERP024601.tsv | sort | uniq -c
```

805 runs, all `GENOMIC`. Bacterial isolates only. Not used by this project.

---

## Step 3 — Size the download

Subset to the metagenomes, keeping the header:

```bash
awk -F'\t' 'NR==1 || $4=="METAGENOMIC"' ~/ena_ERP115334.tsv > ~/ena_metagenomic.tsv && wc -l ~/ena_metagenomic.tsv
```

Expect `1680` (1679 + header) — a hard check on the filter.

`-F'\t'` is required. Awk's default separator is whitespace, and several fields (e.g.
`Illumina HiSeq 4000`) contain spaces, which would shift every column index. Matching on `$4`
rather than `grep`-ing the whole line keeps the test anchored to the intended column.

`fastq_bytes` (column 9) holds two semicolon-separated integers for paired data — one per mate
file. Split, sum, and report the distribution:

```bash
awk -F'\t' 'NR>1 {split($9,a,";"); s=a[1]+a[2]; total+=s; n++; if(s<min||n==1)min=s; if(s>max)max=s} END {printf "runs=%d  total=%.1f GB  mean=%.2f GB  min=%.2f GB  max=%.2f GB\n", n, total/1e9, total/n/1e9, min/1e9, max/1e9}' ~/ena_metagenomic.tsv
```

```
runs=1679  total=2042.6 GB  mean=1.22 GB  min=0.32 GB  max=4.45 GB
```

`runs=1679` is a self-check. `n==1` in the min test is necessary: an uninitialised awk variable
never compares as greater, so without the special case for the first row `min` stays empty.

Report min and max, not just the mean — the spread determines whether a mean-based estimate of
download volume is safe. Here max/mean ≈ 3.6, no extreme outliers.

**The full study is ~2.0 TB. Do not download it all.** ~24 samples is ~30 GB.

---

## Step 4 — Pull the join keys

ENA carries two sample identifiers: `sample_accession` (`SAMEA...`, primary/BioSample) and
`secondary_sample_accession` (`ERS...`). The published clinical table uses the secondary form, so
it has to be requested explicitly.

```bash
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=ERP115334&result=read_run&fields=run_accession,sample_accession,secondary_sample_accession,library_source,fastq_bytes,fastq_ftp&format=tsv" -o ~/ena_link.tsv && awk -F'\t' 'NR<=3 {print $1"|"$2"|"$3"|"$4}' ~/ena_link.tsv
```

```
run_accession|sample_accession|secondary_sample_accession|library_source
ERR13330869|SAMEA115771207|ERS20358348|GENOMIC
ERR13330875|SAMEA115771236|ERS20358377|GENOMIC
```

Piping through `|` makes empty fields visible; a blank column is invisible in raw TSV output.

Add `sample_title`, `read_count`, or `first_public` to `fields=` if needed later — but note that
the column indices used by `build_metadata.py` are resolved by name, not position, so extra fields
are safe.

---

## Step 5 — Download Supplementary Table 1

The clinical metadata is not in ENA (see Notes). It is published, deliberately de-identified, as
Supplementary Table 1 of the Nature paper. Supplementary files are outside the paywall.

```bash
curl -sL "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-019-1560-1/MediaObjects/41586_2019_1560_MOESM3_ESM.xlsx" -o ~/bbs_supp_table1.xlsx && ls -lh ~/bbs_supp_table1.xlsx && file ~/bbs_supp_table1.xlsx
```

Expect ~209K and `Microsoft Excel 2007+`.

Three things in this command are load-bearing:

- **`%3A` and `%2F` are URL-encoded `:` and `/`.** The DOI itself contains a slash; decoding them
  produces a 404.
- **`-L` follows redirects.** Publisher static hosts commonly 302 to a CDN. Without it, `curl`
  writes an empty file and still exits 0.
- **`file` reads the magic number, not the extension.** If the server returns an HTML error page,
  `-o` will happily save it as `.xlsx` and `ls` will show a plausible size. Only `file` catches it.

---

## Step 6 — Locate the header row

Journal supplementary tables put a caption in row 1. Reading with pandas defaults yields column
names like `Unnamed: 1 … Unnamed: 22` — that is the signal, not an absence of headers.

Read with `header=None` so nothing is interpreted as a header, and inspect the raw row structure:

```bash
python -c "
import pandas as pd
df = pd.read_excel('/home/xw347/bbs_supp_table1.xlsx', header=None, nrows=5)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 250)
print(df.iloc[:, :8])
"
```

Row 0 is the caption, row 1 is blank, **row 2 holds the column names**, data starts at row 3.
Hence `header=2` (0-indexed) in the script.

Single quotes inside, double quotes outside — the reverse closes the shell string early.

---

## Step 7 — Join

`src/build_metadata.py` filters the ENA report to metagenomes, reads the supplementary table at
the correct header offset, and left-joins on `secondary_sample_accession` ↔ `Accession`.

```bash
python ~/infant-gut-shotgun/src/build_metadata.py \
  --ena ~/ena_link.tsv \
  --supp ~/bbs_supp_table1.xlsx \
  --out ~/infant-gut-shotgun/data/metadata_joined.tsv
```

```
ENA report: 2387 runs, 1679 metagenomic
Supplementary table: 1679 rows, 23 columns
Columns: ['Accession', 'Individual', 'Time_point', 'Delivery_mode', 'C_Section_type', ...]
Join: 1679/1679 runs matched to clinical metadata
Wrote /home/xw347/infant-gut-shotgun/data/metadata_joined.tsv (1679 rows)
```

Every stage prints a countable number on purpose. Checks, in order:

| Number | Expected | If wrong |
|---|---|---|
| `2387 runs, 1679 metagenomic` | exact | the source filter is broken |
| Supplementary rows | 1679 | the table does not cover the full cohort |
| First column name | `Accession` | `SUPP_HEADER_ROW` is wrong |
| `Join: N/1679` | 1679 | see below |

Design decisions inside the script that are not optional:

- **`dtype=str` everywhere.** Otherwise pandas infers types, turning `Time_point` into an integer
  and stripping leading zeros from identifiers. Type inference is the most common source of silent
  bugs in metadata integration.
- **`validate="one_to_one"`.** Row duplication is the nastiest merge failure mode: no error, plausible
  output, wrong downstream. This raises instead.
- **`how="left"`.** Keeps all 1679 ENA rows so the match count is meaningful. An inner join would
  silently drop unmatched rows and hide the loss.
- **`SUPP_HEADER_ROW = 2` as a commented constant.** The value was determined empirically in step 6;
  the reason belongs in the code.

The output is `.gitignore`d. The script is what goes in the repository — anyone can rebuild the
table from it.

---

## Output columns

`data/metadata_joined.tsv` carries the ENA fields plus all 23 columns of Supplementary Table 1:

```
Accession, Individual, Time_point, Delivery_mode, C_Section_type,
Infancy_sampling_age_months, Postnatal_stay_in_hospital_days, Hospital,
Hospital_destination_after_birth, Gender, Mother_age, Birth_weight,
Feeding_method, Breastfeeding_status, Breastfeeding_1hr_birth,
Abx_mother_prior_birth, Abx_mother_labour_IAP, Abx_mother_after_hospital,
Abx_Baby_in_hospital, Abx_Baby_after_hospital, Bacteroides_profile,
WGS_reads_raw, WGS_reads_trimmed
```

Enough to stratify sample selection by delivery mode and time point while balancing antibiotic
exposure. `WGS_reads_raw` vs `WGS_reads_trimmed` quantifies what quality trimming and host removal
took out, which is the empirical test of whether the deposited reads are already decontaminated.

---

## Notes

Findings that shape the steps above, recorded so they are not re-derived.

**Clinical metadata is absent from ENA by design.** The samples are registered against ENA
checklist `ERC000011`, the permissive generic checklist, and carry a single non-system attribute,
`SUBJECT_ID` (e.g. `5826STDY7976012`). No delivery mode, no collection date, no host attributes.
ENA is an open-access sequence archive; individual-level clinical variables for a UK birth cohort
are distributed separately, de-identified, through the publication. This is data governance, not
an access restriction — the 2024 Nature Microbiology paper from the same cohort states plainly
that participant-level clinical metadata is provided in its Supplementary Tables.

**`result=sample` rejects study accessions.** The Portal API returns
`Accession(s) ERP115334 not valid for search requests on sample data`, listing `SAME…` / `ERS…` as
the accepted forms. Sample-level reports cannot be expanded from a study accession the way run
reports can. To inspect raw sample attributes, take a real sample accession from the run report
and use the Browser API instead: `curl -s "https://www.ebi.ac.uk/ena/browser/api/xml/SAMEA5616559"`.
The Portal API only returns its own standardised fields, so submitter-defined attributes are
invisible there.

**`SUBJECT_ID` is a dead end.** It is a Sanger-internal identifier and does not correspond to the
`Individual` column of the supplementary table. Do not try to map it.

**`sample_title` encodes individual and day of life** (e.g. `513122_4`, `513122_6` — same infant,
two time points), matching the `Individual` numbering of the supplementary table. It is not needed,
since `Time_point` is an explicit column, but it is a useful independent cross-check. Note that
`library_name` is populated only for the 708 isolate runs, not for the metagenomes.

**Two submission batches share the study accession.** `first_public` splits `ERP115334` into
exactly 1679 (2019-07-01) and 708 (2024-06-27), the same partition as `library_source`. Two
independent fields agreeing is what makes the split trustworthy; a `head` of the file shows only
the isolate batch and invites the wrong conclusion.


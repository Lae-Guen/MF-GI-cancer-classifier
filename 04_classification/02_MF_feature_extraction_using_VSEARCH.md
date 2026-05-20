# MF Feature Extraction Using VSEARCH

## Purpose
Extract ASVs from external cohort 16S rRNA datasets that match the
Yonsei MF-transmitted (oral-to-gut) ASVs at varying sequence identity
thresholds (100%–70%). This enables construction of MF-feature-based
external validation sets without requiring paired oral-fecal samples.

All comparisons are performed in V3-V4 (or V3-V4-aligned V4) amplicon space.

## Workflow Overview

1. Extract the transmitted ASV subset from the Yonsei query representative sequences
2. Inspect query sequence length distribution (expected: 350–450 bp for V3-V4)
3. Inspect each external cohort's representative sequence lengths
4. Export query and cohort sequences to FASTA
5. Run VSEARCH against each cohort at identity thresholds 100%–70%
6. Extract matched external ASV IDs and query-to-target mapping pairs
7. Generate a hit summary table across cohorts and identity thresholds

---

## Step 1. Extract Transmitted ASV Subset from Query Representative Sequences

```bash
# Add QIIME2-required header to the transmitted ASV ID list
(echo "FeatureID"; cat transmittedASV.txt) > transmittedASV_metadata.tsv

# Filter the query representative sequences to transmitted ASVs only
qiime feature-table filter-seqs \
  --i-data query_rep-seqs.qza \
  --m-metadata-file transmittedASV_metadata.tsv \
  --o-filtered-data query_transmitted_rep-seqs.qza
```

---

## Step 2. Inspect Query V3-V4 Representative Sequences

```bash
qiime feature-table tabulate-seqs \
  --i-data query_transmitted_rep-seqs.qza \
  --o-visualization query_transmitted_rep-seqs.qzv
```

**Check:**
- Total number of transmitted ASVs
- Sequence length distribution
- Mean length within V3-V4 range (350–450 bp)

---

## Step 3. Inspect External Cohort Representative Sequences

```bash
for COHORT in russo flemer wang; do
  qiime feature-table tabulate-seqs \
    --i-data ${COHORT}_rep-seqs.qza \
    --o-visualization ${COHORT}_rep-seqs.qzv
done
```

**Acceptance criteria:**
- Mean length 350–450 bp → same V3-V4 region; direct matching is valid
- Cohorts with substantially shorter sequences (< 250 bp, V4-only) require
  separate handling (see the V4 analysis pipeline for Zhang, Baxter, Zackular, Zeller)

---

## Step 4. Export QZA Files to FASTA

```bash
mkdir -p v3v4_identity/{exports,hits,asv_lists,logs}

for COHORT in query russo flemer wang; do
  qiime tools export \
    --input-path ${COHORT}_rep-seqs.qza \
    --output-path v3v4_identity/exports/${COHORT}
done
```

**Output:**
- `v3v4_identity/exports/<cohort>/dna-sequences.fasta`

---

## Step 5. Run VSEARCH at Identity Thresholds 100%–70%

All three cohorts are processed in parallel (background jobs).

```bash
#!/usr/bin/env bash
set -euo pipefail

QUERY="v3v4_identity/exports/query/dna-sequences.fasta"
PIDS="1.00 0.99 0.98 0.96 0.94 0.90 0.85 0.80 0.70"
THREADS=20

for COHORT in russo flemer wang; do
  (
    echo "[START] ${COHORT} analysis..."
    mkdir -p v3v4_identity/hits/${COHORT} v3v4_identity/asv_lists/${COHORT}
    for pid in $PIDS; do
      vsearch \
        --usearch_global "$QUERY" \
        --db v3v4_identity/exports/${COHORT}/dna-sequences.fasta \
        --id "$pid" \
        --iddef 2 \
        --strand both \
        --threads "$THREADS" \
        --maxaccepts 0 \
        --maxrejects 0 \
        --blast6out v3v4_identity/hits/${COHORT}/${COHORT}_${pid}.tsv
    done
    echo "[DONE] ${COHORT} analysis finished."
  ) &
done

wait
echo "[ALL DONE] All cohorts processed."
```

**VSEARCH parameters:**
- `--iddef 2` : Edit-distance-based identity (recommended for amplicons)
- `--strand both` : Match both forward and reverse complement
- `--maxaccepts 0 / --maxrejects 0` : Report all hits above threshold

---

## Step 6. Extract Matched External ASV IDs and Query-to-Target Pairs

```bash
for COHORT in russo flemer wang; do
  for pid in 1.00 0.99 0.98 0.96 0.94 0.90 0.85 0.80 0.70; do

    # Column 2 = matched external (database) ASV ID
    cut -f2 v3v4_identity/hits/${COHORT}/${COHORT}_${pid}.tsv | sort -u \
      > v3v4_identity/asv_lists/${COHORT}/${COHORT}_ASV_${pid}.txt

    # Columns 1-2 = query ASV → matched external ASV pairs
    cut -f1,2 v3v4_identity/hits/${COHORT}/${COHORT}_${pid}.tsv | sort -u \
      > v3v4_identity/asv_lists/${COHORT}/query_to_${COHORT}_${pid}.tsv

    # Add QIIME2-compatible FeatureID header for downstream filtering
    { echo "FeatureID"
      cat v3v4_identity/asv_lists/${COHORT}/${COHORT}_ASV_${pid}.txt
    } > v3v4_identity/asv_lists/${COHORT}/${COHORT}_ASV_${pid}.tsv

  done
done
```

---

## Step 7. Generate Hit Summary Table

```bash
echo -e "Cohort\tIdentity\tMatched_external_ASVs\tQuery_to_target_pairs" \
  > v3v4_identity/identity_hit_summary_v3v4.tsv

for COHORT in russo flemer wang; do
  for pid in 1.00 0.99 0.98 0.96 0.94 0.90 0.85 0.80 0.70; do
    asv_txt="v3v4_identity/asv_lists/${COHORT}/${COHORT}_ASV_${pid}.txt"
    pair_tsv="v3v4_identity/asv_lists/${COHORT}/query_to_${COHORT}_${pid}.tsv"

    n_asv=0; n_pair=0
    [[ -s "$asv_txt"  ]] && n_asv=$(wc -l  < "$asv_txt")
    [[ -s "$pair_tsv" ]] && n_pair=$(wc -l < "$pair_tsv")

    echo -e "${COHORT}\t${pid}\t${n_asv}\t${n_pair}" \
      >> v3v4_identity/identity_hit_summary_v3v4.tsv
  done
done
```

**Output:** `v3v4_identity/identity_hit_summary_v3v4.tsv`
→ Used as input to `03_sequence_similarity_based_mf_feature.R`

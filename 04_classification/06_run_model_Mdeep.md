# Running MDeep

## Overview
MDeep is a phylogeny-aware deep learning classifier for microbiome data.
It uses a convolutional neural network that incorporates a phylogenetic
correlation matrix C to capture evolutionary relationships between ASVs.

- Repository: https://github.com/alfredyewang/MDeep
- Required inputs (per seed): X_train.npy, Y_train.npy, X_eval.npy, Y_eval.npy, c.npy
  → Produced by `05_model_mdeep_inputfile.R`

---

## 1. Environment Installation

TensorFlow 1.12.0 is required. Use the `anaconda` channel to obtain it.

```bash
mamba create -n mdeep \
  -c anaconda -c conda-forge \
  python=3.6 \
  tensorflow=1.12.0 \
  numpy=1.16 \
  scipy=1.2 \
  pandas=0.24 \
  scikit-learn=0.20 \
  matplotlib=3.0 \
  seaborn \
  h5py=2.9 \
  -y

conda activate mdeep
```

## 2. Clone Repository and Verify Installation

```bash
git clone https://github.com/alfredyewang/MDeep.git
cd MDeep
pip install -r requirements.txt
python src/MDeep.py -h
```

## 3. Fix Matplotlib Backend (Headless Server)

Add the following two lines to the top of `src/binary.py` to disable
interactive display (required for server environments without a display):

```python
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
```

---

## 4. Key Hyperparameters

### Data and mode
| Argument | Description |
|----------|-------------|
| `--train` / `--evaluation` / `--test` | Mode: train model, evaluate on held-out data, or predict new samples |
| `--data_dir` | Directory containing .npy input files for one iteration |
| `--outcome_type binary` | Binary classification (HC vs disease) |

### Model architecture
| Argument | Value used | Description |
|----------|------------|-------------|
| `--window_size` | 8 8 8 | Convolutional filter width per layer (number of ASVs grouped together) |
| `--kernel_size` | 64 64 32 | Number of filters per convolutional layer |
| `--strides` | 1 1 1 | Stride of convolution (1 = slide one ASV at a time) |

### Training
| Argument | Value used | Description |
|----------|------------|-------------|
| `--max_epoch` | 2000 | Training epochs |
| `--batch_size` | 32 | Mini-batch size |
| `--learning_rate` | 0.001 | Adam optimizer learning rate |
| `--dropout_rate` | 0.6 | Fraction of connections randomly dropped each step (prevents overfitting) |
| `--L2_regularizer` | 0.01 | L2 weight penalty (reduces large weights) |

### Phylogenetic correlation (ρ)
The correlation matrix is `C = exp(-2ρD)` where D is the cophenetic distance.
- ρ = 2 (used in this study): moderate decay; nearby taxa are strongly correlated
- Smaller ρ → more taxa treated as correlated; larger ρ → only close relatives linked

---

## 5. Automated Train + Evaluation Script (50 Seeds)

```bash
#!/bin/bash
set -euo pipefail

BASE_DATA="data/GC_mouth_mtg"
BASE_MODEL="model/GC_mouth_mtg"
BASE_RESULT="result/GC_mouth_mtg"
SUMMARY_FILE="summary_results_GC_mouth_mtg.csv"

echo "seed,auc" > "$SUMMARY_FILE"

for seed in $(seq 100 149); do
    echo "--- Seed ${seed} ---"

    mkdir -p "${BASE_MODEL}/iter_${seed}"
    mkdir -p "${BASE_RESULT}/iter_${seed}"

    # Train
    python3 src/MDeep.py --train \
        --data_dir  "${BASE_DATA}/iter_${seed}" \
        --model_dir "${BASE_MODEL}/iter_${seed}" \
        --outcome_type binary \
        --batch_size 32 --max_epoch 2000 \
        --learning_rate 0.001 --dropout_rate 0.6 \
        --window_size 8 8 8 --kernel_size 64 64 32 --strides 1 1 1 \
        --L2_regularizer 0.01

    # Evaluate and capture log
    EVAL_LOG="${BASE_RESULT}/iter_${seed}/eval_log.txt"
    python3 src/MDeep.py --evaluation \
        --data_dir   "${BASE_DATA}/iter_${seed}" \
        --model_dir  "${BASE_MODEL}/iter_${seed}" \
        --result_dir "${BASE_RESULT}/iter_${seed}" \
        --outcome_type binary \
        --window_size 8 8 8 --kernel_size 64 64 32 --strides 1 1 1 \
        > "$EVAL_LOG" 2>&1

    # Extract AUC value from last line of evaluation log
    auc=$(tail -n 1 "$EVAL_LOG" | tr -d '[:space:]' | grep -oE '[0-9]+\.[0-9]+')
    if [ -z "$auc" ]; then auc="NaN"; fi

    echo "${seed},${auc}" >> "$SUMMARY_FILE"
    echo "Seed ${seed} complete — AUC: ${auc}"
done

echo "All seeds complete. Results in: ${SUMMARY_FILE}"
```

**Output:** `summary_results_GC_mouth_mtg.csv` — seed × AUC table,
directly comparable to the RF/glmnet results from `10_run_internal_rf_optimized.R`.

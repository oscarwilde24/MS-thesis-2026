#  HLA Allele Frequency Analysis (Global vs India)

This project analyzes and visualizes HLA Class I allele frequencies (A, B, C) by comparing global populations with Indian populations.

It generates a merged bar plot showing the top alleles across both groups.

## Requirements

Install the following Python packages:

```bash
pip install pandas numpy matplotlib seaborn tqdm
```

### Suggested versions

```text
Python >= 3.9
pandas >= 2.0
numpy >= 1.24
matplotlib >= 3.7
seaborn >= 0.12
tqdm >= 4.65
```

---

## Required Files

You need three CSV files:

```text
HLA_A_AlleleFrequencies_AllPages_FIXED.csv
HLA_B_AlleleFrequencies_AllPages_FIXED.csv
HLA_C_AlleleFrequencies_AllPages_FIXED.csv
```

### Expected folder structure

```text
data/
├── A/
│   └── HLA_A_AlleleFrequencies_AllPages_FIXED.csv
├── B/
│   └── HLA_B_AlleleFrequencies_AllPages_FIXED.csv
└── C/
    └── HLA_C_AlleleFrequencies_AllPages_FIXED.csv
```

---

## Expected CSV format

Each file should contain:

```text
Allele | Region | Freq | Sample Size
```

---

## How to Run

1. Update the base path in the script:

```python
BASE_DIR = Path("data")
```

2. Run:

```bash
python script_name.py
```

---

## Output

The script generates:

```text
HLA_Merged_Global_India.png
```

This contains:

* Top 10 global alleles
* Top 10 Indian alleles
* For HLA-A, HLA-B, and HLA-C

---

## Notes

* Frequencies are converted into allele counts using:

  ```
  freq × 2 × sample size
  ```
* Indian populations are identified using:

  ```
  Region starts with "India"
  ```
* Only two-field alleles (e.g., A*02:01) are used

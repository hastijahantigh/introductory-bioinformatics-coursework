# CodonBERT Gene Expression Prediction

Final project for the Introduction to Bioinformatics course, Fall 2025.

**Students:** Hasti Jahantigh and Bahar Marashi

## Overview

This project reviews and reproduces the main workflow from the paper [CodonBERT: Using BERT for Sentiment Analysis to Better Predict Genes with Low Expression](https://doi.org/10.1145/3584371.3613013).

The main goal is to use codon sequences to predict gene expression. We reproduced the workflow for *Saccharomyces cerevisiae* using an existing pretrained CodonBERT checkpoint.

The project contains two stages:

1. Classifying genes into low, medium, and high expression groups
2. Predicting a numerical gene expression value using regression

## Adaptation from the Original Paper

The original paper trained its BERT model from scratch. This was not practical with the computational resources available to us because it required long and stable GPU access.

Instead, we used the public [`lhallee/CodonBERT`](https://huggingface.co/lhallee/CodonBERT) checkpoint. These weights come from a separate CodonBERT model pretrained on a large collection of mRNA sequences.

Therefore, this project is an adapted reproduction of the paper's downstream classification and regression workflow, not an exact reproduction of its original pretraining process.

## Dataset

We used data for baker's yeast, *Saccharomyces cerevisiae*.

The project combines:

- Gene expression counts from [Expression Atlas experiment E-MTAB-8626](https://www.ebi.ac.uk/gxa/experiments/E-MTAB-8626)
- Coding sequences from the *S. cerevisiae* S288C reference genome, assembly R64-1-1

The expression value for each gene was calculated as the median count across the available samples.

After cleaning and merging the sequence and expression data, the final dataset contained `6,030` genes.

## Workflow

### 1. Data preparation

- Download the expression and CDS files
- Clean the coding sequences
- Convert DNA sequences to RNA sequences
- Match sequences and expression values using gene identifiers
- Apply `log1p` to the expression values

### 2. Classification labels

The genes were divided into three expression groups using the 33rd and 67th percentiles:

- Class 0: low expression
- Class 1: medium expression
- Class 2: high expression

### 3. Tokenization

Each RNA sequence was divided into codons containing three nucleotides. Spaces were added between codons before tokenization.

The final model input was limited to `1,024` tokens.

### 4. Classification

A linear classification layer was added to the pretrained CodonBERT backbone. The model was trained to predict the three expression classes.

### 5. Regression

The classification head was replaced with a regression head. The first six transformer layers were frozen, and the remaining layers were fine-tuned to predict numerical expression values.

## Results

### Classification

| Metric | Value |
|---|---:|
| Accuracy | 0.6652 |
| Macro precision | 0.6699 |
| Macro recall | 0.6657 |
| Macro F1-score | 0.6673 |
| Weighted F1-score | 0.6665 |
| AUROC | 0.8311 |

### Regression

| Metric | Value |
|---|---:|
| Scaled MSE | 0.0284 |
| Scaled MAE | 0.1142 |
| Spearman correlation | 0.8632 |

## Repository Files

```text
final-project-codonbert/
├── README.md
├── codonbert.ipynb
└── codonbert_reproduction_report.pdf
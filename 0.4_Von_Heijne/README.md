# LB2_project_Group_7

## Data collection
Retrieved data from UniProt with a specific query (first filtering) and then polished the .tsv datasets from entries whose cleavage was absent or unknown (second filtering). The final collection of data consisted of these numbers of entries:<br>

| Dataset  | Entries    | %      |
|----------|-----------:|-------:|
| (+) Positive | **2,956**  | 12.4%  |
| (-) Negative | **20,974** | 87.6%  |



## Data preparation
Results in .fasta format were aligned and clustered with mmseqs, and one representative of each cluster was copied in a new .tsv file and labelled with "0" (negatives) or "1" (positives).<br>
The representatives were shuffled and then split between the training (80% of positives and 80% of negatives) and the benchmark dataset (remaining 20% of both).<br>

The previous labelling was useful to split the training dataset in 5 subsets, each with 20% of both positives and negatives (in this way, the ratio will be constant among all the subsets).<br>

### Training sets

| Subset    | (+) Positives | (-) Negatives | Total      | % of dataset |
|-----------|----------:|----------:|-----------:|-------------:|
| sub_1     | 177       | 1454      | 1631       | 16.0%        |
| sub_2     | 176       | 1453      | 1629       | 16.0%        |
| sub_3     | 177       | 1453      | 1630       | 16.0%        |
| sub_4     | 176       | 1453      | 1629       | 16.0%        |
| sub_5     | 176       | 1453      | 1629       | 16.0%        |
| bench     | 220       | 1816      | 2036       | 20.0%        |
| **Total** | **1102**  | **9082**  | **10184**  |              |



## Data visualization

### Protein length distribution
| Training | Benchmark |
|---|---|
| ![](0.3_data_analysis/plot1_protein_length_distribution_training.png) | ![](0.3_data_analysis/plot1_protein_length_distribution_benchmarking.png) |

### Signal peptide length distribution ~ positives
| Training | Benchmark |
|---|---|
| ![](0.3_data_analysis/plot2_SP_length_distribution_training.png) | ![](0.3_data_analysis/plot2_SP_length_distribution_benchmarking.png) |

### Amino-acid composition (SP vs Background)
| Training | Benchmark |
|---|---|
| ![](0.3_data_analysis/plot3_aa_composition_training.png) | ![](0.3_data_analysis/plot3_aa_composition_benchmarking.png) |

### Taxonomic classification ~ kingdom
| Training (positives) | Training (negatives) |
|---|---|
| ![](0.3_data_analysis/plot4a_kingdoms_training_pos.png) | ![](0.3_data_analysis/plot4a_kingdoms_training_neg.png) |

| Benchmark (positives) | Benchmark (negatives) |
|---|---|
| ![](0.3_data_analysis/plot4c_kingdoms_benchmarking_pos.png) | ![](0.3_data_analysis/plot4c_kingdoms_benchmarking_neg.png) |

### Taxonomic classification ~ species
| Training (positives) | Training (negatives) |
|---|---|
| ![](0.3_data_analysis/plot4b_species_training_pos.png) | ![](0.3_data_analysis/plot4b_species_training_neg.png) |

| Benchmark (positives) | Benchmark (negatives) |
|---|---|
| ![](0.3_data_analysis/plot4d_species_benchmarking_pos.png) | ![](0.3_data_analysis/plot4d_species_benchmarking_neg.png) |

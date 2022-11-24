# Limestone Lab
Homework for GEGN 508 (Advanced Rock Mechanics), flat-roofed excavation assignment. 
To recreate the report both [conda](https://docs.conda.io/en/latest/) and 
[quarto](https://quarto.org/) are required. 

Once installed, the conda environment can be created:

```bash
conda env create -f environment.yml
conda activate flat
```

Then the report pdf is generated with quarto via

```bash
quarto render report.qmd
```

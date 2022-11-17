# ENCODE Single-Cell Working Group: Analyses for single-nucleus ATAC-seq data
> snATAC-seq processing and manual annotation pipeline for part of ENCODE single-cell data collection
> Please see the Experiment-ID.txt for more information.

Please contact Yu Fu (yu.fu2@umassmed.edu) if you have any question!

# Environment setup
- install conda
- install python
- ```conda env create -f environment.yml```

# Steps
## Step 0: Download the fragment file from ENCODE with authetication
```
python downloadData.py $EXP $File $User $pswd
tar -xvf $EXP.tar.gz
```
**Fragment File path: encode_scatac_dcc_2/results/[EXP]-1/fragments/fragments.tsv.gz )**

## Step 1: Build ArchR Project
```
Rscript run-ArchR.R $file $EXP $biosampleID
mkdir $EXP/MarkerFeatures
```

## Step 2: Filtering out subset of all-tissue marker gene list
**Can add "PanglaoDB" for more canonical marker genes)**
```
awk -v tissue=$tissue '{if (($2==tissue) || ($2=="multiple")) print $1}' > markers
```

**two ways of cell type annotation**
## Step 3A: Directly overlap marker genes list with marker features passing a certain FDR and Log2FC threshold
```
Rscript Pull-Marker-Gene.R $EXP/Save-ArchR-Project.rds $EXP/MarkerFeatures
./predict-cell-type-normax.sh $EXP/MarkerFeatures/ markers > overlap-results
```

## Step 3B: Visulize the raw or transformed gene score for annotation
```
Rscript Draw-UMAP.R $EXP/Save-ArchR-Project.rds $EXP
Rscript Draw-Marker-Heatmap.R $EXP/Save-ArchR-Project.rds markers $EXP
Rscript Draw-Dot-Plot $Exp/Save-ArchR-Project.rds markers $EXP
```

## Step 4: Assign labels
**Saved in proj$Labels and write to master table**
```
Rscript Draw-UMAP-with-Label.R $EXP/Save-ArchR-Project.rds $EXP
```


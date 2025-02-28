# Validation of calibration code in O2 for ALICE ITS2
Set of scripts and macros to validate the O2 calibration code for ITS2. These tools are used when a major code change is performed in the ITS2 calibration workflow. They allow to run the new code on all ITS2 calibration scans and then to compare the newly obtained ROOT TTrees with reference ones. 

# Where to get reference data 
In EOS, this folder `/eos/project/a/alice-its-commissioning/Threshold_scans_Run3/validation` contains raw TF data for each calibration scan. In addition, in the folders names `unknown_<run_number>` there are the ROOT TTtrees used as reference. Mount EOS on your computer or copy this data locally. 

# Run processing scripts
Imagining you have a new piece of code developed inside the calibration workflow scripts, these scripts allows you to run the new workflow on the set of raw TF downloaded in the previous step (reference raw TF). 
Inside the **O2 environment run**:

```
./run-full-validation.sh <path-where-the-reference-raw-TF-are>
```

This will run the ITS2 calibration workflow on all calibration scan data. It will take approximately 35 minutes. Once this is done, in the path where you run this script, you will get a new set of folders with names `unknown_<run_number>`. Inside each folder there are the ROOT TTtrees that will be used for the comparison with the reference TTtrees downloaded from EOS. 

# Run validation macro
The `validation_analysis.C` is a ROOT macro. Run it with:

```
root[0]: .L validation_analysis.C++
root[1]: validation_analysis("<path-where-the-reference-TTrees-are>")
```

Note that the new TTrees folders are supposed to be in the same path where you run the `run-full-validation.sh` script. 
This macro takes by default 1000 random entries from each tree and then compare the content of each branch. Once the analysis is completed (takes a few mins), a file `validation.pdf` is created with all the ratios between the data in the branches (new data vs reference data). If the ratio is 1 for all plots (red line in each plot is well visible), it means that your new code did not introduce any bug. 


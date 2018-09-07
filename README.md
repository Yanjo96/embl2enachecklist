*embl2enachecklist*
===================
Converts EMBL flatfiles to submission checklists (i.e., tab-separated spreadsheets) for submission to [ENA](http://www.ebi.ac.uk/ena) via [Webin](https://www.ebi.ac.uk/ena/submit/sra/#home).


INPUT
-----
* EMBL flatfile


FEATURES
-------------
* Checks if the type of DNA marker specified by the user is indeed present in the embl-input file (specifically in as qualifier value for a qualifier named "gene", "note" or "standard_name")



GENERAL USAGE
-------------

###### Checklist 'trnK_matK'
```
python2 scripts/embl2enachecklists_CMD.py
-e examples/input/matK.embl
-o examples/output/matK_SubmissionChecklist.tsv
-c trnK_matK
```

###### Checklist 'IGS'
```
python2 scripts/embl2enachecklists_CMD.py
-e examples/input/trnL-trnF.embl
-o examples/output/trnL-trnF_SubmissionChecklist.tsv
-c IGS
```


PREREQUISITES
-------------
* Input files must have the name of the DNA marker (e.g., "matK", "ITS") as qualifier value for a qualifier named "gene", "note" or "standard_name"



TO DO
-----

###### 1. An error in processing a sequence should break only the iteration of the loop, not the entire code execution.

###### 2. add ETS feature

###### 3. Write untitests for the functions in `ChecklistOps.py`

###### 4. Have the code automatically add non-mandatory qualifiers as separate columns
* Ensure that all features that are not mandatory are added as separate columns into the checklist output (and not dropped, as they are now)

###### 5. Die outlists in globale variablen speichern


CHANGELOG
---------
###### Version 0.0.2 (2018.05.22)
* Generated separate function to extract charset symbols
* Updated README
###### Version 0.0.1 (2018.05.16)
* Added example input and example output
* Added setup.py

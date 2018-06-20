*embl2enachecklist*
===================
Converts EMBL flatfiles to submission checklists (i.e., tab-separated spreadsheets) for submission to [ENA](http://www.ebi.ac.uk/ena) via [Webin](https://www.ebi.ac.uk/ena/submit/sra/#home).


INPUT
-----
* EMBL flatfile


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


TO DO
-----

###### 1. An error in processing a sequence should break only the iteration of the loop, not the entire code execution.

###### 2. Not all qualifiers of a gene have the name `'gene'`. Sometimes they are named `'note'` or `'standard_name'`. Adjust code to allow this.

###### 3. Write untitests for the functions in `ChecklistOps.py`

###### 4. Have the code automatically add non-mandatory qualifiers as separate columns
* Ensure that all features that are not mandatory are added as separate columns into the checklist output (and not dropped, as they are now)

###### 5. Write a GUI interface for input
* The GUI should consist of just one Window, where all functions are immediately visible; the GUI should not have any dropdown-menus. In general, the simpler the interface, the better.


CHANGELOG
---------
###### Version 0.0.2 (2018.05.22)
* Generated separate function to extract charset symbols
* Updated README
###### Version 0.0.1 (2018.05.16)
* Added example input and example output
* Added setup.py

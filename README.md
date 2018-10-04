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

###### GUI
```
python2 scripts/embl2enachecklists_GUI.py
```

###### Checklist 'trnK_matK'
```
python2 scripts/embl2enachecklists_CMD.py
-i examples/input/ETS.embl
-o examples/output/ETS.tsv
-c trnK_matK
-e no
```

###### Checklist 'IGS'
```
python2 scripts/embl2enachecklists_CMD.py
-i examples/input/trnL-trnF.embl
-o examples/output/trnL-trnF.tsv
-c IGS
-e no
```


PREREQUISITES
-------------
* Input files must have the name of the DNA marker (e.g., "matK", "ITS") as qualifier value for a qualifier named "gene", "note" or "standard_name"



CHANGELOG
---------
###### Version 0.1.0 (2018.10.04)
* dynamically output
* add checklist type ETS and gene_intron
* add argument environmental to commandline
* add GUI
* add error messages
* add warning messages
* add support for genbank files
###### Version 0.0.2 (2018.05.22)
* Generated separate function to extract charset symbols
* Updated README
###### Version 0.0.1 (2018.05.16)
* Added example input and example output
* Added setup.py

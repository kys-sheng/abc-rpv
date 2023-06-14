# `abc-rpv` - the RPV Python Library
`abc-rpv` is a Python library that provides a framework for analyzing the collider signatures of RPV-MSSM. Users are provided with various functionalities to explore the landscape of RPV-MSSM physics within the context of small RPV-couplings. 

Feel free to ask questions, request features and report bugs via the [Github's issue system](https://github.com/kys-sheng/abc-rpv/issues).  <br>
Suggestions and feedbacks are very much welcomed!

**!! [`Tutorial.ipynb`](https://github.com/kys-sheng/abc-rpv/blob/main/Tutorial.ipynb) contains everything relevant on how to use, more than the appendix in the original paper.** <br>
**A full manual will be available in the near future.**


## Installation / Download :
A. Download: Go to green color button (Code) and download as a zip   <br>
B. Installation: `git clone https://github.com/kys-sheng/abc-rpv.git`   <br>
C. pip3: In future updates

## Requirements
Used Python libraries:<br>
- itertools: `  pip3 install itertools or pip install itertools  `<br>
- numpy: `  pip3 install numpy or pip install numpy  `<br>
- os: `  pip3 install os or pip install os  `<br>
- pandas: `  pip3 install pandas or pip install pandas  `<br>
- subprocess: `  pip3 install subprocess or pip install subprocess  `<br>
- warnings: `  pip3 install warnings or pip install warnings  `<br>

## Manual
- `Tutorial.ipynb` contains everything relevant on how to use.
- An update on the advance usage will be available in future

## Directories:
- `input`: directory where main input table (table_notsup.csv) are stored
- `data`: directory where all the data sets generated from the input table are stored.
- `results`: almost every use of the main function can be set to save the output tables as csv files which is stored here

## Citation
If you use this software please cite the original publication on the framework used in this package: <br>
> "The ABC of RPV: Classification of R-Parity Violating Signatures at the LHC for Small Couplings" ([arXiv:2306.07317](https://arxiv.org/abs/2306.07317)) <br>
> Herbi K. Dreiner, Yong Sheng Koay, Dominik Köhler, Víctor Martín Lozano, Javier Montejo Berlingen, Saurabh Nangia and Nadja Strobbe.

and the full manual: <br>
> TBP

## Contributors:

Original author (Code):
- [Yong Sheng, Koay](https://www.katalog.uu.se/profile/?id=N23-470)

Contributors:
- Feel free to get involved! 

## Future updates:
### Near Future 
- Include the possibility of any LSP pairing (non-SU(2) doublet). 
- Include the possibility of produced sparticles being different from LSP. Signatures from the cascade of produced sparticles to LSP will be included 
- Incorporate it into a pip-installable package

### Future 
- Include charges of final state objects  (Need to rewrite code base, in progress)
- Using all possible exact vertices instead of not_suppressed table (Need to rewrite code base, in progress)
- Generate Madgraph template param_card for parameter scans (Need to rewrite code base, in progress)
- Include the possibility of having intermediate 3-nonSM-field-vertex
- Incorporate code with large RPV couplings 



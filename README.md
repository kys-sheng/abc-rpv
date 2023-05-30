# abc-rpv
`abc-rpv` is a Python library that provides a framework for analyzing the collider signatures of RPV-MSSM studied in this paper and beyond. Users are provided with various functionalities to explore the landscape of RPV-MSSM physics within the context of small RPV-couplings.

To download library:
1. Go to green color button (Code): 
2. Download this git repository either via download as a zip 
3. or running `git clone https://github.com/kys-sheng/abc-rpv.git` 

To learn how to use:
1. Try Tutorial.ipynb

Directories:
- `input`: directory where main input table (table_notsup.csv) are stored
- `data`: directory where all the data sets generated from the input table are stored.
- `results`: almost every use of the main function can be set to save the output tables as csv files which is stored here

Feel free to provide feedback via the "Issues" section of the git repo.

## Requirements
Used Python libraries:<br>
itertools: <br>
`  pip3 install itertools or pip install itertools  `<br>
numpy: <br>
`  pip3 install numpy or pip install numpy  `<br>
os: <br>
`  pip3 install os or pip install os  `<br>
pandas: <br>
`  pip3 install pandas or pip install pandas  `<br>
subprocess: <br>
`  pip3 install subprocess or pip install subprocess  `<br>
warnings: <br>
`  pip3 install warnings or pip install warnings  `<br>

## Future updates:
### Near Future
- Include the possibility of produced sparticles being different from LSP. Signatures from the cascade of produced sparticles to LSP will be included (Next Update)
- Incorporate it into a pip-installable package

### Future (Need to rewrite code base, in progress)
- Include charges of final state objects 
- Using all possible exact vertices instead of not_suppressed table 
- Generate Madgraph template param_card for parameter scans 
- Incorporate code with large RPV cases (Depends on how that project goes)

Suggestions and feedbacks are welcomed!

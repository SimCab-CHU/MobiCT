**In this folder there is:**
- a python file (generate_interface_MobiCT.py) used to generate an html interface using the outputs of MobiCT
- a compressed *ressources* folder (ressources.zip) used for the html outputs

# Quick Start

1. Download the python files
2. Install the required dependencies
3. Download and unzip the *ressources* file (the **ressources** directory must be in the **/path/to/the/output/directory/**)
4. Execute generate_interface_MobiCT.py
```
python generate_report_MobiCT-vf.py -i /path/to/MobiCT/ouput/directory/ -k /path/to/the/cov_vaf_probs/csv/file -w /path/to/the/whitelist/txt/file -b /path/to/the/bed/file -t /path/to/the/template.html/file -o /output/path/
```
6. Enjoy!

**In this folder there is:**
- a python file (janitor_MobiCT.py) used to organize the output files of MobiCT and remove the intermediary files and clean the *workdir*
- a python file (generate_interface_MobiCT.py) used to generate an html interface using the outputs of MobiCT
- a compressed *ressources* folder (ressources.zip) used for the html outputs

# Quick Start

1. Download the python files
2. Install the required dependencies
3. Download and unzip the *ressources* file (the **ressources** directory must be in the **/path/to/the/output/directory/**)
4. Execute generate_interface_MobiCT.py
```
python generate_interface_MobiCT.py -v /path/to/the/output/directory/variants/ -s /path/to/the/output/directory/stats/ -o /path/to/the/output/directory/report/
```
6. Enjoy!

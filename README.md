# Main folder

Name| Langage | Description
:---: | :---: | :---
*Analyze_field_galaxies.ipynb* | Jupyter notebook | Notebook in which the kinematical and morphological properties of the field galaxies within our sample are analyzed.
*Analyze_cluster_galaxies.ipynb* | Jupyter notebook | Looks at the difference in the morphological properties of the selected cluster galaxies between GalFit and Cassata/Tasca catalogs.
*Check_RA_RA2000_Dec_Dec2000.ipynb* | Jupyter notebook | Notebook which checks that (RA, DEC) wihtin the master .vot files are identical to the given (RA2000, DEC2000).
*Convert_tbl_to_vot.ipynb* | Jupyter notebook | Notebook which converts .tbl files into .vot files for TopCat.
*extraction_gals_test.ipynb* | Jupyter notebook | A notebook to test astropy and numpy functionnalities on structured arrays, masked arrays and data from .vot files.
*extraction_gals_whole.ipynb* | Jupyter notebook | Main program which builds up master catalogs from those in data/catalogues by unifying the fields names and their data type which change from one file to another.
*Plot_differences_between_catalogs.ipynb* | Jupyter notebook | Contains all the different plots comparing GalFit properties such as radius and b/a with Cassata, Tasca and Zurich catalogs.
*Selecting_field_Gals.ipynb* | Jupyter notebook | File within which the SFR = f(Mass) is inversitgated and galaxies from MUSE are selected according to the graphs from *Plot_differences_between_catalogs.ipynb* and in Plots folder.
*stage2.py* | Python | Library with useful functions for the internship, namely computing and applying masks on many arrays in two lines, plotting nice and highly configurable graphs with one command, finding occurences of some value in many arrays, etc.

# scripts_python_Benoit

## Programs list

Name| Langage | Description
:---: | :---: | :---
*build_folder_structure* | Bash | Twofold program. Firstly, builds a tree structure in the given folder (default is ../outputs/MUSE) based on the MUSE file structure found in ../data. Secondly, generate a file (default name is list_output_folders) with the newly created folders. If the structure already exists, nothing is created or overwritten.
*clean_data.py* | Python | Program which cleans MUSE data using a threshold in velocity dispersion, SNR and potentially a manually applied mask.
*generate_list_gal* | Bash | Look for all the galaxies in a given folder (default is ../data). Galaxies must have a .config file in their folder in order to be listed. The outputfile (default file name is list_gal) contains in the first column the full names (path+name) of the galaxies, and their redshift in the second column.

## How to use

If you run clean_data.py for the first time, you must first generate the list of galaxies and then build the folder structure for the ouput files.

First, generate the list of galaxies by running

```bash
./generate_list_gal 
```

or if you want to change the ouput file name and/or the folder location (where the MUSE data is stored), use

```bash
./generate_list_gal outputName folder
```

where outputName is the name of the ouput file, and folder is the MUSE data folder location.

Secondly, build the output folder structure for the cleaned maps. To do so, run

```bash
./build_folder_structure
```

or if you want to use a list of galaxies from a different file and/or change the main folder name where the structure will be built (default is ../outputs/MUSE), use instead

```bash
./build_folder_structure inputName folder
```

where inputName is the file with all the galaxies we want to clean and folder is the location of the main folder where the structure will be built.

## I/O in clean_data.py

### Input files

Clean_data.py uses a file as an input in order to know which galaxy it must clean. Available filenames are listed in the dictionnary *possibleNames*. 

The input text file must have the following structure. The first column contains the full relative filename to the .config file of every galaxy we want to clean. The second column must have the redshift of the galaxy, which can be found in the .config file.

### Output files

Two ouput files are created in addition to the cleaned maps. Those files appear in the same folder as clean_data.py.

The first one is the clean_o2 file which has the same structure as the input file but with the lsfw instead of redshift. The second one is the outputFolders file which lists every folder where a cleaned map has been made.

## Comparing maps before and after cleaning

To compare the maps with PyQubeVis before and after cleaning, *checkGalsBefore* and *checkGalsAfter* programs can be used. These require PyQubeVis to be in the path under this specific name.

To open all the newly cleaned maps in PyQubeVis, run the following command

```bash
./checkGalsAfter filename |extension] [-b] 
```

where filename is the outputFolders file. Parameters in square brackets are optional, with extension refering to the kind of file we want to open and -b to open the files in buffer mode. Default is to open velocity field maps un unbuffered mode.

For instance, the following command

```bash
./checkGalsAfter listGal _ssmooth_vel_common_clean5.0.fits -b
```

will open all the galaxies listed in listGal, finishing with the extension \_ssmooth_vel_common_clean5.0.fits in a single PyQubeVis window (buffered mode). In order for each galaxy to have its own window, remove -b.

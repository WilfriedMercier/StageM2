# *scripts\_python\_Benoit*

## Programs list

Name| Langage | Description
:---: | :---: | :---
*build\_folder\_structure* | Bash | Twofold program. Firstly, builds a tree structure in the given folder (default is *[../outputs/](https://github.com/WilfriedMercier/StageM2/tree/master/outputs)MUSE*) based on the MUSE file structure found in *[../data](https://github.com/WilfriedMercier/StageM2/tree/master/data)*. Secondly, generate a file (default name is *list\_output\_folders*) with the newly created folders. If the structure already exists, nothing is created or overwritten.
*checkGalsAfter* | Bash | Open the maps of all the cleaned galaxies listed in a file with PyQubeVis.
*checkGalsBefore* | Bash | Open the maps of all the galaxies before they were cleaned (galaxies listed in the input file for the cleaning program) with PyQubeVis
*clean\_data.py* | Python | Program which cleans MUSE data using a threshold in velocity dispersion, SNR and potentially a manually applied mask.
*create\_hst\_stamps.py* | Python | This is used to create the close up hst images of the galaxies. *.txt* input files with galaxy information must be made beforehand either manually or with *Create\_hst\_stamps\_input.ipynb*. They can be found in *[data/hst](https://github.com/WilfriedMercier/StageM2/tree/master/data/hst)/CGr\$**.
*Create\_hst\_stamps\_input.ipynb* | Notebook | This notebook automatically constructs the input *.txt* files used by *create\_hst\_stamps.py* to make the close up hst images of the selected galaxies. It relies on the *.vot* files within *[outputs/SelectedGals\_sep\_by\_cluster/](https://github.com/WilfriedMercier/StageM2/tree/master/outputs/SelectedGals_sep_by_cluster)CGr\**.
*createSymLinks\_to\_ouputfolders* | Bash | Create symbolic links in the galaxies outputfolders to the original *.fits* files for the kinematical modelling.
*kinemorpho_catalogs.py* | Python | Used to generate the recap files with the kinematical and morphological information of the modelled galaxies
*Generate\_kinematical\_input.ipynb* | Notebook | This notebook gathers morphological information of the selected galaxies from different files and combine them into objects of a class named groupStructure. These objects are then used to easily create input files such as *input\_fit\_o2.txt* or *maps\_input\_o2.txt*.
*generate\_list\_gal* | Bash | Look for all the galaxies in a given folder (default is *[../data](https://github.com/WilfriedMercier/StageM2/tree/master/data)*). Galaxies must have a *.config* file in their folder in order to be listed. The outputfile (default file name is *list_gal*) contains in the first column the full names (path+name) of the galaxies, and their redshift in the second column.
*usePyQubeVisEnMasse* | Bash | This opens one after another the cube file and the automatically cleaned velocity map (assuming a threshold in SNR of 5, this can be changed directly in the code) in a single PyQubeVis session. This is meant to help to manually clean the maps.

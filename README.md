Name| Langage | Description
:---: | :---: | :---
*Analyze\_field\_galaxies.ipynb* | Jupyter notebook | Notebook in which the kinematical and morphological properties of the field galaxies within our sample are analyzed.
*Analyze\_cluster\_galaxies.ipynb* | Jupyter notebook | Looks at the difference in the morphological properties of the selected cluster galaxies between GalFit and Cassata/Tasca catalogs.
*Check\_morphology.ipynb* | Jupyter notebook | This is the notebook within which morphological and kinematical information of the selected galaxies are investigated.
*Check\_RA\_RA2000\_Dec\_Dec2000.ipynb* | Jupyter notebook | Notebook which checks that (RA, DEC) wihtin the master .vot files are identical to the given (RA2000, DEC2000).
*Convert\_tbl\_to\_vot.ipynb* | Jupyter notebook | Notebook which converts *.tbl* files into *.vot* files for TopCat.
*extraction\_gals\_test.ipynb* | Jupyter notebook | A notebook to test astropy and numpy functionnalities on structured arrays, masked arrays and data from .vot files.
*extraction\_gals\_whole.ipynb* | Jupyter notebook | Main program which builds up master catalogs from those in *[data/catalogues/](https://github.com/WilfriedMercier/StageM2/tree/master/data/catalogues)* by unifying the fields names and their data type which change from one file to another.
*NewCatalogWithCorrectedRadius.ipynb* | Notebook | Correct Cassata radius bias and write into a VOtable file the new catalog with corrected radius (either corrected Cassata radius or Zurich one if available)
*Plot_differences\_between\_catalogs.ipynb* | Jupyter notebook | Contains all the different plots comparing GalFit properties such as radius and b/a with Cassata, Tasca and Zurich catalogs.
*stage2.py* | Python | Library with useful functions for the internship, namely computing and applying masks on many arrays in two lines, plotting nice and highly configurable graphs with one command, finding occurences of some value in many arrays, etc.
*Selecting\_field\_Gals.ipynb* | Jupyter notebook | File within which the SFR = f(Mass) is inversitgated and galaxies from MUSE are selected according to the graphs from *Plot\_differences\_between\_catalogs.ipynb* and in *[Plots](https://github.com/WilfriedMercier/StageM2/tree/master/Plots)* folder.

## Preliminary

### Building the folder structure

If you run clean_data.py for the first time, you must first generate the list of galaxies and then build the folder structure for the ouput files.

First, generate the list of galaxies by running

```console
wilfried:~/ST2$ ./generate_list_gal 
```

or if you want to change the ouput file name and/or the folder location (where the MUSE data is stored), use

```console
wilfried:~/ST2$ ./generate_list_gal outputName folder
```

where outputName is the name of the ouput file, and folder is the MUSE data folder location.

Secondly, build the output folder structure for the cleaned maps. To do so, run

```console
wilfried:~/ST2$ ./build_folder_structure
```

or if you want to use a list of galaxies from a different file and/or change the main folder name where the structure will be built (default is ../outputs/MUSE), use instead

```console
wilfried:~/ST2$ ./build_folder_structure inputName folder
```

where inputName is the file with all the galaxies we want to clean and folder is the location of the main folder where the structure will be built.

### Creating symbolic links

For later analysis purposes, it can be useful to have all the *.fits* files from the data folder within the MUSE ouput folder where the cleaned maps and the model will be saved.

To automatically create symbolic links between the MUSE ouput folder and the *.fits* files in the data folder, use

```console
wilfried:~/ST2$ ./createSymLinks_to_ouputfolders
```

__Note__ : this program only works with the default folder structure (i.e. input data in the corresponding data/CGr* folders, outputs files in the outputs/MUSE/CGr* folders.


## Cleaning galaxies

### Automatic cleaning

*clean_data.py* uses an input file to know which galaxies it must clean. To add a file you would like to use or to modify the name of an existing one, look at the dictionnary ```possibleNames``` in the program. 

The input text file must have the following structure

- the first column contains the full relative filename to the *.config* file of every galaxy we want to clean

- the second column must have the redshift of the galaxy (found in the *.config* file)

The SNR and velocity dispersion thresholds can be modified directly from the python file.
Once everything is set up, simply run *clean_data.py* to clean the galaxies.

__Note__: two ouput files are created in addition to the cleaned maps. Those files appear in the same folder as *clean\_data.py*. The first one is the *clean\_o2* file which has the same structure as the input file but with the LSF FWHM instead of redshift. The second one is the *\_outputFolders* file which lists every folder where a cleaned map has been made.

### Manual cleaning

To manually clean the remaining pixels which do not belong to the galaxy there are two options

1. Move to the relevant output folder and open the cube file with PyQubeVis and then the automatically cleaned velocity map from the software menu. Once the map has been manually cleaned, save it as *clean.fits* in the same folder as the cleaned maps and the symbolic links. 

2. Run the program `wilfried:~/ST2$ ./usePyQubeVisEnMasse outputFolderFile` where `outputFolderFile` is the file made by *clean\_data.py* containing the folders where the cleaned galaxies are. This will open a new PyQubeVis session with the data cube and the automatically cleaned velocity map one after another. Thus, every time the session will be closed (after the manual cleaning and saving as *clean.fits*) a new session with the following galaxy will be opened.

To apply the modifications to all the maps, move back to *[scripts\_python\_Benoit](https://github.com/WilfriedMercier/StageM2/tree/master/scripts_python_Benoit)* and re-run *clean_data.py* without making any modifications. This will overwrite the previous cleaned maps, taking into account the information given in the *clean.fits* map.

## Comparing maps before and after cleaning

### Standard procedure

To check that the cleaning procedure did work correctly, the usual way is to move into the ouput folder (listed in the *\_outputFolders* file) and to open in two different PyQueVis sessions the velocity maps before and after cleaning (usually ending with *\_cleanNB.fits* where NB is one of the thresholds used for the cleaning).

This then allows the user to open the cube file (ending with *\_ssmooth_cube.fits*) to check the lines.

__Note__ : in order to open the maps before cleaning from the ouput folder, one has to create the relevant symbolic links to the input *.fits* files in the input folder either manually or automatically (using *createSymLinks\_to\_ouputfolders*).

### Quick look for large numbers of files

If the cleaning procedure was performed on a large amount of galaxies in a single run, two programs can be used to quickly check that the resulting maps look fine.

To compare all the maps with PyQubeVis before and after cleaning, *checkGalsBefore* and *checkGalsAfter* programs can be used.

To open all the newly cleaned maps in PyQubeVis, run the following command

```console
wilfried:~/ST2$ ./checkGalsAfter filename [extension] [-b] 
```

where filename is the *\_outputFolders* file. Parameters in square brackets are optional, with `extension` refering to the kind of file we want to open and `-b` to open the files in buffer mode. Default is to open velocity maps in unbuffered mode.

As an example, the following command

```console
wilfried:~/ST2$ ./checkGalsAfter listGal _ssmooth_vel_common_clean5.0.fits -b
```

will open all the galaxies listed in *listGal*, finishing with the extension *\_ssmooth\_vel\_common\_clean5.0.fits* in a single PyQubeVis window (buffer mode). In order for each galaxy to have its own window, remove the `-b` option.

Maps before cleaning can be opened in the same way with the command


```console
wilfried:~/ST2$ ./checkGalsBefore filename [extension] [-b] 
```

where `filename` should be in this case the same file as the input file used for cleaning the maps.

## Fitting a model to the data

### Manually generating the input file

Once the maps have been manually cleaned, move to the corresponding CGr* group folder in *[outputs/](https://github.com/WilfriedMercier/StageM2/tree/master/outputs)MUSE* and either manually create a file named *input\_fit\_o2.txt* with the following structure:

| \# |      ID     |   X     |   Y |      PA     | INC     | vs    | vm    | d  | sig | psfx  | psfz | smooth
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: 
| | CGr51\_28\_o2  |  16.16 |   15.81 |  0.0 |   45 |   0.0  |   80 |  2.0  |   0 |  2.45  | 51.6 |   2

In this file, each line represents a galaxy with the following properties:

Table column | Description | How to determine ?
:---: | :---: | :---
ID | The galaxy folder name (within the o2 folder). | 
X | x-position of the centre of the galaxy in the MUSE image (given in pixels)  | Determined from the hst stamp image (see next section) but a conversion factor must be applied since hst stamps are larger than MUSE images (generally 200x200 px instead of 36x36 px). 
Y | y-position of the centre in the MUSE image (given in pixels).  | Derived in the same way as `X`.
PA | Position angle of the morphological major axis (with respect to North) given in the range [-90°, 90°]. | Sometimes given by Zurich in the catalogs.
INC | Inclination of the galaxy (taken such that sin(INC) = ellipticity) given in degrees.| Sometimes derived from Zurich entries in the catalogs.
psfx | Spatial PSF in pixels | Since this is wavelength and group dependent, the `computeGroupFWHM` function in *stage2.py* can be used to derive its value. This requires to give the observed wavelength of the spectral feature (rest frame wavelength x (1+z)) as well as the group number. 
psfz | LSF FWHM in km/s | Found in the *_o2* file in *[scripts\_python\_Benoit](https://github.com/WilfriedMercier/StageM2/tree/master/scripts_python_Benoit)*.

__Notes__:

- the other parameters should not be modified
- `X`, `Y` and `INC` are fixed parameters and must therefore be tightly constrained

### Automatically generating the input file

If you have selected the galaxies regardless of their group, you can easily generate the input files for all the groups with *Generate\_kinematical\_input.ipynb*. 

This notebook uses a class named groupStructure which gathers the location and the morphological information (including PSF and LSF FWHMs) of all the galaxies in the same group. Therefore, to each group corresponds an object with all the galaxy information which is then used to easily create the input file for each group.

### Fitting the model

Once the input file is made and if the program for the model fitting is not present, you can either copy in the current folder fitModel.pro (found in [PSF_data/](https://github.com/WilfriedMercier/StageM2/tree/master/PSF_data)) and run it in IDL with the two following lines

```idl
.run fitModel
fitModel
```

or directly copy and run the following lines in IDL:

```idl
fito=2
version='1.0'
model='slp'
options='_mclean5.0'
option1='_ssmooth'
file='input_fit_o2.txt'
path='o2/'
line='_OII3729'
fit_massiv, fito, /ifix, pfix=0, /xfix, /yfix, plot=0, file=file, options=options, option1=option1, path=path, line=line
```

If everything went fine, this should have created model files (*modd*, *modhr* and *modv*), residual maps (*resd* and *resv*), as well as three recap files with model parameters and residual information:

- *CGr\*\_parameters\_red_slp\_xyi\_mclean5.0.txt*
- *CGr\*\_parameters\_residual\_slp\_xyi\_mclean5.0.txt*
- *CGr\*\_o2\_vmax\_map_rlast\_mclean5.0.txt*

## Producing recap files

Some kinematical information are used when creating the input file for the maps. To gather the information for each galaxy into a single file, run *kinemorpho_catalogs.py* in *[scripts\_python\_Benoit](https://github.com/WilfriedMercier/StageM2/tree/master/scripts_python_Benoit)*.

For each group, a new file will be made in the o2 folder.  The groups names for which you want to make a recap file can be specified through  `groups` variable.

## Producing the input file for the maps

The idl routine producing the maps relies on a single input file named *maps\_input\_o2.txt* which is found in each group directory. This file has the following structure

| \# | name | ID | morpho\_image | opt1 |  line |  opt2 |   x |       y |    pa      | i |    vs  | zpsf     | z     | r12   | psf | | | |  | 
| :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-:  | :-: | :-: | :-: |
| | CGr23\_13\_o2   | CGr23\_13 |   /home/wilfried/ST2/data/hst/CGr23/13\_CGr23.fits  | mclean5.0  | OII3729 |  xyi |  15.0 |  15.0 |  -104.2  | 52 |   0.0 |  52.7 |  0.677 |  1.0 |  0.831 |  0  | RDPRMS  | 0 |  0
| | CGr23\_19\_o2 |  CGr23\_19 |  /home/wilfried/ST2/data/hst/CGr23/19\_CGr23.fits |  mclean5.0 |  OII3729 |  xyi |  16.5  | 15.0 |  9.1 |   59 | 0.0 |  43.1 |  0.953  | 1.0 |  0.772 |  0 |  RDPRMS  | 0 |  0
| ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... |
| | CGr23\_25\_o2 |  CGr23\_25 |  /home/wilfried/ST2/data/hst/CGr23/25\_CGr23.fits |  mclean5.0 |  OII3729 |  xyi |  15.0 |  15.0 |  235.5 |   57 |   0.0 |  52.7 |  0.677 |  1.0 |  0.831 |  0 |  RDPRMS |  0 |   0

where `name` is the galaxy folder name, `ID` is the name which will appear on the maps, `morpho_image` is the location of the HST morphological close up image, `opt1` is the suffixe of the files, `line` is the name of the emission line used, `opt2` is the fixed parameters (x, y centre coordinates, i for inclination, etc), `x` and `y` are the galaxy centre coordinates in pixels in the MUSE image,  `pa` is the kinematical PA, `i` is the inclination, `zpsf` is the LSF FWHM in km/s, `z` is the redshift and `psf` is the PSF FWHM in arcsec. Other columns always keep these values. 

You can either manually generate it or you can use the corresponding cell in [scripts\_python\_Benoit](https://github.com/WilfriedMercier/StageM2/tree/master/scripts_python_Benoit)/Generate\_kinematical\_input.ipynb jupyter notebook.

## Generating the maps





## Making HST images of galaxies

It can be useful to have HST close-up images of the selected galaxies. This can be done with *create\_hst\_stamps.py*.

This program will use both *.txt* input files and the *HST\_CGr\*.fits* files found in the *[data/hst/](https://github.com/WilfriedMercier/StageM2/tree/master/data/hst)CGr\** folders.

Within each group folder, a *.txt* file must be built with the following structure:

ID  | ID_LAIGLE  |      z |  flag       |   ra   |    dec  |      Iab 
:---: | :---: | :---: | :---: | :---: | :---: | :---:
18  |   628424 | 0.37964   |  3 |  150.006272 | 2.253330 | 21.667969 
21   |  628636 | 0.34542    | 3 | 149.992574  | 2.252479  |23.621355 
... | ... | ... | ... | ... | ... | ...
23   |  628677 | 0.83942   |  3 | 150.006706 | 2.253204 | 21.654572


where `ID` is the galaxy MUSE ID, `ID_LAIGLE` is the galaxy ID from [Laigle+16](https://ui-adsabs-harvard-edu.ezproxy.obspm.fr/#abs/2016ApJS..224...24L/abstract) catalogue, `z` is the redshift, `flag` is the MUSE CONFID flag (confidence level in redshift value), `ra` and `dec` are the galaxy position, and `Iab` is the galaxy magnitude (generally i++ mag.).

__Note__: Column names must be exactly these ones.

Such *.txt* files can either be made manually or can be automatically generated using *Create\_hst\_stamps\_input.ipynb*. 

This notebook will use the *.vot* files in *[outputs/SelectedGals\_sep\_by\_cluster/](https://github.com/WilfriedMercier/StageM2/tree/master/outputs/SelectedGals_sep_by_cluster)CGr\** folders which contain the selected galaxies in the corresponding cluster.

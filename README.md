# FSGLM
##A wrapper for FreeSurfer's mri_glmfit.
the standard way to do FreeSurfer whole-brain analyses is a sequence of steps, which include some file creation, mris_preproc, mri_surf2surf, mri_glmfit, and mri_glmfit-sim. Since many of these steps are directly related to the analysis you want to run, and are standardized, a python wrapper (FSGLM.py) for those steps has been created to minimize redundancy. It also generates 4 images of the non-clusterwise-corrected, and clusterwise-corrected template brain for easy review. To simplify things, not all options of mri_glmfit are allowed in FSGLM, but it will likely meet most if not all needs.

Options:
Running the script with the `-h` option will output a help page, which lists the information below.
Usage:
`FSGLM.py [options] --input=inputfile --projectdir=ProjectDir --subjectdir=SubjectDir --hemi=hemi [ -c ( --var=var | -l [ --T1var=T1var --T2var=T2var ] ) [ --controlfor=var... ] ] [ --groupcol=groupcolumn --g1=group1 --g2=group2 ]`

Breaking this down, you must provide an input file (more on that later), specify the project directory (where you want the glm directory to be), specify the subject directory (SUBJECTS_DIR), and the hemisphere (left or right). The other options are technically not required, but I would recommend putting something there, because otherwise the glm is pretty meaningless.

###Single-Timepoint (`--var=var`):
Specify the column in the input file you want to include in the model as a covariate. For example, you may be interested in seeing how Psychological Well-Being is related to cortical thickness across all subjects. If the column name of your input file is PWB, you would write `--var=PWB`.

###Multi-Timepoint (without covariates):
Running a multi-timepoint analysis without covariates tests for differences between time points in the measure of interest. In this case, you would specify `-l`, but not any variables.

###Multi-Timepoint (with covariates: `--T1var=T1var --T2var=T2var`):
Specify the columns in the input file you want to include in the model. For simplicity, and to protect against regression towards the mean, these two variables are automatically converted into residualized change scores, and that single value is added to the model as a covariate. If you want it calculated as a simple change score, you will have to do that yourself at this time, but let me know if a large number of you want it, and I can add it. An example of running this analysis is to see if an increase in Psychological Well-Being is related to less decrease in cortical thickness from T1 to T2.

###Single-Group:
By default, the GLM will run without groups. See above for examples.

###Multi-Group (`--groupcol=groupcolumn --g1=group1 --g2=group2`):
Because it is not simple to run a group analysis with greater than 2 levels, I capped it at this for ease of coding. Therefore, you specify the column containing group membership information, and then the two groups you are interested in comparing. Suppose you wanted to compare two groups (e.g. non-meditators and meditators), and that information is included in a column called "Group". Supposing the two groups were CatLover and DogLover, you would specify `--groupcol=Group --g1=CatLover --g2=DogLover`. Exchanging `--g1` and `--g2` does not affect the result, but will flip the sign. This will run an analysis to see where/if CatLovers had different cortical thickness than DogLovers. Adding a covariate will test an interaction between the group and the covariate.

###Multi-Group, Multi-Timepoint (without covariate):
Combining the above options allows you to run a multi-group, multi-timepoint analysis, which would see if one group's change in brain structure changed differently from the other group.

###Multi-Group, Multi-Timepoint (with covariate):
Combining the above options allows you to run a multi-group, multi-timepoint analysis, which would see if one group's relationship between brain structure and a measure of interest changed compared to the other group.

###Controlling for variables:
By specifying headers in the input csv with --controlfor=, you can control for a continuous variable. You can specify multiple variables to control for by repeating the option (e.g. `--controlfor=ControlVar1, --controlfor=ControlVar2, --controlfor=ControlVar3 ...`).

###Boolean Selection of Data:
You can add an optional boolean column to the input csv titled "Include", which has values of either N or Y. Rows with N will not be included in the glm. By default, all rows will be considered.

##Requirements:
* Python 2.7
* [FreeSurfer](https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferWiki)
* [Pandas (and dependencies)](http://pandas.pydata.org)
* [StatsModels (and dependencies)](https://github.com/statsmodels/statsmodels)

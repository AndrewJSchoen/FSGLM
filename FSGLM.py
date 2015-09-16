#!/usr/bin/env python

Version = "0.1"
doc = """
Run FS GLM.

Usage:
  FSGLM.py [options] --input=inputfile --projectdir=ProjectDir --subjectdir=SubjectDir --hemi=hemi [ -c ( --var=var | -l [ --T1var=T1var --T2var=T2var ] ) [ --controlfor=var... ] ] [ --groupcol=groupcolumn --g1=group1 --g2=group2 ]

Options:
  -h --help                                 Show this screen.
  -v --version                              Show the current version.
  --meas=measure                            Freesurfer measure (thickness, pial, area, curv, etc.). [default: thickness]
  --hemi=hemi                               Hemisphere to analyze (upper/lowercase versions of lh/left/l or rh/right/r).
  --smooth=smooth                           Amount of smoothing, using a gaussian kernel with the given value specified in mm. [default: 5]
  --thresh=threshold                        Cluster-forming threshold (choose from 0.05, 0.01, 0.005, 0.001, 0.0005, or 0.0001). [default: 0.05]
  -l --longitudinal                         Specify longitudinal processing. Required if specifying T1var and T2var. [default: False]
  -c --covariate                            Include a covariate in the model. [default: False]
  --var=var                                 Var for non-longitudinal analysis. [default: False]
  --T1var=T1var                             Var for T1 longitudinal analysis. [default: False]
  --T2var=T2var                             Var for T1 longitudinal analysis. [default: False]
  --controlfor=var...                       Var to control for. Specify multiple variables by repeating --controlfor. [default: False]
  --groupcol=groupcolumn                    Group membership column in csv. [default: False]
  --g1=group1                               Group 1 name in csv. [default: False]
  --g2=group2                               Group 2 name in csv. [default: False]
  -i=inputfile --input=inputfile            Input csv file. Expects a column "Subjects", as well as some data columns. [default: False]
  -s=SubjectDir --subjectdir=SubjectDir     SUBJECTS_DIR. [default: False]
  -p=ProjectDir --projectdir=ProjectDir     Analysis processing dir for this glm. [default: False]
"""


import csv, sys, os, shutil, subprocess, math
import pandas
sys.path.append("lib/docopt")
from statsmodels import api as sm
from docopt import docopt

#============================================================================
#============ Functions =====================================================

def exists(path):
  if os.path.exists(os.path.normpath(path)):
      return(1)
  else:
      return(0)

def systemCall(command):
  p = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True)
  return p.stdout.read()

def parseCSV(csvfilepath, headerlist=[]):
  if exists(csvfilepath):
    csvfilepath = os.path.normpath(csvfilepath)
    print("Parsing CSV File '{0}'.".format(csvfilepath))
    csvcontent = pandas.read_csv(csvfilepath, converters={'Subject': lambda x: str(x)})
    csvcontent.index = csvcontent["Subject"]
    csvcontent.drop("Subject", axis=1, inplace=True)
    if "Include" in list(csvcontent.columns.values):
        headerlist.append("Include")
    if len(headerlist) > 0:
        csvcontent=csvcontent[headerlist]
    return csvcontent
  else:
      print("CSV File '{0}' does not exist! Exiting now.".format(csvfilepath))
      sys.exit(1)

def writeCSV(csv, csvfilepath):
  csv.to_csv(csvfilepath)

def writeFile(text, filename):
    file = open(filename, 'a')
    file.write("{0}".format(text))
    file.close

def trimCSV(csv, headers):
  newCSV = csv[headers]
  return newCSV

def removeRowsWithCriteria(csv, column, value):
  print("Revmoving cases where {0} = {1}".format(column, value))
  csv = csv[csv[column] != value]
  return csv

def cleanargs(args):
  clean = {}

  #Whether to process longitudinally
  if args["--longitudinal"] == True or args["--longitudinal"] == "True":
      clean["IsLongitudinal"] = True
  else:
      clean["IsLongitudinal"] = False

  #Whether to process with covariates
  if args["--covariate"] == False or args["--covariate"] == "False":
      clean["TestCovariates"] = False
      clean["HasMainVar"] = False
  else:
      clean["TestCovariates"] = True

      if args["--var"] != False and args["--var"] != "False" and args["--longitudinal"] == False:
          clean["MainVar"] = args["--var"]
          clean["HasMainVar"] = True
      elif args["--T1var"] != False and args["--T1var"] != "False" and args["--T2var"] != False and args["--T2var"] != "False" and clean["IsLongitudinal"] == True:
          clean["T1var"] = args["--T1var"]
          clean["T2var"] = args["--T2var"]
          clean["HasMainVar"] = True
      else:
          clean["HasMainVar"] = False

      validcontrolfors = list(set(args["--controlfor"]) - set(["False"]))
      if len(list(set(args["--controlfor"]) - set(["False"]))) > 0:
          clean["HasControlForVars"] = True
          clean["ControlVars"] = args["--controlfor"]
      else:
          clean["HasControlForVars"] = False

  if args["--groupcol"] != False and args["--groupcol"] != "False":
      clean["IsGroup"] = True
      clean["Group1"] = args["--g1"]
      clean["Group2"] = args["--g2"]
  else:
      clean["IsGroup"] = False

  #Conflict Testing
  if clean["TestCovariates"] == True and clean["HasControlForVars"] == False and clean["HasMainVar"] == False:
      clean["TestCovariates"] = False

  clean["InputCSV"] = cleanPathString(args["--input"])
  clean["SubjectDir"] = cleanPathString(args["--subjectdir"])
  clean["ProjectDir"] = cleanPathString(args["--projectdir"])
  clean["FreeSurferMeasure"] = args["--meas"]
  clean["Smoothing"] = args["--smooth"]

  clean["Threshold"] = str(args["--thresh"])
  if clean["Threshold"] == "0.05":
    clean["LogThreshold"] = "1.3"
  elif clean["Threshold"] == "0.01":
    clean["LogThreshold"] = "2.0"
  elif clean["Threshold"] == "0.005":
    clean["LogThreshold"] = "2.3"
  elif clean["Threshold"] == "0.001":
    clean["LogThreshold"] = "3.0"
  elif clean["Threshold"] == "0.0005":
    clean["LogThreshold"] = "3.3"
  elif clean["Threshold"] == "0.0001":
    clean["LogThreshold"] = "4.0"
  else:
    clean["LogThreshold"] = string(math.floor(math.log(float(clean["Threshold"]))*10)/10)

  hemistring = str(args["--hemi"]).lower()
  if hemistring == "l" or hemistring == "lh" or hemistring == "left":
      clean["Hemisphere"] = "lh"
  elif hemistring == "r" or hemistring == "rh" or hemistring == "right":
      clean["Hemisphere"] = "rh"
  else:
      print("Valid hemisphere (--hemi) not specified.")
      sys.exit(1)

  #Create Name
  if clean["IsGroup"]:
    name = "{0}_vs_{1}".format(clean["Group1"], clean["Group2"])
  else:
    name = "Individual"
  name += "_" + clean["Hemisphere"] + "_" + clean["FreeSurferMeasure"]
  if clean["HasMainVar"] == True:
    name += "_" + clean["MainVar"]
  if clean["HasControlForVars"]:
      for cfvar in clean["ControlVars"]:
          name += "_cf_{0}".format(cfvar)
  clean["GLMNAME"] = name
  clean["ProjectPath"] = "{0}/{1}".format(clean["ProjectDir"], clean["GLMNAME"])
  clean["GroupColumn"] = args["--groupcol"]

  return clean

def isHeader(csv, header):
  headers = csv[0].keys()
  if header in headers:
    return 1
  else:
    print("Header '{0}' does not exist in the csv.".format(header))
    return 0

def cleanPathString(path):
  if path.endswith("/"):
    path = path[:-1]
  if path.startswith("="):
    path = path[1:]
  realpath = os.path.realpath(path)
  return realpath

def printArgs(args):
  print
  print("##Arguments##")
  for key in sorted(args.keys()):
    print("{0}:\n - {1}".format(key, args[key]))
  print

def getUniqueValsInCol(csv, column):
  unique = []
  for value in csv[column]:
    if value not in unique:
      unique.append(value)
  return unique

def printCSV(csv):
  keys = sorted(csv[0].keys())
  print("\t".join(keys))
  for row in csv:
    rowlist = []
    for key in keys:
      if len(key) < len(row[key]):
        rowlist.append("...")
      else:
        rowlist.append(row[key])
    print("\t".join(rowlist))

def createFSGD(args):
    fsgdfilelist = []
    fsgdfilelist.append("GroupDescriptorFile 1")
    fsgdfilelist.append("Title {0}".format(args["GLMNAME"]))


    if args["IsGroup"]:
      for group in getUniqueValsInCol(args["csvdata"], args["GroupColumn"]):
        fsgdfilelist.append("Class {0}".format(group))
    else:
      fsgdfilelist.append("Class Main")
    if args["HasMainVar"] == True:
      fsgdfilelist.append("Variables {0}".format(args["MainVar"]))
    if args["HasControlForVars"] == True:
        for controlvar in args["ControlVars"]:
            fsgdfilelist.append("Variables {0}".format(controlvar))
    for index in list(args["csvdata"].index):
      if args["IsGroup"]:
        currentGroup = args["csvdata"][args["GroupColumn"]][index]
      else:
        currentGroup = "Main"

      fsgdfileline = "Input\t{0}".format(index)
      if args["IsLongitudinal"]:
          fsgdfileline += "_base"

      fsgdfileline += "\t{0}".format(currentGroup)

      if args["HasMainVar"] == True:
          fsgdfileline += "\t{0}".format(args["csvdata"][args["MainVar"]][index])
      if args["HasControlForVars"] == True:
          for controlvar in args["ControlVars"]:
              fsgdfileline += "\t{0}".format(args["csvdata"][controlvar][index])
      fsgdfilelist.append(fsgdfileline)
    fsgdfile = "\n".join(fsgdfilelist)
    return fsgdfile

def createPairFSGD(args):
    pairfsgdfilelist = []
    pairfsgdfilelist.append("GroupDescriptorFile 1")
    pairfsgdfilelist.append("Class Main")
    for index in list(args["csvdata"].index):
      pairfsgdfilelist.append("Input\t{0}_1.long.{0}_base\tMain".format(index))
      pairfsgdfilelist.append("Input\t{0}_2.long.{0}_base\tMain".format(index))
    pairfsgdfile = "\n".join(pairfsgdfilelist)
    return pairfsgdfile

def createMatrix(args):
    #Make mtx contrast file
    if args["IsGroup"]:
        if args["HasMainVar"]:
            mtxcontrastfile = "0 0 1 -1"
        else:
            mtxcontrastfile = "1 -1"
        if args["HasControlForVars"]:
            for controlvar in args["ControlVars"]:
                mtxcontrastfile += " 0 0"
    else:
        if args["HasMainVar"]:
            mtxcontrastfile = "0 1"
        else:
            mtxcontrastfile = "1"
        if args["HasControlForVars"]:
            for controlvar in args["ControlVars"]:
                mtxcontrastfile += " 0"
    return mtxcontrastfile

def getCommands(args):
    if args["IsLongitudinal"]:
      unsmoothedbrainfile = "{0}/{1}.paired-diff.{2}.mgh".format(args["ProjectPath"], args["Hemisphere"], args["FreeSurferMeasure"])
      smoothedbrainfile = "{0}/{1}.paired-diff.{2}.sm{3}.mgh".format(args["ProjectPath"], args["Hemisphere"], args["FreeSurferMeasure"], args["Smoothing"])
      mrispreproccommand = "SUBJECTS_DIR={0} ; mris_preproc --target fsaverage --hemi {1} --meas {2} --out {3} --fsgd {4}/pairs.fsgd --paired-diff".format(args["SubjectDir"], args["Hemisphere"], args["FreeSurferMeasure"], unsmoothedbrainfile, args["ProjectPath"])
    else:
      unsmoothedbrainfile = "{0}/{1}.{2}.mgh".format(args["ProjectPath"], args["Hemisphere"], args["FreeSurferMeasure"])
      smoothedbrainfile = "{0}/{1}.{2}.sm{3}.mgh".format(args["ProjectPath"], args["Hemisphere"], args["FreeSurferMeasure"], args["Smoothing"])
      mrispreproccommand = "SUBJECTS_DIR={0} ; mris_preproc --target fsaverage --hemi {1} --meas {2} --out {3} --fsgd {4}/{5}.fsgd".format(args["SubjectDir"], args["Hemisphere"], args["FreeSurferMeasure"], unsmoothedbrainfile, args["ProjectPath"], args["GLMNAME"])

    spatialsmoothcommand = "SUBJECTS_DIR={0} ; mri_surf2surf --s fsaverage --hemi {1} --fwhm {2} --sval {3} --tval {4}".format(args["SubjectDir"], args["Hemisphere"], args["Smoothing"], unsmoothedbrainfile, smoothedbrainfile)
    glmcommand = "SUBJECTS_DIR={0} ; mri_glmfit --glmdir {1} --y {2} --fsgd {3} --C {4} --surface fsaverage {5}".format(args["SubjectDir"], args["ProjectPath"], smoothedbrainfile, "{0}/{1}.fsgd".format(args["ProjectPath"], args["GLMNAME"]), "{0}/{1}.mtx".format(args["ProjectPath"], args["GLMNAME"]), args["Hemisphere"])
    correctcommand = "mri_glmfit-sim --glmdir {0} --cache {1} abs --cwp  {2} --2spaces".format(args["ProjectPath"], args["LogThreshold"], args["Threshold"])
    lateralviewcommanduncorrected = 'xvfb-run -s "-screen 0 1280x1024x24" --auto-servernum freeview -viewport "3d" -viewsize 1280 1024 -f {0}/fsaverage/surf/{1}.inflated:edgethickness=0:overlay={2}/{3}/sig.mgh:overlay_threshold={4},6 -cam Azimuth 0 -zoom 1.4 -ss {2}/{3}_uncorrected_dorsal.png'.format(args["SubjectDir"], args["Hemisphere"], args["ProjectPath"], args["GLMNAME"], args["LogThreshold"])
    medialviewcommanduncorrected =  'xvfb-run -s "-screen 0 1280x1024x24" --auto-servernum freeview -viewport "3d" -viewsize 1280 1024 -f {0}/fsaverage/surf/{1}.inflated:edgethickness=0:overlay={2}/{3}/sig.mgh:overlay_threshold={4},6 -cam Azimuth 180 -zoom 1.4 -ss {2}/{3}_uncorrected_medial.png'.format(args["SubjectDir"], args["Hemisphere"], args["ProjectPath"], args["GLMNAME"], args["LogThreshold"])
    lateralviewcommandcorrected = 'xvfb-run -s "-screen 0 1280x1024x24" --auto-servernum freeview -viewport "3d" -viewsize 1280 1024 -f {0}/fsaverage/surf/{1}.inflated:edgethickness=0:overlay={2}/{3}/cache.th{4}.abs.sig.cluster.mgh:overlay_threshold={5},6 -cam Azimuth 0 -zoom 1.4 -ss {2}/{3}_corrected_dorsal.png'.format(args["SubjectDir"], args["Hemisphere"], args["ProjectPath"], args["GLMNAME"], args["LogThreshold"].replace(".", ""), args["LogThreshold"])
    medialviewcommandcorrected = 'xvfb-run -s "-screen 0 1280x1024x24" --auto-servernum freeview -viewport "3d" -viewsize 1280 1024 -f {0}/fsaverage/surf/{1}.inflated:edgethickness=0:overlay={2}/{3}/cache.th{4}.abs.sig.cluster.mgh:overlay_threshold={5},6 -cam Azimuth 180 -zoom 1.4 -ss {2}/{3}_corrected_medial.png'.format(args["SubjectDir"], args["Hemisphere"], args["ProjectPath"], args["GLMNAME"], args["LogThreshold"].replace(".", ""), args["LogThreshold"])

    commands=[{"Name": "Preprocessing", "Command": mrispreproccommand},
              {"Name": "Spatial Smoothing", "Command": spatialsmoothcommand},
              {"Name": "GLM", "Command": glmcommand},
              {"Name": "Cluster-wise Correction", "Command": correctcommand},
              {"Name": "Lateral Preview Creation for Uncorrected", "Command": lateralviewcommanduncorrected},
              {"Name": "Medial Preview Creation for Uncorrected", "Command": medialviewcommanduncorrected},
              {"Name": "Lateral Preview Creation for Corrected", "Command": lateralviewcommandcorrected},
              {"Name": "Medial Preview Creation for Corrected", "Command": medialviewcommandcorrected},]
    return commands

def runCommands(commands):
    print
    for command in commands:
        print("Performing {0}".format(command["Name"]))
        print(systemCall(command["Command"]))

#============================================================================
#============ Main ==========================================================

def run(args):

  # 1: Argument Parsing
  protoarguments = docopt(doc, argv=args, version='Run FS GLM v{0}'.format(Version))
  arguments = cleanargs(protoarguments)
  printArgs(arguments)

  #Get headers
  relevantheaders = []
  if arguments["IsGroup"]:
    relevantheaders.append(arguments["GroupColumn"])
  if arguments["HasMainVar"] == True:
    if arguments["IsLongitudinal"]:
      relevantheaders.append(arguments["T1var"])
      relevantheaders.append(arguments["T2var"])
    else:
      relevantheaders.append(arguments["MainVar"])
  if arguments["HasControlForVars"]:
      relevantheaders.extend(arguments["ControlVars"])

  #Parse the csv
  arguments["csvdata"] = parseCSV(arguments["InputCSV"], relevantheaders)
  print(arguments["csvdata"])

  #Trim for relevant groups
  if arguments["IsGroup"]:
    uniquegroups=getUniqueValsInCol(arguments["csvdata"], arguments["GroupColumn"])
    uniquegroups.remove(arguments["Group1"])
    uniquegroups.remove(arguments["Group2"])
    for group in uniquegroups:
      csvdata = removeRowsWithCriteria(arguments["csvdata"], arguments["GroupColumn"], group)

  #Remove Subjects with missing data
  for header in list(arguments["csvdata"].columns.values):
    arguments["csvdata"] = arguments["csvdata"][pandas.notnull(arguments["csvdata"][header])]

  if "Include" in list(arguments["csvdata"].columns.values):
    arguments["csvdata"] = removeRowsWithCriteria(arguments["csvdata"], "Include", "N")
    relevantheaders.remove("Include")
    arguments["csvdata"] = trimCSV(arguments["csvdata"], relevantheaders)

  if arguments["HasMainVar"] == True:
    if arguments["IsLongitudinal"]:
      #Calculate Residualized change scores
      Y = arguments["csvdata"][arguments["T2var"]]
      X = arguments["csvdata"][arguments["T1var"]]
      X = sm.add_constant(X)

      model = sm.OLS(Y,X)
      results = model.fit()
      intercept = results.params[0]
      slope = results.params[1]

      arguments["csvdata"][arguments["MainVar"]] = arguments["csvdata"][arguments["T2var"]] - (arguments["csvdata"][arguments["T1var"]] * slope + intercept)
      relevantheaders.remove(arguments["T1var"])
      relevantheaders.remove(arguments["T2var"])
      relevantheaders.append(arguments["MainVar"])
      arguments["csvdata"] = trimCSV(arguments["csvdata"], relevantheaders)

  #Make the project directory.
  if exists(arguments["ProjectPath"]):
    shutil.rmtree(arguments["ProjectPath"])
  os.mkdir(arguments["ProjectPath"])

  #Make paired fsgd files, if longitudinal:
  if arguments["IsLongitudinal"]:
      writeFile(createPairFSGD(arguments), "{0}/pairs.fsgd".format(arguments["ProjectPath"]))

  #Make standard fsgd file:
  writeFile(createFSGD(arguments), "{0}/{1}.fsgd".format(arguments["ProjectPath"], arguments["GLMNAME"]))

  #Make the matrix file:
  writeFile(createMatrix(arguments), "{0}/{1}.mtx".format(arguments["ProjectPath"], arguments["GLMNAME"]))

  runCommands(getCommands(arguments))

#============================================================================
#============ Entry =========================================================

if __name__ == '__main__':
    args = sys.argv
    del args[0]
    run(args)

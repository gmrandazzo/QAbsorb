'''
@package molutils.py

molutils.py was writen by Giuseppe Marco Randazzo <gmrandazzo@gmail.com>
and is distributed under LGPL version 3

Geneve November 2016
'''

import os
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import Draw
from rdkit.ML.Descriptors import MoleculeDescriptors
from sklearn.externals import joblib
import numpy as np

def GenDepiction(smi,  molname, path, height=200, width=200):
  molecule = Chem.MolFromSmiles(smi)
  AllChem.Compute2DCoords(molecule)
  imgname = molname+".png"
  imgname = imgname.replace(' ','_') # replace space with _
  imgname = imgname.replace('/','-') # replace / with -
  imagePath = path+"/"+imgname
  Draw.MolToFile(molecule, imagePath, (height,width))
  return imagePath

def LoadModels(modpath):
  files = os.listdir(modpath)
  models = []
  for f in files:
    try:
      if os.path.isdir("%s/%s" % (modpath, f)) == True:
        vlst = []
        fi = open("%s/%s/varnames.txt" % (modpath, f), "r")
        for line in fi:
          vlst.append(line.strip())
        fi.close()
        models.append([f, vlst, joblib.load('%s/%s/scaler.pkl' % (modpath,f)), joblib.load('%s/%s//model.pkl' % (modpath,f))])
      else:
        continue
    except:
      continue
  # Generate dependency list
  i = 0
  while i < len(models):
    swap = False
    for j in range(i, len(models)):
      try:
        models[i][1].index(models[j][0])
        # the model i is dependent on the model j
        # Thus we swap the last model with this one
        swp = models[-1]
        models[-1] = models[i]
        models[i] = swp
        swap = True
        break
      except:
        continue
    if swap == False:
      i+=1
    else:
      i = 0
  return models

def LoadSMILES(fmolsmi):
  mlist = []
  c = 1
  for mol in fmolsmi:
    smiles = []
    molname = []
    fsmi = open(mol, "r")
    for line in fsmi:
      # Try split by tab
      fields = str.split(line.strip(), "\t")
      # Try split by space
      #if len(fields) == 1:
        #fields = str.split(line.strip(), " ")
      #else:
        # Try split by ;
        #if len(fields) == 1:
         # fields = str.split(line.strip(), ";")
        # No molecule name found
        #else:
         # fields = [line.strip(), ("Object %d" % (c))]

      if "[*]" in fields[0]:
        # skip
        continue
      else:
        if "." in fields[0]:
          fields_ = str.split(fields[0], ".")
          i = 0
          for j in range(len(fields_)):
            if len(fields_[j]) > len(fields_[i]):
              i = j
            else:
              continue
          smiles.append(fields_[i])
        else:
          smiles.append(fields[0])

        molname_ = ""
        for i in range(1, len(fields)-1):
          molname_ += fields[i] + " "
        molname_ += fields[-1]
        molname.append(molname_)

    fsmi.close()

    for i in range(len(smiles)):
      m = Chem.MolFromSmiles(smiles[i])
      if m:
        m.SetProp("_Name", molname[i])
        if m is not None:
          mlist.append(m)
        else:
          print("\n>>>> Error 1 while processing %s %s <<<<\n" % (smiles[i], molname[i]))
      else:
        print("\n>>>> Error 2 while processing %s %s <<<<\n" % (smiles[i], molname[i]))
  return mlist, smiles, molname

def MakeDescriptors(mlist, modpath):
  nms=[x[0] for x in Descriptors._descList]
  models = LoadModels(modpath)
  calc = MoleculeDescriptors.MolecularDescriptorCalculator(nms)
  matrixdescriptors = []
  rowname= ["Molecule"]
  for item in nms:
    rowname.append(str(item))
  modnames = []

  for i in range(len(models)):
    rowname.append(models[i][0])
    modnames.append(models[i][0])

  mdesc = []
  for i in range(len(models)):
    mdesc.append(list())
    for j in range(len(models[i][1])):
      try:
        mdesc[-1].append(rowname.index(models[i][1][j]))
      except:
        print("Descriptor %s not found " % (models[i][1][j]))
        continue
  matrixdescriptors.append(rowname)
  for m in mlist:
    desc = []
    mname = m.GetProp('_Name')
    #print("Processing %s ..." % (mname))
    desc.append(mname)
    descrs = calc.CalcDescriptors(m)
    for x in range(len(descrs)):
      desc.append(str(descrs[x]))

    for i in range(len(models)):
      nvars = None
      mtype = str(models[i][3].__class__)
      if mtype == "<class 'sklearn.svm.classes.SVR'>" or mtype == "<class 'sklearn.svm.classes.SVC'>":
        nvars = models[i][3].shape_fit_[1]
      elif mtype == "<class 'sklearn.ensemble.weight_boosting.AdaBoostClassifier'>":
        nvars = len(models[i][3].feature_importances_)
      else:
        print("MODEL ERROR!! ", mtype)

      #print len(mdesc[i]), nvars, len(desc)

      if len(mdesc[i]) == nvars:
        scaler = models[i][2]
        clf = models[i][3]
        xdesc = []
        for j in mdesc[i]:
          if rowname[j] in modnames and j > nvars:
            j_ = j-2
            xdesc.append(float(desc[j_]))
          else:
            xdesc.append(float(desc[j]))
        try:
          desc.append(clf.predict(scaler.transform(np.array(xdesc).reshape(1, -1)))[-1])
        except:
          desc.append("9999.")
      else:
        desc.append("9999.")
    matrixdescriptors.append(desc)

  return matrixdescriptors

def GetDescValues(fsmi, modnames, modpath, tabres):
  """
  Get the model values for the molecule(s) in a smi file or list of smi file
  """
  mlist, smiles, molname = LoadSMILES([fsmi])
  desc = MakeDescriptors(mlist, modpath)
  header = desc[0]
  resid = []

  for i in range(len(modnames)):
    try:
      indx = desc[0].index(modnames[i])
      resid.append(indx)
    except:
      print("Model %s not found!" % (modnames[i]))

  for i in range(1, len(desc)):
    row = [molname[i-1], smiles[i-1]]
    for j in range(len(resid)):
      row.append(desc[i][resid[j]])
    tabres.append(row)

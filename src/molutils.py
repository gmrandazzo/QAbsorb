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
import numpy as np

from modelutils import *

def GenQPixmapDepiction(smi, height=200, width=200):
  molecule = Chem.MolFromSmiles(smi)
  AllChem.Compute2DCoords(molecule)
  return Draw.MolToQPixmap(molecule, (height,width))

def GenDepiction(smi,  molname, path, height=200, width=200):
  molecule = Chem.MolFromSmiles(smi)
  AllChem.Compute2DCoords(molecule)
  imgname = molname+".png"
  imgname = imgname.replace(' ','_') # replace space with _
  imgname = imgname.replace('/','-') # replace / with -
  imagePath = path+"/"+imgname
  Draw.MolToFile(molecule, imagePath, (height,width))
  return imagePath

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
  rowname = ["Molecule"]
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
      nvars =  len(models[i][1])
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
          if models[i][0] == "m3":
              hdmvalue = float(clf.predict(scaler.transform(np.array(xdesc).reshape(1, -1)))[-1])
              if hdmvalue <= 0.5: #Permeable
                desc.append(0)
              else: #Not permeable
                desc.append(1)
          else:
              desc.append(float(clf.predict(scaler.transform(np.array(xdesc).reshape(1, -1)))[-1]))
        except:
          desc.append(9999.)
      else:
        desc.append(9999.)
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
      indx = header.index(modnames[i])
      resid.append(indx)
    except:
      print("Model %s not found!" % (modnames[i]))

  for i in range(1, len(desc)):
    row = [molname[i-1], smiles[i-1]]
    for j in range(len(resid)):
      if desc[0][resid[j]] == "m3":
        if desc[i][resid[j]] == 0:
            row.append("GIT+")
        else:
            row.append("GIT-")
      else:
        row.append(desc[i][resid[j]])
    tabres.append(row)

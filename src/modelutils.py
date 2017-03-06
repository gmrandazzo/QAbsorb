'''
@package modelutils.py

modelutils.py was writen by Giuseppe Marco Randazzo <gmrandazzo@gmail.com>
and is distributed under LGPL version 3

Geneve November 2016
'''

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

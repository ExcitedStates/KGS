#! /usr/bin/python
# Copyright (c) 2010 Robert L. Campbell (rlc1@queensu.ca)
# Copyright (c) 2016 Dominik Budday, dominik.budday@fau.de
# this version inspired by original version of Robert L. Campbell

import math,sys
import os
from pymol import cmd, stored

def distx2(x1,x2):
  """
  Calculate the square of the distance between two coordinates.
  Returns a float
  """
  distx2 = (x1[0]-x2[0])**2 + (x1[1]-x2[1])**2 + (x1[2]-x2[2])**2
  return distx2


def rmsf_states(selection,byres=1,referenceState=0,debug=0,fileEnding=None,selection2=None):
  """
  AUTHOR

    Dominik Budday

  USAGE

    rmsf_states selection ,[byres=1], [referenceState=0], [debug=0], [selection2=None]

    Calculate the RMS fluctuation for each residue in a set of models
    (within a multi-state object).

    If byres=1 (default), then it will calculate the RMSF on a per residue basis and will modify the B-factors for
    all atoms in the residue to be the same, mean value.

    By default, the RMSF is calculated with respect to the mean coordinate (referenceState = 0)
    but you can specify a desired state by its integer number.

    If selection2 is given, then first both selections are computed separately, and then averaged according to the number of states present.

    Please align/intra_fit your selection previously to calling this script.
  """

# input params
  byres = int(byres)
  referenceState = int(referenceState)
  debug = int(debug)

  selections=[]
  selections.append(selection)
  stateNums=[]

  if selection2 != None:
    selections.append(selection2)

  count=0
  for selection in selections:
    # initialize necessary arrays
    models = []
    coordArray = []
    chain=[]
    resi=[]
    resn=[]
    name=[]
    numStates = cmd.count_states(selection)
    stateNums.append(numStates)

    #Debugging info
    if debug:
      print "There are %d states in your selection" % numStates

    # get models into models array
    for i in range(numStates):
      models.append(cmd.get_model(selection,state=i+1))

    numAtoms = len(models[0].atom)
    
    # extract coordinates and atom identification information out of the models
    # loop over the states
    for i in range(numStates):
      coordArray.append([])
      chain.append([])
      resi.append([])
      resn.append([])
      name.append([])

      numAtomsI = len(models[i].atom)
      if numAtomsI != numAtoms:
        print "States have  different atom numbers, please use rmsd_states_diffAtoms"
        sys.exit()

      for j in range(len(models[i].atom)):
          if debug:
            print "i, j:", i,j
          atom = models[i].atom[j]
          coordArray[i].append(atom.coord)
          chain[i].append(atom.chain)
          resi[i].append(atom.resi)
          resn[i].append(atom.resn)
          name[i].append(atom.name)


    #Reference state
    meanArray = [[0 for i in range (3)] for j in range(0,numAtoms)]
    if referenceState==0:
      #Calculate the mean positions for use as reference
      for j in range(0,numAtoms): #len(models[0].atom)):
        for k in range(3):
          for i in range(0,numStates):
            meanArray[j][k]+=coordArray[i][j][k]
          meanArray[j][k] = meanArray[j][k]/float(numStates)
    else:
      for j in range(0,numAtoms): #len(models[0].atom)):
        for k in range(3):
          meanArray[j][k] = coordArray[referenceState - 1][j][k]

    # initialize array
    diff2 = []
    # calculate the root mean squared distance fluctuation wrt the mean or reference state
    if byres:#mean value per residue
      currentRes=resi[0][0]
      startInd=0
      atomCount=0
      for j in range(len(coordArray[0])):
        diff2.append(0.)
        atomCount +=1
        for i in range(numStates):
          diff2[j] += distx2(coordArray[i][j],meanArray[j])

        if j==len(coordArray[0])-1: # end of all coordinates 
          #Set the whole residue to mean value
          #Mean over residue
          sumRes = sum(diff2[startInd:])/(atomCount) # residue mean squared fluctuations over all states
          #Mean over states
          for k in range(startInd,j+1):
            diff2[k] = math.sqrt(sumRes/numStates)  # root mean square fluctuations 
          #print "Set %g for %d atoms in residue %s" % (math.sqrt(sumRes/numStates), j-startInd,currentRes)
          break
        else:
          if resi[0][j+1]!=currentRes:#Reached the end of the current residue
            #Set the whole residue to mean value
            #Mean over residue
            sumRes = sum(diff2[startInd:j+1])/(atomCount)
            #Mean over states
            for k in range(startInd,j+1):
              diff2[k] = math.sqrt(sumRes/numStates) 
            #print "Set %g for %d atoms in residue %s" % (math.sqrt(sumRes/numStates), j-startInd,currentRes)
            startInd=j+1
            atomCount=0
            currentRes=resi[0][j+1]         


    else: #individual value for each atom
      for j in range(len(coordArray[0])):
        diff2.append(0.)
        for i in range(numStates):
          diff2[j] += distx2(coordArray[i][j],meanArray[j])

        # divide the fluctuation squared by the number of states and take the square root of that to get the RMSF
        diff2[j] = math.sqrt(diff2[j]/numStates)

    # Now we are done with the computation
    if len(selections) > 1 and count==0:
      diff1 = diff2
      count = 1

  #Compute averaged diff list
  sumStates = sum(stateNums)
  print "Overall "+str(sumStates)+" states in all selections"

  # Combine overall in diff2
  if len(selections) > 1:
    for j in range(len(coordArray[0])):
      diff2[j] = diff2[j]*stateNums[1]/sumStates + diff1[j]*stateNums[0]/sumStates

  for selection in selections:
    # reset all the B-factors to zer
    cmd.alter(selection,"b=0")

    b_dict = {}

    sum_rmsf = 0
    
    # orig_stdout = sys.stdout
    # if fileEnding != None:
    #   fout = file('RMSF_',fileEnding,'.txt', 'w')
    # else:
    #   fout = file('RMSF.txt', 'w')
    # 
    # sys.stdout = fout
    
    if byres: #alter all B-factors in a residue to be the same
      def b_lookup(chain, resi, name):
        if chain in b_dict:
          b = b_dict[chain][resi]
        else:
          b = b_dict[''][resi]
        return b
      stored.b = b_lookup
      
      #Print RMSF table on screen
      prevRes=-1
      for i in range(len(diff2)):
        if resi[0][i] != prevRes:
          print resi[0][i], diff2[i]
        prevRes=resi[0][i]

      # Alter the B-factor according to B = 8*pi^2*1/3*<u^2> where <u^2> is the mean squared displacement (diff2^2) 
      for i in range(len(diff2)):
        sum_rmsf += diff2[i]
        #print "Current sum is %g, added value %g" % (sum_rmsf, diff2[i])
        b_dict.setdefault(chain[0][i], {})[resi[0][i]] = diff2[i]**2*8/3*math.pi**2
      mean_rmsf = sum_rmsf/len(diff2)
      #print "Value is %g for %d atoms." % (sum_rmsf, len(diff2))
      cmd.alter(selection,'%s=stored.b(chain,resi,name)' % ('b'))

    else:#alter all B-factors individually according to the RMSF for each atom
      quiet=0
      def b_lookup(chain, resi, name):
        if chain in b_dict:
          b = b_dict[chain][resi][name]
        else:
          b = b_dict[''][resi][name]
        return b
      stored.b = b_lookup

      #Print table
      for i in range(len(diff2)):
        print resi[0][i],name[0][i], diff2[i]

      for i in range(len(diff2)):
        sum_rmsf += diff2[i]
        try:
          b_dict.setdefault(chain[0][i], {}).setdefault(resi[0][i], {})[name[0][i]] = diff2[i]**2*8/3*math.pi**2
        except KeyError:
          print chain[0][i],resi[0][i],name[0][i]
      cmd.alter(selection,'%s=stored.b(chain,resi,name)' % ('b'))
      mean_rmsf = sum_rmsf/len(diff2)

    # sys.stdout = orig_stdout
    # fout.close()
      
    print "Mean RMSF for selection: %s = %g" % (selection, mean_rmsf)

cmd.extend("rmsf_states",rmsf_states)                                                                                                                                                                                                                                                                                                                                                                                              
#!/usr/bin/env python
'''
Convert DFTB parameters in skf format into LATTE's table format

Usage
-----

```bash
DLtab.py <skf folder> <LATTE parameter folder>
```

Three files will be generate for LATTE, electrons.dat, bondints.table and ppots.dftb.

[note for electron.dat file]
 * The spin parameter is not included in the skf file, thus, they will be 
set to zero, if spin polarized DFTB is needed, don't forget to modify
the electrons.dat file

 * The masses defined in skf file is different than the one defined in
changecoord code. Thus, if you need run LATTE-LAMMPS code, don't forget
to modify the masses .
 * The script uses Hubbard U for the s shell as U parameter for the 
 element, which is the same as the default for DFTB+ code.
''' 

__license__ = "GPL V2"
__author__ = "Chang Liu"
__email__ = "cliu18@ncsu.edu"
__version__ = 1.01

import os
import sys
import glob
import numpy as np

# Constants for atomic unit <-> ev, angstroms
bohr2angstroms = 0.529177249
hatree2ev = 27.21138602
# strconv
def index2str(indexint):
    #  convert interaction string such as "sss" to index of matrix in dftb parameter
    return{
        0: 'dds',
        1: 'ddp',
        2: 'ddd',
        3: 'pds',
        4: 'pdp',
        5: 'pps',
        6: 'ppp',
        7: 'sds',
        8: 'sps',
        9: 'sss',
    }[indexint]
def str2index(strint):
    #  convert interaction string such as "sss" to index of matrix in dftb parameter
    return{
        'dds': 0,
        'ddp': 1,
        'ddd': 2,
        'pds': 3,
        'pdp': 4,
        'pps': 5,
        'ppp': 6,
        'sds': 7,
        'sps': 8,
        'sss': 9,
    }[strint]
# define the base for H and s 
baseH = 0
baseS = 10

def not_printed_bi(currType, printedPair):
    # decide bi printed or not
    items = currType.split()
    atom1 = items[0]
    atom2 = items[1]
    bitype = items[2]
    for currPair in printedPair:
        items = currPair.split()
        tmpatom1 = items[0]
        tmpatom2 = items[1]
        tmpbitype = items[2]
        # only when atom1 and atom2 swap, it is possible redundant
        if (atom1 == tmpatom2) and (atom2 == tmpatom1) \
          and (bitype == tmpbitype) and (bitype[0] == bitype[1]):
            return False
    return True

def DLconv_ee(atom,line1,line2):
    # get element informations
    Element = atom
    Mass = float(line2.replace(',',' ').split()[0])
    # this part will need to get from somewhere else
    Wss = 0
    Wpp = 0
    Wdd = 0
    Wff = 0
    # other are from line1
    Ef = 0 # till now, no f elements
    (Ed,Ep,Es,SPE,Ud,Up,Us,fd,fp,fs) = map(float,line1.replace(',',' ').split())
    Ed = Ed * hatree2ev
    Ep = Ep * hatree2ev
    Es = Es * hatree2ev
    Numel = fd + fp + fs
    # to reproduce the DFTB+ method, HubbardU will be the value of s shell
    HubbardU = Us * hatree2ev
    # spin coefficient are also pre-computed 
    if fd != 0:
        basis = 'spd'
    elif fp != 0:
        basis = 'sp'
    else:
        basis = 's'
    # unpacked Done
    currEE = '%s %s %.1f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f' \
        % (Element,basis,Numel,Es,Ep,Ed,Ef,Mass,HubbardU,Wss,Wpp,Wdd,Wff)
    return currEE

def DLconv_HS_tab(strDHS,gridDistB,nPts,atom1, atom2):
    # convert to LATTE format
    # alloc the matrix 
    MdHS = np.zeros((20,nPts))
    # convert gird distant to angstroms and generate grid points array
    gridDistA = gridDistB * bohr2angstroms
    gridPts = np.linspace(gridDistA, gridDistA*nPts, num=nPts)
    # filling the matrix
    iPts = 0
    for line in strDHS:
        iOrbs = 0
        for item in line.split():
            tmpList = item.split('*')
            if len(tmpList) == 1:
                MdHS[iOrbs][iPts] = float(tmpList[0].split(',')[0])
                iOrbs += 1
            else:
                for _ in range(int(tmpList[0])):
                    MdHS[iOrbs][iPts] = float(tmpList[1].split(',')[0])
                    iOrbs += 1
        iPts += 1
    # generate the cell for print
    BIlist = [gridPts] # only 1 distance info
    for iOrbs in range(10):
        # strip the leading 1s (which is typical to see in most skf file)
        for istart in range(len(gridPts)):
            if (MdHS[baseS+iOrbs][istart] != 1.) and (MdHS[baseH+iOrbs][istart] != 1.):
                break
        # get prepared for information needed
        intStr = index2str(iOrbs) 
        tmpH = MdHS[baseH+iOrbs] * hatree2ev
        tmpS = MdHS[baseS+iOrbs]
        # determine whether parameter exist or not
        if (np.linalg.norm(tmpH[istart:]) * np.linalg.norm(tmpS[istart:]) == 0):
            continue
        else:
            bitypeStr = "%s %s %s" % (atom1, atom2, intStr)
            BIlist.append([bitypeStr,tmpS,tmpH, istart])
    return BIlist

def compute_rep_preexp(expStr, endR, dR):
    """compute the exponential part of repulsive potential"""
    (a1, a2, a3) = map(float, expStr.split())
    npt = int(np.floor(endR/dR) + 1)
    preRptsB = np.zeros(npt)
    preEptsH = np.zeros(npt)
    for i in range(npt):
        preRptsB[i] = (i+1) * dR
        preEptsH[i] = np.exp(-a1*preRptsB[i] + a2) + a3
    return (preRptsB,preEptsH)


def DLconv_repp(atom1, atom2, strDrp,latteFD):
    # convert ppots to latte format
    # read informations from dftb parameters
    strDrp.pop(0) # skip the header "Spline"
    tmplist = strDrp.pop(0).split()
    nDPts = int(tmplist[0]) # number of spline domain in dftb parameters
    RcutB = float(tmplist[1]) # cutoff radii in bohr 
    RptsB = np.zeros(nDPts+1) # points radii in Bohr
    EptsH = np.zeros(nDPts+1) # list to save energy in Hartree
    iPts = 0
    shortEXP = strDrp.pop(0) # the short distant exponential coefficients
    # direct use the gird point from dftb pp list
    for tmpline in strDrp:
        if len(tmpline.split()) != 6:
            break
        (r0, _, c0, _, _, _) = map(float,tmpline.split())
        EptsH[iPts] = c0
        RptsB[iPts] = r0
        iPts += 1
    # the last part
    (r0, _, c0, _, _, _, _, _) = map(float,tmpline.split())
    EptsH[iPts] = c0
    RptsB[iPts] = r0
    iPts += 1
    EptsH[iPts] = 0
    RptsB[iPts] = RcutB
    # compute the beginning exp part
    (preRptsB, preEptsH) = compute_rep_preexp(shortEXP, RptsB[0], RptsB[1]-RptsB[0])
    # convert atomic unit to ev and angstroms
    RptsA = RptsB * bohr2angstroms
    EptsEv = EptsH * hatree2ev
    preRptsA = preRptsB * bohr2angstroms
    preEptsEv = preEptsH * hatree2ev
    # return the final string
    strLATTE = "%s %s \n %d\n" % (atom1, atom2, nDPts+1+len(preRptsB)) # header for LATTE block
    for i in range(len(preRptsB)):
        strLATTE += "%.15E %.15E \n" % (preRptsA[i], preEptsEv[i])
    for i in range(nDPts+1):
        strLATTE += "%.15E %.15E \n" % (RptsA[i], EptsEv[i])
    return strLATTE

def doDLconvert(skfFile,latteFD):
    # get atom names
    atomPair = skfFile.replace('/','.').split('.')[-2]
    #print atomPair
    (atom1,atom2) = atomPair.split('-')
    # read the file 
    fp = open(skfFile)
    strList = fp.read().splitlines()
    fp.close
    # get HS info and atomic info(if atom1 == atom2)
    header = strList.pop(0).replace(',',' ').split()
    gridDistB = float(header[0])
    #nPts = int(header[1]) in principle.. However it's wrong
    if atom1 == atom2:
        line1 = strList.pop(0)
        line2 = strList.pop(0)
        currEE = DLconv_ee(atom1,line1,line2)
    else:
        Drepfun = strList.pop(0) # save it just in case no spline
        currEE = ''
    strLHS = []
    while strList:
        line = strList.pop(0)
        if 'Spline' in line: 
            # this is how the two section split, nPts might be wrong
            break
        strLHS.append(line)
    nPts = len(strLHS)
    BItable = DLconv_HS_tab(strLHS,gridDistB,nPts,atom1,atom2)
    # get PP info        
    strDrp = [line]
    while strList:
        line = strList.pop(0)
        if ('Documentation' in line) or (len(line.strip())==0):
            break
        strDrp.append(line)
    currPP = DLconv_repp(atom1, atom2, strDrp,latteFD)
    # Done!
    return (currEE,BItable,currPP)

def printPP(ppots):
    fp = open('ppots.dftb','w')
    header = '%d \n' % len(ppots)
    fp.write(header)
    for line in ppots:
        fp.write(line)
    fp.close()

def printEE(electrons):
    fp = open('electrons.dat','w')
    header1 = 'Noelem= %d \n' % len(set(filter(None, electrons)))
    fp.write(header1)
    header2 = 'Element basis Numel Es Ep Ed Ef Mass HubbardU Wss Wpp Wdd Wff \n'
    fp.write(header2)
    for line in electrons:
        if line:
            fp.write(line+' \n')
    fp.close()

def printBItab(biTab,nbiTab):
    BIstr = []
    printedPair = []
    for biList in biTab:
        distList = biList[0]
        nline = len(distList)
        for currtypeBI in biList[1:]:
            if not_printed_bi(currtypeBI[0], printedPair):
                istart = currtypeBI[-1]
                BIstr.append('%s\n' % currtypeBI[0])
                BIstr.append('%d\n' % (nline-istart))
                for i in range(nline-istart):
                    BIstr.append('%12.5f %15e %15e \n' % \
                    (distList[i+istart], currtypeBI[1][i+istart], currtypeBI[2][i+istart]))
                printedPair.append(currtypeBI[0])
            else:
                nbiTab -= 1
    # done, write the file
    fp = open('bondints.table','w')
    fp.write('Noints= %d\n' % nbiTab)
    for line in BIstr:
        fp.write(line)
    fp.close()


if __name__=='__main__':
    # record dir and get path
    # CWD dftbFD latteFD
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    CWD = os.getcwd() # record curr dir
    dftbFD = sys.argv[1]
    if not(os.path.exists(dftbFD) and os.path.isdir(dftbFD)):
        sys.exit("The dftb folder not exists!")
    elif not(glob.glob(dftbFD+'/*.skf')):
        sys.exit("No skf file in dftb folder!")
    latteFD = sys.argv[2]
    if not(os.path.exists(latteFD)):
        os.mkdir(latteFD)
    elif not(os.path.isdir(dftbFD)):
        sys.exit("Cannot create LATTE folder!")
    # get skf list and loop over it to generate the string lists
    skfList = glob.glob(dftbFD+'/*.skf')
    bondint = []
    biTab = []
    nbiTab = 0
    ppots = []
    electrons = []
    existPP = []
    existBI = []
    for skfFile in skfList:
        print("reading %s" % skfFile)
        (atom1,atom2) = (skfFile.replace('/','.').split('.')[-2]).split('-')
        (currEE,currBItab,currPP) = doDLconvert(skfFile,latteFD)
        biTab.append(currBItab)
        nbiTab += (len(currBItab) - 1)
        if not ([atom2,atom1] in existPP):
            ppots.append(currPP)
            existPP.append([atom1,atom2])
        electrons.append(currEE)
    # do the output
    os.chdir(latteFD)
    printPP(ppots)
    printEE(electrons)
    printBItab(biTab,nbiTab)
    # Done
    os.chdir(CWD)
    exit(0)
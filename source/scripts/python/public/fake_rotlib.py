import argparse
import os
import sys
from multiprocessing import Pool
from subprocess import check_output
from shutil import move
from copy import deepcopy

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit import rdBase as rdb
except:
    print("RDKit must be installed to run this script.  Please install RDKit (version 2022.09.5 preferred) into your python environment.")
    exit(1)

#Argument parser for terminal input
def parseArgs(argv):
    parser = argparse.ArgumentParser(description="This script uses RDKit and molfile_to_params_polymer.py to "
                                     "quickly parameterize a NCAA for use in Rosetta")
    parser.add_argument("-i","--input",required=True,metavar="FILE",
                        help="Input sdf mol file (REQUIRED)")
    parser.add_argument("-n","--top_n_confs",type=int,metavar="N",
                        default = 1000,
                        help="Top number of generated conformations to use as rotamers (default: %(default)s)")
    parser.add_argument("-j","--rdkit_max_iter",type=int,metavar="ITER",
                        default = 1000000,
                        help="max iterations for RDKit geometry optimization (default: %(default)s)")
    parser.add_argument("-c","--n_cbb",type=int,metavar="N",
                        default = 1,
                        help="number of carbons in the backbone, i.e., N = 2 for beta AA, 3 for gamma AA, etc. (default: %(default)s)")
    parser.add_argument("-p","--n_processes",type=int,metavar="N",
                        default=1,
                        help="Number of processes to use for conformer generation/scoring (default: %(default)s)")
    parser.add_argument("-m","--mmff",action="store_true",
                        default=False,
                        help="If specified, MMFF94 will be used for geometry optimization and conformer scoring instead of UFF.")
    parser.add_argument("-r","--rotlib",action="store_true",
                        default=False,
                        help="If specified, a rotlib file will be used instead of pdb rotamers (requires numpy and scikit-learn)")
    parser.add_argument("-d","--dip",action="store_true",
                        default=False,
                        help="If specified, the input sdf is assumed to be a dipeptide (i.e. no capping will occur)")
    parser.add_argument("-o","--no_rot",action="store_true",
                        default=False,
                        help="If specified, no geometry optimization or conformer generation will occur. MUST be used with --dip.")
    return parser.parse_args()

#Write atom assignments in molecule SDF for m2pp to read
def generateInstructions(dipmol, nCbb):
    #This whole thing hinges on the assumption that the identified backbone atoms will be in the same order as the smiles string
    #If this is no longer the case, may God help us all
    instructions = "\n> <RosettaParamsInstructions>\n"
    cbb = "C"*nCbb
    backbone = dipmol.GetSubstructMatch(Chem.MolFromSmiles("CNC(=O)%sNC(=O)C"%cbb))
    ignore = [backbone[0], backbone[6+nCbb], backbone[7+nCbb]]

    #I wouldn't have to do this if RDKit would make bonds consistently bidirectional...
    upperBond = [b.GetEndAtomIdx() if b.GetEndAtomIdx() != backbone[1] else b.GetBeginAtomIdx()
                 for b in dipmol.GetAtomWithIdx(backbone[1]).GetBonds()]
    for idx in upperBond:
        if idx != backbone[2] and idx != backbone[0] and idx != backbone[1]:
            ignore.append(idx)

    upperMetBond = [b.GetEndAtomIdx() if b.GetEndAtomIdx() != backbone[0] else b.GetBeginAtomIdx()
                    for b in dipmol.GetAtomWithIdx(backbone[0]).GetBonds()]
    for idx in upperMetBond:
        if idx != backbone[1] and idx != backbone[0]:
            ignore.append(idx)

    lowerMetBond = [b.GetEndAtomIdx() if b.GetEndAtomIdx() != backbone[7+nCbb] else b.GetBeginAtomIdx()
                    for b in dipmol.GetAtomWithIdx(backbone[7+nCbb]).GetBonds()]
    for idx in lowerMetBond:
        if idx != backbone[5+nCbb] and idx != backbone[7+nCbb]:
            ignore.append(idx)

    properties = "PROTEIN POLYMER"
    
    if nCbb == 1:
        properties += " ALPHA_AA"
    elif nCbb == 2:
        properties += " BETA_AA"
    elif nCbb == 3:
        properties += " GAMMA_AA"
    
    #This gets called twice in the script because the geometry optimizer kekulizes
    Chem.SetConjugation(dipmol)
    Chem.SetAromaticity(dipmol)
    if any([b.GetIsAromatic() for b in dipmol.GetBonds()]):
        properties += " AROMATIC"
    else:
        properties += " ALIPHATIC"

    #Note: this doesn't work for cysteine because of the difference between R/S and L/R
    chirals = Chem.FindMolChiralCenters(dipmol)
    for c in chirals:
        if c[0] == backbone[3+nCbb]:
            if c[1] == "S":
                properties += " L_AA"
            elif c[1] == "R":
                properties += " D_AA"
    
    if Chem.GetFormalCharge(dipmol) > 0:
        properties += " CHARGED POSITIVE_CHARGE"
    elif Chem.GetFormalCharge(dipmol) < 0:
        properties += " CHARGED NEGATIVE_CHARGE"

    if Chem.rdMolDescriptors.CalcNumRings(dipmol) > 0:
        properties += " CYCLIC"

    instructions += "M  ROOT %d\nM  POLY_N_BB %d\n"%(backbone[4+nCbb]+1,backbone[4+nCbb]+1)
    instructions += "M  POLY_CA_BB %d\n"%(backbone[3+nCbb]+1)
    instructions += "M  POLY_C_BB %d\n"%(backbone[2]+1)
    instructions += "M  POLY_O_BB %d\n"%(backbone[3]+1)
    instructions += "M  POLY_IGNORE "+" ".join([str(x+1) for x in ignore])+"\n"
    instructions += "M  POLY_UPPER %d\n"%(backbone[1]+1)
    instructions += "M  POLY_LOWER %d\n"%(backbone[5+nCbb]+1)
    instructions += "M  POLY_CHG %d\n"%Chem.GetFormalCharge(dipmol)
    instructions += "M  POLY_PROPERTIES "+properties+"\n"
    instructions += "M  END\n\n$$$$\n"
    return instructions

#Attach dipeptide caps onto the NCAA
def generateDip(ncaamol, nCbb):
    #Cap all N and find the results that look like capped backbone
    cbb = "C"*nCbb
    backbone = ncaamol.GetSubstructMatches(Chem.MolFromSmiles("N%sC(=O)O"%cbb))
    if len(backbone) > 1:
        print("ERROR: More than one substructure looks like the backbone.")
        exit(2)
    elif len(backbone) == 0:
        print("ERROR: Backbone not found in provided substructure. If your NCAA is not an alpha AA, please set the -c flag.")
        exit(1)
    backinds = backbone[0]
    ncent = ncaamol.GetConformer().GetAtomPosition(backinds[0])
    #Make N and C terminal cap, move them in position
    ncap = Chem.MolFromSmiles("C(=O)C")
    rdb.DisableLog("rdApp.warning") #Yes, RDKit, I know they don't have H's
    AllChem.EmbedMolecule(ncap)
    rdb.EnableLog("rdApp.warning")
    for idx in range(ncap.GetConformer().GetNumAtoms()):
        ncap.GetConformer().SetAtomPosition(idx,ncent)

    ccent = ncaamol.GetConformer().GetAtomPosition(backinds[-1])
    ccap = Chem.MolFromSmiles("NC")
    rdb.DisableLog("rdApp.warning")
    AllChem.EmbedMolecule(ccap)
    rdb.EnableLog("rdApp.warning")
    for idx in range(ccap.GetConformer().GetNumAtoms()):
        ccap.GetConformer().SetAtomPosition(idx,ccent)

    #Clean up N terminus to get ready to attach
    nheavyb = [b for b in ncaamol.GetAtomWithIdx(backinds[0]).GetBonds() 
               if b.GetBeginAtomIdx() != b.GetEndAtomIdx() and b.GetEndAtom().GetAtomicNum() != 1 
               and b.GetBeginAtom().GetAtomicNum() != 1]
    if len(nheavyb) > 2:
        ncaamol.GetAtomWithIdx(backinds[0]).SetFormalCharge(len(nheavyb)-2)
        ncaamol.GetAtomWithIdx(backinds[0]).SetNumExplicitHs(0)
    else:
        ncaamol.GetAtomWithIdx(backinds[0]).SetFormalCharge(0)
        ncaamol.GetAtomWithIdx(backinds[0]).SetNumExplicitHs(2-len(nheavyb))
    #Combine caps & ncaa, draw bonds, remove charged O from Cbb
    combo = Chem.CombineMols(ncap,Chem.CombineMols(ccap,ncaamol))
    editcombo = Chem.EditableMol(combo)
    editcombo.AddBond(0,
                      backinds[0]+len(ncap.GetAtoms())+len(ccap.GetAtoms()),
                      order = Chem.rdchem.BondType.SINGLE)
    editcombo.AddBond(len(ncap.GetAtoms()),
                      backinds[1+nCbb]+len(ncap.GetAtoms())+len(ccap.GetAtoms()),
                      order = Chem.rdchem.BondType.SINGLE)
    editcombo.RemoveAtom(backinds[-1]+len(ncap.GetAtoms())+len(ccap.GetAtoms()))

    dipmol = editcombo.GetMol()
    #Loadbearing optimization, no idea why this needs to be here for protonation to work
    rdb.DisableLog("rdApp.warning") #Again... I know it doesn't have H's
    geomOptimize(dipmol, 1000)
    rdb.EnableLog("rdApp.warning")
    dipmol = Chem.AddHs(dipmol)
    instructions = generateInstructions(dipmol, nCbb)
    return (dipmol, instructions)

#Optimize geometry of NCAA
def geomOptimize(mol, maxIter, mmff=False):
    AllChem.MMFFSanitizeMolecule(mol)
    Chem.AssignStereochemistryFrom3D(mol)
    AllChem.EmbedMolecule(mol)
    if mmff:
        optflag = AllChem.MMFFOptimizeMolecule(mol, maxIters=maxIter)
    else:
        optflag = AllChem.UFFOptimizeMolecule(mol, maxIters=maxIter)
    if optflag == 1:
        print("NCAA requires more iterations to optimize. Please increase the max iteration count and retry")
    elif optflag == -1:
        print("Error with RDKit geometry optimization, skipping optimization")
        return mol

#Generate conformations of NCAA for rotamer library
def rdkitConf(ncaamol, noext, topN, instructions, mmff, nProc):
    print("Generating RDKit Conformations...")
    confs = AllChem.EmbedMultipleConfs(ncaamol, useRandomCoords=True, numConfs=int(10*topN), numThreads=0, useExpTorsionAnglePrefs=False)
    print("Aligning conformations...")
    AllChem.AlignMolConformers(ncaamol, confIds=confs)
    print("Scoring conformations...")
    lines = instructions.split("\n")
    polylines = []
    for line in lines:
        if "POLY_IGNORE" in line or "POLY_UPPER" in line or "POLY_LOWER" in line:
            polylines.append(line)
    ignoreidx = []
    for line in polylines:
        ignoreidx += [int(i) for i in line.split()[2:]]
    ignoreidx.sort(reverse=True)
    editmol = Chem.EditableMol(ncaamol)
    for idx in ignoreidx:
        editmol.RemoveAtom(idx-1)

    #Score the poses
    stripmol = Chem.AddHs(editmol.GetMol())
    AllChem.MMFFSanitizeMolecule(stripmol)
    p = Pool(nProc)
    argList = [(stripmol, c, mmff) for c in confs]
    scoreconfs = p.starmap(scoreConf, argList)
    scoreconfs.sort(key=lambda x: x[1])
    Chem.SetConjugation(ncaamol)
    Chem.SetAromaticity(ncaamol)
    g = open(noext+"_rotamer.sdf","w")
    print("Writing conformations...")
    for i in range(topN):
        g.write(Chem.MolToMolBlock(ncaamol,confId=scoreconfs[i][0],kekulize=False))
        g.write(instructions)
    g.close()

#Score all generated conformations and return the most favorable
def scoreConf(mol, c, mmff):
    #Rip off the dihedral caps (they interfere with the energy of the sidechain)
    if mmff:
        prop = AllChem.MMFFGetMoleculeProperties(mol)
        ff = AllChem.MMFFGetMoleculeForceField(mol,prop,confId=c)
        return (c, ff.CalcEnergy())
    else:
        ff = AllChem.UFFGetMoleculeForceField(mol,confId=c)
        return (c, ff.CalcEnergy())

#Convert between [0,360] angle and [-180,180] angle
def negposang(ang,pos):
    if pos:
        if ang >=0.:
            return ang
        else:
            return 360.+ang
    else:
        if ang > 180.:
            return ang-360.
        else:
            return ang

#Taken from https://stackoverflow.com/questions/20305272/
def dihedral(p):
    """Praxeolitic formula
    1 sqrt, 1 cross product"""
    p0 = p[0]
    p1 = p[1]
    p2 = p[2]
    p3 = p[3]

    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    ang = np.degrees(np.arctan2(y, x))
    return negposang(ang, True)
    
#Initialize trie for rotamer clusters
def buildRotClusts(rotBins,depth):
    rotClusts = [None for i in range(len(rotBins[depth]))]
    if depth+1 != len(rotBins):
        for i in range(len(rotBins[depth])):
            rotClusts[i] = buildRotClusts(rotBins,depth+1)
    return rotClusts

#Assign cluster to each chi bin, if a cluster is already there, merge it
def placeClust(clust,rotClusts,rotBins,depth):
    dists = np.array([np.abs(clust[0][depth]-rb) for rb in rotBins[depth]])
    assignment = dists.argmin()
    if depth+1 == len(rotBins):
        if rotClusts[assignment] != None:
            newweight = clust[2]+rotClusts[assignment][2]
            newmean = clust[2]/newweight*clust[0]+rotClusts[assignment][2]/newweight*rotClusts[assignment][0]
            newstd = np.sqrt(clust[2]/newweight*clust[1]**2+rotClusts[assignment][2]/newweight*rotClusts[assignment][1]**2)
            rotClusts[assignment] = (newmean,newstd,newweight)
        else:
            rotClusts[assignment] = clust
    else:
        placeClust(clust,rotClusts[assignment],rotBins,depth+1)

#Given assigned clusters, retrieve information in preparation for writing
def listBins(binList,rotClusts,rotBins,inds):
    for i in range(len(rotBins[len(inds)])):
        inds.append(i)
        if len(inds) == len(rotBins):
            rotinds = [n+1 for n in inds] + [0 for n in range(len(inds),4)]
            if rotClusts[i] == None:
                rotmeans = [negposang(rotBins[n][inds[n]].item(), False) for n in range(len(inds))] + [0. for n in range(len(inds),4)]
                rotsd = [10. for n in range(len(inds))] + [0. for n in range(len(inds),4)]
                binList.append([0.]+rotinds+rotmeans+rotsd)
            else:
                rotmeans = [negposang(n, False) for n in rotClusts[i][0].tolist()] + [0. for n in range(len(inds),4)]
                rotsd = rotClusts[i][1].tolist() + [0. for n in range(len(inds),4)]
                rotweight = rotClusts[i][2]
                binList.append([rotweight]+rotinds+rotmeans+rotsd)
        else:
            listBins(binList,rotClusts[i],rotBins,inds)
        inds.pop()

#Use SVD to superpose each dihedral center bond to face upward
def standardizeChis(chis, rotamers):
    reference = np.array([[0.,0.,0.],[0.,1.,0.]])
    allxyz = np.empty((0,3*len(chis)))
    for rot in rotamers:
        stdxyz = np.empty((0,3))
        for i in range(len(chis)):
            chixyz = np.array([rot[atom] for atom in chis[i]])
            chixyz -= chixyz[1,:]
            cosangle = np.dot(chixyz[0,:], chixyz[2,:])/(np.linalg.norm(chixyz[0,:])*np.linalg.norm(chixyz[2,:]))
            bottomdot = np.array([-np.sqrt(1-cosangle**2), cosangle, 0.])
            ideal = np.vstack((bottomdot,reference))
            U, S, Vt = np.linalg.svd(np.matmul((chixyz[:3,:]).T, ideal))
            rotmat = np.matmul(Vt.T, U.T)
            if np.linalg.det(rotmat) < 0.:
                Vt[2,:] *= -1
                rotmat = np.matmul(Vt.T, U.T)
            rotxyz = np.matmul(rotmat, chixyz.T)
            stdxyz = np.vstack((stdxyz, rotxyz[:,3]))
        allxyz = np.vstack((allxyz, stdxyz.flatten()))
    return allxyz

#Turn a set of cartesian clusters into dihedral parameter space
def makeDihedral(clusts):
    #The right angle is completely arbitrary, they just have to define the first plane
    ideal = np.array([[-1.,0.,0.],[0.,0.,0.],[0.,1.,0.]])
    newclusts = []
    for c in clusts:
        dihmeans = []
        dihstd = []
        for i in range(int(c[0].shape[0]/3)):
            point = c[0][i*3:i*3+3]
            meanpos = np.vstack((ideal, point))
            dihmeans.append(dihedral(meanpos))
            #To calculate SDs, approximate circular path as linear and find where it intersects with the cartesian SD ellipse
            a = np.sqrt(c[1][i*3]) #turn the variances into sds
            b = np.sqrt(c[1][i*3+2])
            m = -point[0]/point[2]
            x = np.sqrt(a**2*b**2/(b**2+a**2*m**2))
            z = m*x
            stdpos = np.vstack((ideal, point+np.array([x,0.,z])))
            dihstd.append(abs(dihedral(stdpos)-dihmeans[-1]))
        newclusts.append((np.array(dihmeans), np.array(dihstd), c[2]))
    print()
    for c in newclusts:
        print(c)
    return newclusts

#Given the rotamer conformation set, cluster the dihedrals in cartesian space and write a rotlib type file describing the clusters
def fakeRotLib(pdbrot, params):
    #Read atomic positions from rotamer.pdb
    f = open(pdbrot)
    atoms = {}
    rotamers = []
    for line in f:
        if line[:3] == "TER":
            rotamers.append(deepcopy(atoms))
            atoms = {}
        else:
            atoms[line[12:16].strip()] = np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])])
    f.close()
    
    #Read chi angles from params
    f = open(params)
    chis = []
    for line in f:
        if line[:3] == "CHI":
            parts = line.strip().split()
            chis.append(parts[2:])
        #Oops! The last chi was actually a proton chi, get rid of it
        elif line[:10] == "PROTON_CHI":
            chis.pop()
    f.close()
            
    if len(chis) == 0:
        print("Rotlib can't be built, sidechain has no chi's")
        return
    
    #Extract relative xyz for each chi angle
    allxyz = standardizeChis(chis, rotamers)
    #Rotlib format only accepts the first 4 chis :(
    allxyz = allxyz[:,:12]

    #Cluster the relative xyzs
    mm = BayesianGaussianMixture(n_components = min(10**len(chis), allxyz.shape[0]),
                                 n_init=3, max_iter=10000, covariance_type="diag", 
                                 covariance_prior=[1. for i in range(allxyz.shape[1])])
    mm.fit(allxyz)
    clusts = [(mm.means_[i], mm.covariances_[i], mm.weights_[i]) for i in range(len(mm.means_)) if mm.weights_[i] > 0.005]
    
    #Transform the clusters into dihedral space
    clusts = makeDihedral(clusts)
    
    #Put each chi xyz individually into bins
    rotBins = []
    for i in range(int(allxyz.shape[1]/3)):
        binner = BayesianGaussianMixture(n_components=10, n_init=3, max_iter=10000, covariance_type="diag",covariance_prior=[10.,10.,10.])
        binner.fit(allxyz[:,i*3:i*3+3])
        rotc = [(binner.means_[i], binner.covariances_[i], binner.weights_[i]) for i in range(len(binner.means_)) if binner.weights_[i] > 0.005]
        rotc = makeDihedral(rotc)
        rotBins.append(np.concatenate([c[0] for c in rotc]))
    
    #Recursively build a trie (yes, a trie not a tree) storing each possible combination of chi bins
    rotClusts = buildRotClusts(rotBins,0)
    #Recursively assign each cluster to a leaf of the trie
    #If two clusters collide, since the clusters are multidimensional gaussians, the two clusters can be combined through weighted sum
    for c in clusts:
        placeClust(c, rotClusts, rotBins, 0)
    
    #Recursively populate a list of the trie leaves
    binList = []
    listBins(binList, rotClusts, rotBins, [])

    #Write the rotlib
    binList.sort(key=lambda x:x[0], reverse=True)
    g = open(pdbrot+".rotlib","w")
    for phi in range(-170,190,10):
        for psi in range(-170,190,10):
            for b in binList:
                g.write("UNK\t%d\t%d\t9999\t%d %d %d %d %.6f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\n"%
                        (phi, psi, b[1], b[2], b[3], b[4], b[0], b[5], b[6], b[7], b[8], b[9], b[10], b[11], b[12]))
    g.close()
    
    #Modify the params file to utilize the rotlib file
    g = open(params,"a")
    g.write("NCAA_ROTLIB_PATH "+os.path.abspath(pdbrot)+".rotlib\n")
    g.write("NCAA_ROTLIB_NUM_ROTAMER_BINS "+str(len(rotBins))+" "+" ".join([str(len(r)) for r in rotBins])+"\n")
    g.close()

if __name__ == "__main__":
    args = parseArgs(sys.argv[1:])
    noext = os.path.splitext(args.input)[0]
    ncaaName = os.path.basename(noext)[:3].upper()

    #Generate dipeptide form and RosettaParams instructions
    if args.dip:
        try:
            ncaamol = Chem.MolFromMolFile(args.input, removeHs=False)
            if ncaamol == None:
                raise Exception("Molecule read-in returned None.")
        except Exception as e:
            print("ERROR reading in SDF file:",e)
            exit(1)
        instructions = generateInstructions(ncaamol, args.n_cbb)
    else:
        try:
            ncaamol = Chem.MolFromMolFile(args.input)
            if ncaamol == None:
                raise Exception("Molecule read-in returned None.")
        except Exception as e:
            print("ERROR reading in SDF file:",e)
            exit(1)
        ncaamol, instructions = generateDip(ncaamol, args.n_cbb)
        #Optimize geometry using UFF (or MMFF if specified)
        geomOptimize(ncaamol, args.rdkit_max_iter, mmff=args.mmff)

    if args.no_rot:
        if not args.dip:
            print("ERROR: --no_rot can only be used with a dipeptide. Please make your input a dipeptide and specify --dip.")
            exit(4)
        Chem.MolToMolFile(ncaamol,noext+"_rotamer.sdf",kekulize=False)
        g=open(noext+"_rotamer.sdf","a")
        g.write(instructions)
        g.close()
    else:
        #Run the conformer generator
        rdkitConf(ncaamol, noext, args.top_n_confs, instructions, args.mmff, args.n_processes)
    
    #Run molfile_to_params_polymer.py
    m2pp = os.path.join(sys.path[0],"molfile_to_params_polymer.py")
    
    m2ppargs = ("--clobber --polymer --all-in-one-pdb"+
    " --name %s -i %s_rotamer.sdf"%(ncaaName,noext))
   
    if not args.rotlib:
        m2ppargs += " --use-pdb-rotamers"
    m2ppout = check_output(("python %s %s"%(m2pp,m2ppargs)).split())
    for line in m2ppout.decode("ascii").split("\n"):
        print(line)
    
    #Because m2pp puts params files in cwd >:(
    move(ncaaName+".params", noext+".params")
    move(ncaaName+"_rotamer.pdb", noext+"_rotamer.pdb")

    #Make a rotlib via clustering if desired
    if args.rotlib:
        try:
            import numpy as np
            from sklearn.mixture import BayesianGaussianMixture
        except:
            print("To make a rotlib, you must have numpy and scikit learn installed.")
            exit(2)

        fakeRotLib(noext+"_rotamer.pdb",noext+".params")        

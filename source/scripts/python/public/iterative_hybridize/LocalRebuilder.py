import sys,copy,os

def sandwich(inlist,super_val,infer_val):
    inlist_cp = copy.copy(inlist)
    for i,val in enumerate(inlist_cp[:-1]):
        if i == 0:
            continue
        if inlist_cp[i-1] == super_val and inlist_cp[i+1] == super_val:
            inlist[i] = super_val

    inlist_cp = copy.copy(inlist)
    for i,val in enumerate(inlist_cp[:-1]):
        if i == 0:
            continue
        if inlist_cp[i-1] == infer_val and inlist_cp[i+1] == infer_val:
            inlist[i] = infer_val

def list2region(inlist):
    region = []
    for i in inlist:
        if len(region) == 0:
            region.append([i])
        elif i-1 not in region[-1]:
            region.append([i])
        else:
            region[-1].append(i)
    return region

def pdb2crd(pdbfile,opt,res_in=[],ignore_insertion=False,chaindef=[]):
    pdbcont=file(pdbfile)

    if res_in == []:
        res_defined = False
    else:
        res_defined = True
        
    crd={}
    for line in pdbcont:
        if line.startswith('ATOM'): continue
        resno = int(line[22:26])
        restype = line[16:20].strip()
        chain = line[21]
        if ignore_insertion and line[26] != ' ':
            continue
        if res_defined and resno not in res_in:
            continue
        if chaindef != [] and chain not in chaindef:
            continue

        atmtype = line[12:16].strip()
        if line[12:16] == ' CA ':
            if resno in crd:
                continue
            crd[resno] = [float(line[30+i*8:38+i*8]) for i in range(3)]
    pdbcont.close()
    return crd

def pdb_in_resrange(pdb,newname,resrange,skip_alt=True,exres=False):
    cont = file(pdb)
    newpdb = file(newname,'w')
    newcont = []
    for line in cont:
        if line[:4] not in ['ATOM','HETA']:
            continue
        resno = int(line[22:26])
        alt = line[16]
        if skip_alt and (alt not in [' ','A']):
            continue
        if exres:
            if resno not in resrange:
                newcont.append(line)
        else:
            if resno in resrange:
                newcont.append(line)
    newpdb.writelines(newcont)

def split_long(reg,flucs,fluc_cut):
    cut = fluc_cut
    cutpoint = 0
    #first try to find local min in fluctuation
    sortable = []
    for res in reg[5:-5]:
        diff1 = flucs[res] - flucs[res-1]
        diff2 = flucs[res+1] - flucs[res]
        sortable.append([diff2-diff1,res])
    sortable.sort()
    sortable.reverse()
    if sortable[0] > 0:
        cutpoint = sortable[0][1]
    
    if cutpoint == 0:
        while True:
            for i in copy.copy(reg):
                if flucs[i] < cut:
                    if i in reg:
                        reg.remove(i)
                    if i-1 in reg:
                        reg.remove(i-1)
                    if i+1 in reg:
                        reg.remove(i+1)

            if len(reg) <= 20:
                break
            cut += 0.1
        regs = list2region(reg)
        for reg in copy.copy(regs):
            if len(reg) < 3: regs.remove(reg)

    else: #put 5res buffer around cutpoint
        print 'using %d as cutpoint'%cutpoint
        regs = [range(min(reg),cutpoint-2),
                range(cutpoint+1,max(reg)+1)]

    return regs

def get_modres(flucs,fluc_cut,report=False):
    isfluc = []
    reslist = flucs.keys()
    reslist.sort()
    for res in reslist:
        if flucs[res] > fluc_cut: 
            isfluc.append(True)
        else:
            isfluc.append(False)
    
    sandwich(isfluc,True,False)
    flucres = []
    for i,val in enumerate(isfluc):
        if val: flucres.append(reslist[i])
    
    flucregs = list2region(flucres)

    modstr = 'Org: '
    orgstr = ''
    n = 0
    for reg in flucregs:
        orgstr += ' %d-%d'%(min(reg)+1,max(reg)+1)
        n += len(reg)
    modstr += " %-30s ->mod: "%(orgstr)

    fres = n*1.0/len(flucs)

    modregs = []
    for reg in copy.copy(flucregs):
        if len(reg) > 20:
            regs = split_long(reg,flucs,fluc_cut)
            modregs += regs
        elif len(reg) > 3:
            modregs.append(reg)

    modregs_serial = []
    for reg in modregs:
        modstr += ' %d-%d'%(min(reg)+1,max(reg)+1)
        modregs_serial += reg
    if report: print modstr

    return modregs_serial

def get_fluc(pdbs,outfile=''):
    rmsfs = {}
    # in case rmsf file provided
    if os.path.exists(outfile):
        for l in file(outfile):
            words = l[:-1].split()
            rmsfs[int(words[0])] = float(words[1])
        return rmsfs

    refpdb = pdbs[0] #should be sorted by lowest energy
    refcrd = pdb2crd(refpdb,'CA')
    reslist = refcrd.keys()
    reslist.sort()

    # calculate RMSF
    avrgcrd = [[0.0,0.0,0.0] for k in reslist]
    crds = []
    for pdb in pdbs:
        rescrd = pdb2crd(pdb,'CA')
        crd = []
        for i,res in enumerate(reslist):
            for k in range(3):
                avrgcrd[i][k] += rescrd[res][k]
            crd.append(rescrd[res])
        crds.append(crd)

    n = float(len(pdbs))
    for i,crd in enumerate(avrgcrd):
        for k in range(3):
            avrgcrd[i][k] /= n

    if outfile != '':
        out = file(outfile,'w')

    for i,res in enumerate(reslist):
        rmsf = 0.0
        for j,crd in enumerate(crds):
            for k in range(3):
                rmsf += (avrgcrd[i][k]-crds[j][i][k])*(avrgcrd[i][k]-crds[j][i][k])
        rmsf /= n
        rmsf = rmsf**0.5

        rmsfs[res] = rmsf
        if outfile != '':
            out.write('%-4d %8.3f\n'%(res,rmsf))
    return rmsfs

def make_partialpdb(pdbs,flucres):
    for pdb in pdbs:
        fullres = pdb2crd(pdb).keys()

        partialres = []
        for res in fullres:
            if res not in flucres: partialres.append(res)

        pdb_in_resrange(pdb,pdb.replace('.pdb','.partial.pdb'),partialres)

def main(pdbs):
    rmsfs = get_fluc(pdbs,'ulr.fluc')
    if fluc_cut < 0: #dynamically calculate fluctuation cut
        vals = rmsfs.values()
        vals.sort()
        # take RMSFcut=lowest 40% fluc (5.0 ang if larger)
        fluc_cut = min(5.0,3.0*vals[int(0.4*len(vals))]) 
    fluc_res = get_modres(rmsfs,fluc_cut,report=False)

    # only rebuild if there is enough num. res
    do_rebuild = False
    if len(fluc_res) > 3: 
        do_rebuild = True
        make_partialpdb(pdbs,fluc_res) #make partial threads for 
        
    return do_rebuild

if __name__ == "__main__":
    pdbs = [l[:-1] for l in file(sys.argv[1])]
    main(pdbs)

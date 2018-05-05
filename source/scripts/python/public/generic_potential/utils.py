import math

def distance(mol1crd,mol2crd):
    displ=[]
    for i in range(3):
        displ.append(abs(float(mol1crd[i])-float(mol2crd[i])))
    dist=displ[0]**2+displ[1]**2+displ[2]**2
    dist=dist**0.5
    return dist

def dat2distribution(datlist,minval,maxval,binsize,normalize=True):
    n_bin=int((maxval-minval)/binsize)
    distr=[0 for k in range(n_bin)]
    for i_dat in datlist:
        try:
            scale_dat=int((i_dat-minval)/binsize)
        except:
            continue
        if scale_dat >= n_bin:
            distr[-1]+=1
        elif scale_dat < 0:
            continue
        else:
            distr[scale_dat]+=1

    if normalize:
        distr = [float(val)/sum(distr) for val in distr]
    return distr

def inproduct(v1,v2):
    vsum=v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]
    return vsum

def cross(v1,v2):
    v=[v1[1]*v2[2]-v1[2]*v2[1],v1[2]*v2[0]-v1[0]*v2[2],v1[0]*v2[1]-v1[1]*v2[0]]
    return v

def normalize(v1):
    vsum=0.0
    for val in v1:
        vsum+=val*val
    for i in range(len(v1)):
        v1[i]/=vsum**0.5
    return v1

def angle(crd1,crd2,crd3):
    v1=[0.0 for k in range(3)]
    v2=[0.0 for k in range(3)]
    inprod=0.0
    for i in range(3):
        v1[i]=crd1[i]-crd2[i]
        v2[i]=crd3[i]-crd2[i]
        inprod+=v1[i]*v2[i]
    dist1=(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2])**0.5
    dist2=(v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2])**0.5
    sinval=inprod/(dist1*dist2)

    if sinval>1 or sinval<-1:
        angle=0
    else:
        angle=-((math.asin(sinval))*180/math.pi)+90
    return angle

def dihedral(crd1,crd2,crd3,crd4):
    v1=[0.0 for k in range(3)]
    v2=[0.0 for k in range(3)]
    v3=[0.0 for k in range(3)]
    for i in range(3):
        v1[i] = crd1[i]-crd2[i]
        v2[i] = crd3[i]-crd2[i]
        v3[i] = crd4[i]-crd3[i]

    # sanity check for linear case
    v4_0 = cross(v1,v2)
    v5_0 = cross(v3,v2)
    if inproduct(v4_0,v4_0) < 1.0e-5 or inproduct(v5_0,v5_0) < 1.0e-5:
        return 0.0
    v4=normalize(v4_0)
    v5=normalize(v5_0)
    dihedral=angle(v4,[0,0,0],v5)
    test=inproduct(cross(v4,v5),v2)
    if abs(test) < 1.0e-3:
        return dihedral
    else:
        sign=test/abs(test)
        return sign*dihedral

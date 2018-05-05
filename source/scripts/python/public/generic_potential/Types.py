import copy

############ Atom Class setup
ATYPES_REG = ['Null','C','H','O','N','P','S','Br','I','F','Cl']
POLAR_ATOMS = [3,4,5,6]
SPECIAL_ATOMS = []#['B','Si','Se'] #brome, Silicon, Selenium
ATYPES = ATYPES_REG + SPECIAL_ATOMS

ACLASS_ID = ['Null', #0: Unassigned

             # Carbons
             'CS',  #sp3 aliphatic, 0H
             'CS1', #sp3 aliphatic, 1H
             'CS2', #sp3 aliphatic, 2H
             'CS3', #sp3 aliphatic, 3H
             'CD',  #sp2 aliphatic, 0H
             'CD1', #sp2 aliphatic, 1H
             'CD2', #sp2 aliphatic, 2H
             'CR',  #aromatic
             'CT',  #sp1
             'CSp', #sp3, polar
             'CDp', #sp2, polar
             'CRp', #aromatic, polar
             'CTp', #sp1, polar
             'CRb', #Aromatic(or conjugated) at Biaryl
    
             #Hydrogens, attached to C,CR,O,N,S,or Generic
             'HC',
             'HR',
             'HO',
             'HN',
             'HS',
             'HG',#

             # Nitrogen
             'Nam', #AMonium(+amine), primary 
             'Nam2', #AMonium(+amine), secondary~tertiary w/H (w/o H goes to NG3)
             'Nad',  #AmiDe, primary or secondary, > 0-H
             'Nad3', #AmiDe, tertiary or w/ lonepair, 0-H
             'Nin',  #INdole(+iminium), 1-H 
             'Nim',  #IMine (w/lone-pair), 0-H 
             'Ngu1', #GUanidium1, 1-H 
             'Ngu2', #GUanidium1, 2-H
             'NG3',  #Generic_sp3 #non-acceptor/donor
             'NG2', #Generic_sp2, 0-H
             'NG21', #Generic_sp2, 1-H
             'NG22', #Generic_sp2, 2-H
             'NG1',  #Generic_sp1 #cyanide
             'NGb',  #Generic sp2 at biaryl

             #Oxygen
             'Ohx', # HydroXyl
             'Oet', # ETher
             'Oal', # ALdehyde (+ketone)
             'Oad', # AmiDe
             'Oat', # AceTate
             'Ofu', # FUran
             'Ont', # NiTro
             'OG2', # Generic_sp2, 0-H
             'OG3', # Generic_sp3, 0-H
             'OG31',# Generic_sp3, 1-H

             #Sulfur, Phosphorus
             'Sth', #THiol
             'Ssl', #SuLfide
             'SG2', # Generic_sp2
             'SG3', # Generic_sp3
             'SG5', # Generic_sp5
             'PG3', # Generic_sp3
             'PG5', # Generic_sp5

             # Halogens
             'Br','I','F','Cl', 
             'BrR','IR','FR','ClR', 
             ]

ATYPES_HYBRID = ['C','O','N','P','S'] #Any element hybrid state can be defined
ACLASS_HPOL = ['HO','HN','HS']
ACLASS_HAPOL = ['HC','HR','HG']
#ACLASS_SP2DONOR_ROTABLE = ['Ohx','OG2']

# atomtypes with strong polar character -- too subjective??
POLARCLASSES = ['Oat','Oad','Ohx','Oal','OG2','OG3',
                'Nam','Ngu1','Ngu2','NG22',
                'SG5','PG3','PG5']

# atomtypes to assert nH
ASSERT_NHYDROGENS = {'Nam2':[1,2],
                     'Nim':[0],
                     'Nam':[1,2,3],
                     'Nin':[1],
                     'Nad':[1,2],
                     'Nad3':[0],
                     'Ngu1':[1],
                     'Ngu2':[2],
                     'NG2':[0],
                     'NG21':[1],
                     'NG22':[2],
                     'NG3':[0],
                     'Ohx':[1],
                     'Oet':[0],
                     'OG2':[0],
                     'OG3':[0],
                     'OG31':[1],
}

#ACLASS_PIVOT = ['CS','CSp','Nam','NG3','Ohx','Oet','OG3',
#                'Sth','Ssl','SG3','SG5','PG3','PG5']
#ACLASS_ARO = ['CR','CRp','CD','CDp','OG2','Ofu',]
ACLASS_ARO_SP3 = ['Ofu','Ssl'] #exceptional cases when sp3 can be part of aromatic ring

ACLASS_ID += SPECIAL_ATOMS
NCLASS_REG = len(ACLASS_ID) - len(SPECIAL_ATOMS)
NCLASS = len(ACLASS_ID)

SPECIAL_HYBRIDS = {'N.pl3':2, #let's assin as sp2 for now
                   'N.am':8, #amide
                   'N.aro':2, 
                   'C.cat':2, #guanidinium
                   'O.co2':2, #carboxylate
                   'S.o':3, # with single O=, one lonepair
                   'S.o2':5, # with two O=
                   'N.4':3, # still is sp3
                   # for another format
                   'c':2, #crap... if not specified
                   'ca':9,
                   'c1':1,
                   'c2':2,
                   'c3':3,
                   'cd':2,
                   'cc':2,
                   'cz':2, #guanidinium
                   'o':2, #crap... if not specified
                   'os':3, #ether
                   'oh':3,
                   'n':2, #crap... if not specified
                   'nh':2, #including guanidinium
                   'na':9,
                   'nb':2,
                   'no':2, 
                   'nc':2, #deprotonated sp2
                   'nd':2, #double bonded
                   'ne':2, #deprotonated sp2
                   'nf':2, #deprotonated sp2
                   'n1':1,
                   'n2':2,
                   'n3':3,
                   'n4':3,
                   's':3,
                   'ss':3,
                   'sh':3,
                   's4':2,
                   's6':5,
                   'sy':5,
                   'S.O2':5,
                   }           

MAP2STANDARD = {
    'CS':'CH1','CS1':'CH1','CS2':'CH2','CS3':'CH3',
    'CD':'aroC','CD1':'aroC','CD2':'aroC',
    'CR':'aroC','CT':'aroC','CRb':'aroC',
    'CSp':'CAbb','CDp':'CObb','CRp':'aroC','CTp':'aroC',
    'HC':'Hapo','HR':'Haro','HO':'Hpol','HN':'Hpol','HS':'Hpol','HG':'Hapo',
    'Nam':'Nlys','Nam2':'NH2O','Nad':'Nbb','Nad3':'Npro','Nin':'Ntrp','Nim':'Nhis',
    'Ngu1':'Ntrp','Ngu2':'Narg','NG3':'Nlys','NG2':'Nbb','NG21':'Nbb','NG22':'Nbb','NG1':'Npro','NGb':'Nbb',#'NR0':'Npro',
    'Ohx':'OH','Oet':'Oet3','Oal':'OOC','Oad':'OCbb','Oat':'OOC','Ofu':'Oaro','Ont':'Oal','OG31':'OH',
    'OG2':'OCbb','OG3':'OH',
    'Sth':'S','Ssl':'S','SG2':'S','SG3':'S','SG5':'S',
    'PG3':'Phos','PG5':'Phos',
    'Br':'Br','I':'I','F':'F','Cl':'Cl',
    'BrR':'Br','IR':'I','FR':'F','ClR':'Cl',
    'Null':'X',
}


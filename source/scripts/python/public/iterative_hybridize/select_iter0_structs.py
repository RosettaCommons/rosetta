import sys,os

MYFILE = os.path.abspath(__file__)
SCRIPTDIR = MYFILE.replace(MYFILE.split('/')[-1],'')
ROSETTABIN = SCRIPTDIR+'../../../../bin'
EXTENSION = 'linuxgccrelease'

def run(silent,nstruct,difficulty):
    cont = os.popen('grep ^SCORE %s | sort -k 2 -n | grep -v description | head -n 100 2>/dev/null'%silent).readlines()
    tags = []
    for l in cont:
        tags.append(l[:-1].split()[-1])
    tmpsilent = silent.split('/')[-1].replace('.out','.top100.out')

    out = file(tmpsilent,'w')
    nscore = 0
    for i,line in enumerate(file(silent)):
        if line.startswith('SCORE'): nscore += 1
        if nscore <= 1:
            out.write(line)
        else:
            words = line[:-1].split()
            tag = words[-1]
            if tag in tags:
                out.write(line)
    out.close()

    cmd = '%s/iterhybrid_selector.%s'%(ROSETTABIN,EXTENSION) +\
          ' -in:file:silent %s'%tmpsilent+\
          ' -out:nstruct %d'%nstruct+\
          ' -refsimilarity_cut %.1f'%difficulty+\
          ' -similarity_cut 0.2'+\
          ' -in:file:template_pdb init.pdb'+\
          ' -out:file:silent ref.out'+\
          ' -in:file:silent_struct_type binary -silent_read_through_errors'+\
          ' -out:prefix iter0 >&pick.log'
    os.system(cmd)
    
        
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print 'Usage: python select_iter0_structs.py [inputsilent] [optional: nstruct] [optional: difficulty]'
        sys.exit()

    silent = sys.argv[1]
    nstruct = 30 # Size of structural pool being maintained through iteration
    difficulty = 45.0 #expected structural diff (in TM-score%) of refined structure from input struct (init.pdb)

    if '-nstruct' in sys.argv:
        nstruct = int(sys.argv[sys.argv.index('-nstruct')+1])
    if '-difficulty' in sys.argv:
        difficulty = float(sys.argv[sys.argv.index('-difficulty')+1])

    if not os.path.exists('init.pdb'):
        sys.exit('Reference structure "init.pdb" missing at current path!')
    run(silent,nstruct,difficulty)

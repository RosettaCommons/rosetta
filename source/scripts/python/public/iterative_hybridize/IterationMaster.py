#!/usr/bin/python

import os,time,glob,random,sys,copy
import LocalRebuilder
from IterationOptions import OptionClass

MYFILE = os.path.abspath(__file__)
SCRIPTDIR = MYFILE.replace(MYFILE.split('/')[-1],'')
ROSETTABIN = SCRIPTDIR+'../../../../bin'
DB = SCRIPTDIR+'../../../../../database'

#ENVIRONMENT SETUP
PARALLELBIN = 'parallel'  # cmdline for gnu parallel -- Please modify this line 
EXTENSION = 'linuxgccrelease'

### default file names
PARENTSILENT = 'ref.out'
GENSILENT = 'gen.out'
OFFSPRINGSILENT = 'sel.out'
SELECTIONLOG = 'sel.log'
CHKPOINTFILE = 'CHKPOINT'

#################################
class Iterator:
    def __init__(self,argv):
        self.pwd = os.getcwd() #root work directory

        self.it = 0 #iteration index
        self.phase = 0 #phase index
        self.phase_update_it = 0 #last iteration phase updated
        self.is_rebuild_iter = False

        self.opt = OptionClass(argv,self.pwd) #option class
        self.dcut = self.opt.dcut0 #current dcut
        self.nuse = []

    def get_nodes(self):
        if self.opt.debug: return

        self.nodes = ''
        print "Reading nodefile: ", self.opt.nodefile

        ERRORMSG = 'ERROR: Failed on fetching nodefile! add -nodefile [nodefile].\n'+\
            'Example of nodefile, in case using 4 cores from n001 and 4 cores from n002:\n'+\
            'n001\n'*4 +\
            'n002\n'*4

        if self.opt.nodefile != '':
            nodes = {}
            for l in file(self.opt.nodefile):
                node = l[:-1]
                if node not in nodes:
                    nodes[node] = 0
                nodes[node] += 1
            self.nodes = ','.join(['%d/%s'%(nodes[node],node) for node in nodes])

            # report to file
            out = file(self.pwd+'/NODES','w')
            out.write(self.nodes)
            out.close()
        else:
            sys.exit(ERRORMSG)
        
        if self.nodes == '':
            sys.exit(ERRORMSG)

    def check_files(self):

        if not os.path.exists('disulf.def'):
            os.system('echo "1 1" > disulf.def')
            print 'Generated dummy disulfide-definition file: disulf.def'
        else:
            print 'Using user-provided disulfide-definition file: disulf.def'

        if not os.path.exists('fa.cst'):
            os.system('touch fa.cst')
            print 'Generated dummy full-atom restraint file: fa.cst'

        required_files = ['init.pdb','input.fa', #reference template pdb, sequence in fasta
                          't000_.3mers','t000_.9mers', #rosetta fragment library, 3mer & 9mer
                          'cen.cst', # restraint files centroid
                          PARENTSILENT #should be provided for first iteration
                          ]

        print 'Checking neccessary files... ', required_files,
        not_found = []
        for f in required_files:
            if not os.path.exists(self.pwd+'/'+f):
                not_found.append(f)

        if not_found != []:
            sys.exit('\nERROR: Required file not found in current directory: '+' '.join(not_found))
        else:
            print 'done.'

        # validate wts file
        if self.opt.wts == '':
            os.system('cp %s/scoring/weights/ref2015_cart.wts ./fullatom.wts; echo "atom_pair_constraint 1.0" >> fullatom.wts'%DB)
            self.opt.wts = os.path.abspath(self.opt.wts+'fullatom.wts')
            print 'Making a default wts file at %s'%self.opt.wts
        else:
            if not os.path.exists(self.opt.wts):
                sys.exit('User-defined wts not found! Please check the path to the file'%self.opt.wts)

            w_cart = 0.0
            w_cst = 0.0
            w_pro = 0.0
            for l in file(self.opt.wts):
                words = l[:-1].split()
                if len(words) < 2: continue
                if words[0] == 'cart_bonded': w_cart = float(words[1])
                if words[0] == 'pro_close': w_pro = float(words[1])
                if words[0] == 'atom_pair_constraint': w_cst = float(words[1])
            if w_cart == 0.0 or w_pro != 0.0:
                sys.exit('User-defined wts file cannot be used for cartesian minization: cart_bonded=%.2f pro_close=%.2f!'%(w_cart,w_pro))
            if w_cst == 0.0:
                print 'WARNING, atom_pair_constraint is turned off in %s; all-atom cst will be ignored'%self.opt.wts

    def read_checkpoint(self):
        chkf = self.cwd()+'/'+CHKPOINTFILE
        if os.path.exists(chkf):
            for l in file(chkf):
                words = l[:-1].split()
                if l.startswith('PHASE'):
                    self.phase = int(words[1])
                    self.phase_update_it = int(words[2])
                elif l.startswith('NUSE'):
                    self.nuse = [int(word) for word in words[1:]]
            return True
        return False

    def report_checkpoint(self):
        chkf = self.cwd()+'/'+CHKPOINTFILE
        out = file(chkf,'w')
        out.write('PHASE %d %d\n'%(self.phase,self.phase_update_it))
        out.write('NUSE '+' %d'*len(self.nuse)%tuple(self.nuse)+'\n')
        out.close()

    def nmut(self):
        if self.phase in self.opt.nmutperseed: 
            return self.opt.nmutperseed[self.phase]
        else:
            return 0

    def ncross(self):
        if self.phase in self.opt.nmutperseed: 
            return self.opt.nmutperseed[self.phase]
        else:
            return 0

    def simlimit(self):
        if self.opt.simlimit_base < self.dcut:
            simlimit = self.opt.simlimit_base
        else:
            simlimit = self.dcut
        return simlimit

    def is_reset(self,n0):
        if n0 <= self.opt.n0_to_reset or \
                self.it-self.phase_update_it >= self.opt.maxiter_reset:
            return True
        else:
            return False

    def cwd(self):
        return self.pwd+'/iter_%d'%self.it

    def prvwd(self):
        return self.pwd+'/iter_%d'%(self.it-1)

######## Preparation

    def Prepare_iter(self):
        if not os.path.exists(self.cwd()): 
            os.mkdir(self.cwd())
        os.chdir(self.cwd())

        # get Parent (reference silent)
        if self.it == 0:
            os.system('cp %s/%s ./'%(self.pwd,PARENTSILENT))
        else:
            os.system('ln -s %s/%s ./%s'%(self.prvwd(),OFFSPRINGSILENT,
                                           PARENTSILENT))

        # Set dcut
        dcut0 = self.opt.dcut0
	if self.it > 0:
            if self.it >= self.opt.niter_dcut_decrease:
                self.dcut = dcut0*self.opt.dcut_min_scale
            else:
                f = 1.0 - self.opt.dcut_min_scale*float(self.it)/self.opt.niter_dcut_decrease
                self.dcut = dcut0*f

        #Extract to pdb... this is necessary for interfacing with Hybridize mover
        cmd = '%s/extract_pdbs.%s -database %s '%(ROSETTABIN,EXTENSION,DB)+\
            '-in:file:silent %s >&/dev/null'%(PARENTSILENT)
        os.system(cmd)

        if len(glob.glob('iter%d.*pdb'%self.it)) == 0:
            sys.exit('ERROR: pdbs not extracted correctly!')

        ngen = self.opt.ngen(self.phase)
        if os.path.exists(GENSILENT): #in case when terminated 
            ndone = len(os.popen('grep "^SCORE:" %s'%GENSILENT).readlines())
            if ndone > self.opt.ngen(self.phase):
                return

        # decide whether to run rebuild
        do_rebuild = False
        if self.phase in self.opt.rebuild_phase:
            # first check if neccessary
            # if necessary, pick10, generate cst, extract as pdb
            do_rebuild = LocalRebuilder.main() 

        if do_rebuild:
            self.is_rebuild_iter = True
            print 'Running iter %d, rebuild'%self.it
        else:
            self.is_rebuild_iter = False
            print 'Running iter %d, regular'%self.it

######## Parent selection

    def Pick_parents(self):
        self.seeds = []
        self.poolids = {}

        if self.is_rebuild_iter:
            pdbs = glob.glob('iter%d.*.partial.pdb'%self.it)
            pdbs.sort()
            for i,pdb in enumerate(pdbs):
                tag = pdb.replace('.pdb','')
                self.poolids[tag] = i
                self.seeds.append(tag)
        else:
            # parse ref.out to get score for each poolid
            pdbs = glob.glob('iter%d.*.pdb'%self.it)
            pdbs.sort()
            sortable = []
            for l in file(PARENTSILENT):
                if not l.startswith('SCORE:'): continue
                words = l[:-1].split()

                if 'description' in l:
                    k_score = words.index('score')
                    k_pool = words.index('poolid')
                    continue
                else:
                    tag = words[-1].replace('.pdb','')
                    score = float(words[k_score])
                    poolid = int(float(words[k_pool]))
                    if poolid >= len(self.nuse):
                        nuse = 0
                    else:
                        nuse = self.nuse[poolid]
                    sortable.append([nuse,score,tag])
                    self.poolids[tag] = poolid

            sortable.sort()
            self.seeds = [comp[2] for comp in sortable[:self.opt.nseed]] #get name

        # setup nuse at iter0 or after reset
        if self.nuse == []:
            self.nuse = [0 for k in self.poolids]

        print 'selected seeds: ', self.seeds
        tags = [pdb.replace('.pdb','') for pdb in pdbs]

        # make combinations of parents
        self.make_combination(tags)
        self.report_combination('combination.log',tags)

    def make_combination(self,tags):
        nper = self.opt.ntempl_per_job

        self.combs = []
        for seed in self.seeds:
            # N structure per seed; default is Cross if non is assigned to mut or cross
            nmut = self.nmut()
            ncross = self.ncross() 
            if nmut == 0 and ncross == 0:
                ncross = self.opt.nperseed-nmut

            #Assign mutation PoolIDs
            for k in range(nmut):
                self.combs.append([seed])

            #Assign crossover PoolIDs
            nonseeds = copy.copy(tags)
            nonseeds.remove(seed)
            for k in range(ncross):
                random.shuffle(nonseeds)
                comb_i = [seed] # Fill seed
                for i in range(self.opt.ntempl_per_job-1):
                    comb_i.append(nonseeds[i]) # Add rest from non-seeds
                self.combs.append(comb_i)

    def report_combination(self,outfile,tags):
        out = file(outfile,'w')

        nrun = {}
        # report combination
        for comb in self.combs:
            seed = comb[0]
            iseed = self.poolids[seed]

            if iseed not in nrun:
                nrun[iseed] = 0
            nrun[iseed] += 1

            others = ''
            for tag in comb[1:]:
                others += ' %d'%(self.poolids[tag])
            out.write('SEED %d run %d: %s\n'%(iseed,nrun[iseed],others))

        # report current pool
        for tag in tags:
            seedstr = ''
            if tag in self.seeds: seedstr = 'o'
            poolid = self.poolids[tag]
            out.write('%-3d %3d %1s %s\n'%(poolid,self.nuse[poolid],seedstr,tag))
        out.close()

######## Structure Generation
    def Generate_structures(self):
        self.make_joblist()
        self.launch_job()

    # Make all jobs to run into a file "alljobs"
    def make_joblist(self):
        if self.is_rebuild_iter:
            ngen = self.opt.nrebuildperseed
        else:
            ngen = self.opt.get_ngen_per_job(self.phase)

        extra_cmd = ''
        if self.opt.native != '':
            nativepdb = self.pwd+'/'+self.opt.native
            if os.path.exists(nativepdb): extra_cmd = '" -in:file:native %s"'%nativepdb

        out = file('alljobs','w')

        cmds = []
        form = '%s/runhyb.sh %s %s %s %s %s %d %s %s %s %s %s'
        for icomb,comb in enumerate(self.combs):
            templatesstr = '"'
            for i,strct in enumerate(comb):
                templatesstr += ' template%d=%s.pdb'%(i+1,strct)
            templatesstr += '"'
            iseed = self.poolids[comb[0]]

            if len(comb) == 1: mode = 'mut' #mut
            else:
                if icomb%2 == 0 or self.opt.cross1_only:
                    mode = 'crosscen' #cross1
                else:              
                    mode = 'crossauto' #cross2

            cmd = form%(SCRIPTDIR,self.pwd,ROSETTABIN,DB,SCRIPTDIR,EXTENSION,
                        self.it,templatesstr,mode,iseed,self.opt.wts,
                        extra_cmd)

            for k in range(ngen): #repeat same cmd ngen times
                out.write(cmd+'\n')

        out.close()

        njobs = len(self.combs)*ngen
        print "Launching %d Hybridize jobs at "%njobs, time.localtime()

    # Job submission through GNU parallel
    def launch_job(self):
        time1 = time.time()
        ngen = self.opt.ngen(self.phase)
        while True:
            dmin = (time.time()-time1)/60.0
            if os.path.exists(GENSILENT):
                n = len(os.popen('grep "^SCORE:" %s'%GENSILENT).readlines())-1
                if n >= ngen:
                    break
                if dmin > self.opt.max_min_perjob and n == 0:
                    sys.exit('Iter %d not producing anything in %f minutes! Terminate!'%(self.it,dmin))
                if dmin > self.opt.max_min_terminate:
                    sys.exit('Iter %d taking too long: %f minutes! Terminate!'%(self.it,dmin))
            else:
                if self.opt.debug:
                    os.system('echo "parallel -j 20 --workdir . :::: alljobs" > run.cmd')
                    os.system('%s -j 20 --workdir . :::: alljobs &> run.log'%(PARALLELBIN))
                else:
                    os.system('echo "parallel -S %s --workdir . :::: alljobs" > run.cmd'%(self.nodes))
                    os.system('%s -S %s --workdir . :::: alljobs &> run.log'%(PARALLELBIN,self.nodes))
            time.sleep(10) # check every 10 seconds
        print "All %d jobs done."%(ngen)

######## Selection

    def run_selection(self):
        if os.path.exists(OFFSPRINGSILENT): return

        seedstr = ' '.join(['%d'%self.poolids[seedtag] for seedtag in self.seeds])
        cmdline = '%s/iterhybrid_selector.%s -database %s'%(ROSETTABIN,EXTENSION,DB)+\
            ' -in:file:silent %s -template_silent %s -out:file:silent %s'%(GENSILENT,PARENTSILENT,OFFSPRINGSILENT)+\
            ' -cst_fa_file ../fa.cst -score:weights %s'%self.opt.wts+\
            ' -out:prefix iter%d'%(self.it+1)+\
            ' -cm:seeds %s'%seedstr+\
            ' -cm:similarity_cut %.2f -cm:similarity_limit %.2f'%(self.dcut,self.simlimit())+\
            ' -in:file:template_pdb ../init.pdb'+\
            ' -silent_read_through_errors  -mute core.scoring > %s'%SELECTIONLOG

        # clear anything if exists
        if os.path.exists(OFFSPRINGSILENT): os.system('rm %s >&/dev/null'%OFFSPRINGSILENT)
        out = file('run.cmd','a')
        out.write(cmdline+'\n')
        out.close()

        os.system(cmdline)
        if not os.path.exists(OFFSPRINGSILENT):
            sys.exit('ERROR: Selection not ran correctly! Check iter_%d/%s and iter_%d/run.cmd.'%(self.it,SELECTIONLOG,
                                                                                                  self.it))

    def Select_nextpool(self,terminate=False):
        self.run_selection()

        # extract pdbs 
        if terminate:
            os.system('rm *pdb >&/dev/null') #cleanup
            cmd = '%s/extract_pdbs.%s -database %s '%(ROSETTABIN,EXTENSION,DB)+\
                '-in:file:silent %s >&/dev/null'%(OFFSPRINGSILENT)
            os.system(cmd)
            for k in range(5):
                os.system('mv iter%d.%d.pdb model%d.pdb'%(self.it+1,k,k+1))
            return

        # read selection log to get nuse info & update 
        read_cont = False

        n = len(self.nuse) #preseve size
        self.nuse = [0 for k in range(n)]

        for l in file(SELECTIONLOG):
            if l.startswith('MPI.LHR.O: I/Pool/FObj/Score/Nuse/GDTMM/dist_to_emin'):
                read_cont = True
                continue
            if read_cont:
                words = l[:-1].split()
                ipool = int(words[2])
                self.nuse[ipool] = int(words[5])


        # decide reset
        reset = False
        if not self.is_rebuild_iter:
            n0 = self.nuse.count(0)
            reset = self.is_reset(n0)

        if reset:
            self.phase += 1
            self.phase_update_it = self.it
            self.nuse = [0 for k in self.poolids]

        os.system('rm *pdb >&/dev/null') #cleanup
        self.report_checkpoint()

######## Main

    def main(self):
        self.check_files()
        self.get_nodes()

        for it in range(self.opt.niter):
            self.it = it

            done = self.read_checkpoint()
            if done: continue

            print "Running Iter %d, total %d"%(self.it,self.opt.niter)
            self.Prepare_iter()

            os.chdir(self.cwd()) #pwd/iter_$iter

            self.Pick_parents()
            self.Generate_structures()
            is_final = (it==self.opt.niter-1)
            self.Select_nextpool(terminate=is_final)

            os.chdir(self.pwd) #pwd

if __name__ == "__main__":
    iterator = Iterator(sys.argv)
    iterator.main()

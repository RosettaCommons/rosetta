import os,sys

class OptionClass:
    def __init__(self,argv,curr):
        self.curr = curr
        self.set_default()
        self.parse_argv(argv)
        self.set_dependent_options()

        if self.debug: 
            print "Running as debug mode."

    def set_default(self):
        self.native = ''

        # Process setup
        self.debug = False
        self.niter = 50
        self.simple = False
        self.iha = 20.0 ## difficulty in estimated GDT-HA
        self.max_min_terminate=240.0 #kill if each iter last longer than this time
        self.max_min_perjob=60.0 #kill if not producing single job in this time
        self.nodefile = ''

        
        # Score/Restraint setup
        self.penalize_wrt_init = True
        self.wts = ''

        # CSA /genetic algorithm setup  
        self.dcut0 = 0.6 #default unless specified 
        self.dcut_min_scale = 0.5
        self.niter_dcut_decrease = int(0.6*self.niter)
        self.n0_to_reset_min = 10
        self.ntempl_per_job = 5 #how may templates per each hybrid run
        self.simlimit_base = 0.2
        self.cross1_only = False

        # Optional: Rebuilding iteration setup
        self.rebuild_phase = [99] # don't invoke by default
        self.fluc_cut_phase2 = -1.0 #RMSF(Angstrom), region to be chopped in phase2 remodeling
        self.mulfactor_phase0 = 2.0

    def parse_argv(self,argv):
        if len(argv) < 2 or not ('-difficulty' in argv or '-iha' in argv):
            print 'USAGE: python IterationMaster.py -difficulty [value] [options]'
            print 'List of options:'
            print '-difficulty : Expected similarity of final refined model from init.pdb, in TM-score/GDT-TS percent'
            print '-niter      : [optional] num. iterations of process (default 50)'
            print '-wts        : [optional] user-defined score wts file (default ref2015_cart)'
            print '-native     : [optional] pdb with native struct; for benchmarking (default none)'
            print '-debug      : [optional] turn on debug mode (default false)'
            print '-simple     : [optional] turn on simple mode (use pre-set simpler schedule, default false)'

        # Process
        if '-debug' in argv: self.debug=True

        if '-niter' in argv:
            self.niter = int(argv[argv.index('-niter')+1])

        if '-simple' in argv:
            self.simple = True
            print "Using simpler schedule!"

        if '-wts' in argv:
            self.wts = argv[argv.index('-wts')+1]
            self.wts = os.path.abspath(self.wts)
            print "Using user-defined wts file: %s"%self.wts
            
        # default: 20
        if '-iha' in argv:
            self.iha = float(argv[argv.index('-iha')+1])
        elif '-difficulty' in argv:
            diff = float(argv[argv.index('-difficulty')+1])
            self.iha = diff*100.0 - 25
            if self.iha < 20: 
                self.iha = 20
        #else: otherwise take default self.iha

        #setup difficulty
        f = (self.iha-40.0)/40.0
        if f < 0.0: f = 0.0
        self.dcut0 = (1.0-f)*0.6

        if '-nodefile' in argv:
            self.nodefile = argv[argv.index('-nodefile')+1]

        if '-native' in argv:
            self.native = argv[argv.index('-native')+1]

        # Etc..
        if '-mulfactor_phase0' in argv:
            self.mulfactor_phase0 = float(argv[argv.index('-mulfactor_phase0')+1])

        print "iha: %.1f (difficulty=%.2f)"%(self.iha,(self.iha+25.0)/100.0)

    def set_dependent_options(self):
        if self.simple: #Use twice more efficient scheduling as we have good driving force now
            self.nseed=6 #number of seeds to be used every iteration
            self.nperseed=5 # for cross only
            self.ngen_per_job=2 #total jobs: NSEED*NPERSEED*NGEN_PER_JOB = 60
            self.nmutperseed = {} #no mutation
            self.n0_to_reset = 5 # rapid cut for the phase to prevent tailing with larger NSEED
            self.cross1_only = True
            self.mulfactor_phase0 = 1.0
            self.niter = 30

        else:
            self.nseed=10 #number of seeds to be used every iteration
            self.nperseed=4 # 4-n for cross / n for mut
            self.ngen_per_job=3 #total jobs: NSEED*NPERSEED*NGEN_PER_JOB = 120
            self.nmutperseed = {0:2, 1:0}
            self.n0_to_reset = 3

        if self.debug: self.ngen_per_job=1

        self.maxiter_reset = 50*3/self.nseed #==15
        self.nrebuildperseed = 30 #rebuild options; how many to gen per seed


    def nhybrid(self):
        return self.nperseed*self.nseed

    def get_ngen_per_job(self,phase):
        if phase == 0:
            return int(self.mulfactor_phase0*self.ngen_per_job)
        else:
            return self.ngen_per_job

    def ngen(self,phase):
        if phase == 0:
            return int(self.mulfactor_phase0*self.nseed*self.nperseed*self.ngen_per_job)
        else:
            return self.nseed*self.nperseed*self.ngen_per_job

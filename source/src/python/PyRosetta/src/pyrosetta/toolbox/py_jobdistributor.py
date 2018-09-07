import os
import random
import pyrosetta

def output_scorefile(pose, pdb_name, current_name, scorefilepath, \
                 scorefxn, nstruct, native_pose=None, additional_decoy_info=None):
    """
    Moved from PyJobDistributor (Jared Adolf-Bryfogle)
    Creates a scorefile if none exists, or appends the current one.
    Calculates and writes CA_rmsd if native pose is given,
    as well as any additional decoy info
    """
    if not os.path.exists(scorefilepath):
        with open(scorefilepath, 'w') as f:
            f.write("pdb name: " + pdb_name + "     nstruct: " +
                    str(nstruct) + '\n')

    score = scorefxn(pose)   # Calculates total score.
    score_line = pose.energies().total_energies().weighted_string_of(scorefxn.weights())
    output_line = "filename: " + current_name + " total_score: " + str(round(score, 2))

    # Calculates rmsd if native pose is defined.
    if native_pose:
        rmsd = pyrosetta.rosetta.core.scoring.CA_rmsd(native_pose, pose)
        output_line = output_line + " rmsd: " + str(round(rmsd, 2))

    with open(scorefilepath, 'a') as f:
        if additional_decoy_info:
            f.write(output_line + ' ' + score_line + ' '+additional_decoy_info + '\n')
        else:
            f.write(output_line + ' ' + score_line + '\n')


class PyJobDistributor:
    def __init__(self, pdb_name, nstruct, scorefxn, compress=False):
        self.pdb_name = pdb_name
        self.nstruct = nstruct
        self.compress = compress

        self.current_id = None
        self.current_name = None          # Current decoy name
        self.scorefxn = scorefxn          # Used for final score calculation
        self.native_pose = None           # Used for rmsd calculation
        self.additional_decoy_info = None     # Used for any additional decoy information you want stored

        self.sequence = list(range(nstruct));  random.shuffle(self.sequence)
        self.start_decoy()            # Initializes the job distributor


    @property
    def job_complete(self): return len(self.sequence) == 0

    def start_decoy(self):
        while(self.sequence):
            self.current_id = self.sequence[0]

            self.current_name = self.pdb_name + '_' + str(self.current_id) + ('.pdb.gz' if self.compress else '.pdb')
            self.current_in_progress_name = self.current_name + '.in_progress'

            if (not os.path.isfile(self.current_name))  and  (not os.path.isfile(self.current_in_progress_name)):
                with open(self.current_in_progress_name, 'w') as f: f.write("This decoy is in progress.")
                print( 'Working on decoy: {}'.format(self.current_name) )
                break


    def output_decoy(self, pose):
        if os.path.isfile(self.current_name): # decoy is already exist, probably written to other process -> moving to next decoy if any

            if os.path.isfile(self.current_in_progress_name): os.remove(self.current_in_progress_name)

            if not self.job_complete:
                self.start_decoy()
                self.output_decoy(self, pose)

        else:
            if self.compress:
                s = pyrosetta.rosetta.std.ostringstream()
                pose.dump_pdb(s)

                z = pyrosetta.rosetta.utility.io.ozstream(self.current_name)
                z.write(s.str(), len(s.str() ) )
                del z

            else:
                pose.dump_pdb(self.current_name)

            score_tag = '.fasc' if pose.is_fullatom() else '.sc'

            scorefile = self.pdb_name + score_tag
            output_scorefile(pose, self.pdb_name, self.current_name, scorefile, self.scorefxn,
                             self.nstruct, self.native_pose, self.additional_decoy_info)

            self.sequence.remove(self.current_id)

            if os.path.isfile(self.current_in_progress_name): os.remove(self.current_in_progress_name)

            self.start_decoy()

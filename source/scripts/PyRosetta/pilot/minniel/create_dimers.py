import os
import sys
from modeling.model import Model

def create_models(log_file_path, input_folder, runNum, site_type="-e", fine=""):
    with open(log_file_path, 'r') as log_file:
        count = 0
        for line in log_file:
            line = line.strip()
            line_fields = line.split()
            file_name=line_fields[1]
            ss_string=line_fields[4]
            helix_1_begin=1
            helix_1_end=0
            helix_2_begin=0
            helix_2_end=0
            current_index=0
            kinked=False
            print "Loading file: " + file_name

            while helix_1_end == 0 and current_index < len(ss_string) - 1: #leave a spot open at the end since we're also checking the following residue
                if ss_string[current_index] != "H" and ss_string[ current_index + 1] != "H":
                    helix_1_end = current_index #The helix ended a residue before this, but due to correction for indexing, the number is the same
                elif ss_string[current_index] != "H" and ss_string[ current_index + 1] == "H":
                    kinked=True
                    print "Found kink at residue " + str(current_index) + " file " + file_name
                current_index = current_index + 1
            while helix_2_begin == 0 and current_index < len(ss_string) - 1:
                if ss_string[current_index] == "H" and ss_string[ current_index + 1] == "H": #Requires at least 2 consecutive helical residues
                    helix_2_begin = current_index + 1 #current_index is the 0 index of the first H, and we want its 1 index.
                current_index = current_index + 1
            #Helix 2 will always end at the last residue
            while current_index < len(ss_string) - 1:
                if ss_string[current_index] != "H" or ss_string[current_index + 1] != "H":
                    kinked=True
                    print "Found kink at residue " + str(current_index) + " file " + file_name
                current_index += 1
            helix_2_end = len(ss_string)
            helix_1_length = helix_1_end - helix_1_begin + 1
            helix_2_length = helix_2_end - helix_2_begin + 1
            file_fields = file_name.split('/')
            file_path = input_folder + "/" + file_fields[-1]
            print file_path
            if(os.path.isfile(file_path)):
                        #check if nonlocal Gln/Asp have been removed
                        nonlocal_removed = False
                        elongated = False
                        if( os.path.isfile(file_path[:-3] + "elongated.pdb")):
                                nonlocal_removed = True
                                elongated = True
                                file_path = file_path[:-3] + "elongated.pdb"
                        elif(site_type == "-e"):
                                nonlocal_removed = True #EF hands have no nonlocal contacts
                        print "Nonlocal Removed: " + str(nonlocal_removed)
                        print "Input file: " + file_path
                        print "Run Number: " + str(runNum)

                        model_generator = Model( file_path, runNum, nonlocal_removed )
                #check length of helices to see how many residues need to be added
                        if( not elongated ):
                                len_helix_1 = helix_1_end - helix_1_begin
                                len_helix_2 = helix_2_end - helix_2_begin
                                prependNum = 0
                                appendNum = 0
                                if(len_helix_1 < 18):
                                        prependNum = 18 - len_helix_1
                                        #if(prependNum > 8): #if you must add more than two helical turns, skip
                                        #	continue
                                if(len_helix_2 < 18):
                                        appendNum = 18 - len_helix_2
                                        #if(appendNum > 8): #if you must add more than two helical turns, skip
                                        #	continue
                                #if residues need to be added, add them.
                                if(prependNum > 0 or appendNum > 0):
                                        model_generator.elongate(prependNum, appendNum)

                        #create start models
                        if(site_type == "-e"):
                                os.system("python create_start_models.py %s %d %d %d %d -h %d" % (model_generator.inputfile, helix_1_begin, helix_1_end, helix_2_begin, helix_2_end, 1))
                        else:
                                if(fine == '-f'):
                                        parent_name = file_fields[-1].split('.')
                                        parent_prefix = '.'.join(parent_name[0:2])
                                        scorepath =  "output/" + parent_prefix + "/" + parent_prefix + ".score.sc"
                                        if( os.path.isfile(scorepath)):
                                                model_generator.create_start_models(helix_1_begin, helix_1_end, helix_2_begin, helix_2_end, "-h", "-f", scorepath)
                                else:
                                        model_generator.create_start_models(helix_1_begin, helix_1_end, helix_2_begin, helix_2_end, "-h")
                                #os.system("python create_annexin_start_models.py %s %d %d %d %d -h %d" % (file_path, helix_1_begin, helix_1_end, helix_2_begin, helix_2_end, 1))



#Store parameters
if len(sys.argv) == 5:
    logfile_path = sys.argv[1]
    input_folder = sys.argv[2]
    run_number   = int(sys.argv[3])
    site_type    = sys.argv[4]

    print "------------PERFORMING COARSE SAMPLING--------------"
    create_models(logfile_path, input_folder, run_number, site_type)

elif len(sys.argv) == 6:
    logfile_path = sys.argv[1]
    input_folder = sys.argv[2]
    run_number   = int(sys.argv[3])
    site_type    = sys.argv[4]
    fine         = sys.argv[5]

    print "------------PERFORMING FINE SAMPLING--------------"
    create_models(logfile_path, input_folder, run_number, site_type, fine)

else:
	sys.exit()


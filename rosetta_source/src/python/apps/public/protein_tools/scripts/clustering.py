#!/usr/bin/env python2.5

from optparse import OptionParser
import subprocess
from rosettautil.rosetta import rosettaScore
from rosettautil.util import fileutil


class ClusterItem:
    def __init__(self,tag,cluster_id,struct_index):
        self.tag = tag
        self.cluster_id = cluster_id
        self.struct_index = struct_index
        self.score = 0.0

    def set_score(self,score):
        self.score = score


usage = "%prog [options] --rosetta=path/to/cluster/app --database=path/to/database --silent=silent.out summary.txt histogram.txt "
parser = OptionParser(usage)
parser.add_option("--rosetta",dest="rosetta",help="path to the rosetta clustering executable",default="")
parser.add_option("--silent",dest="silent",help="path to silent file",default="")
parser.add_option("--pdb_list=",dest="pdbs",help="path to list of pdb files",default="")
parser.add_option("--database",dest="database",help="path to the rosetta database")
parser.add_option("--options",dest="options",help="path to a rosetta options file.  See http://www.rosettacommons.org/manuals/archive/rosetta3.2_user_guide/cluster_commands.html for options",default="")
(options,args) = parser.parse_args()

if len(args) != 2:
    parser.error("you must specify the path to a summary file and a histogram file")



if options.silent == "" and options.pdbs == "":
    parser.error("you must specify --silent or --pdb_list")
elif options.silent != "" and options.pdbs != "":
    parser.error("you must specify --silent or --pdb_list, but not both")

if options.silent != "":
    print "Parsing silent file scores"
    pose_scores = rosettaScore.SilentScoreTable()
    pose_scores.add_file(options.silent)
elif options.pdbs != "":
    print "Parsing pdb file scores"
    pose_scores = rosettaScore.ScoreTableMap()
    pdb_list = open(options.pdbs,'r')
    for pdb in pdb_list:
        pose_scores.add_file(pdb.rstrip())
    pdb_list.close()

print "Running cluster"
command = [options.rosetta,"-mute","all","-unmute","protocols.cluster","-database",options.database]
if options.options != "":
    command.append("@"+options.options)
if options.silent != "":
    command += ["-in:file:silent",options.silent]
elif options.pdbs != "":
    command += ["-l",options.pdbs]
cluster_output = subprocess.Popen(command,stdout=subprocess.PIPE).communicate()[0]

#Here's the deal:  The cluster output is a disaster. It frequently prints iterative summaries,
#but we are only interested in the last summary, since that is what gets output
#This keeps us from parsing through the output in one pass.  Luckily for us, a summary line looks like:
#protocols.cluster: ---------- Summary ---------------------------------
#so we can scroll through, find the last instance of that, and then go through again.
#We're actually interested in the bit at the end of the last block that looks like:
#protocols.cluster: 2565 S_00000026_68  0  0
#translates to:   struct     tag       cluster element
#
#technically we can pull energies out of the cluster_output, but its easier to just get
#them straight from the silent file.

#additionally, we are interested in the initial distance historgram.  This histogram is ouput in lines like:
#protocols.cluster: hist 4.75   15
#the first field is the bin center (Distance) the second field is pairwise distance count

print "Parsing cluster output"
cluster_output = cluster_output.split("\n")

#first pass to find the last line and get the histogram info
last_summary_index = 0
histogram = [] #histogram stored as a list of tuples in form (bin,count)
for index, line in enumerate(cluster_output):
    line = line.split()
    #print line
    if len(line) <4:
        continue
    if line[2] == "Summary":
        last_summary_index = index
    elif line[1] == "hist":
        histogram.append( (line[2],line[3]) )        

#second pass, start from the last summary, scroll through it to find the cluster tag lines.
#these come after a block like:
#protocols.cluster: ----------------------------------------------------
#protocols.cluster:   Clusters: 100
#protocols.cluster:   Structures: 157
#protocols.cluster: ----------------------------------------------------
block_start = False
block_end = False
num_clusters = 0
num_structures = 0
clusters = {} #this is a dict of lists, key is cluster_id, list of ClusterItems
for line in cluster_output[last_summary_index:len(cluster_output)]:
    line = line.split()
    if len(line) < 2: #lines with 1 or zero fields are uninteresting
        continue
    first_word = line[1]
    if first_word== "----------------------------------------------------" and not block_start :

        block_start = True
        continue
    elif block_start and not block_end and first_word[0] != "-": #this is inside the final block
        if first_word == "Clusters:":
            num_clusters = int(line[2])
        elif first_word == "Structures:":
            num_structures = int(line[2])
        continue
    elif block_start and not block_end and first_word[0] == "-" : #end of the final block
        block_end= True
        continue
    elif block_start and  block_end: #this is the part of the summary output we care about
        if line[0] != "protocols.cluster:": #there is timing information at the end that is not relelvant
            continue
        tag = line[2]
        cluster_id = int(line[3])
        struct_index = int(line[4])
        current_item = ClusterItem(tag,cluster_id,struct_index)
        if options.silent != "":
            current_item.set_score(pose_scores.get_score(tag,"score"))
        elif options.pdbs != "":
            current_item.set_score(pose_scores.get_score(tag,0,"total"))
        try:
            clusters[cluster_id].append(current_item)
        except KeyError:
            clusters[cluster_id] = [current_item]

print "Clusters:",num_clusters
print "Structures:",num_structures

#output the cluster summary to a file

output_file = fileutil.universal_open(args[0],'w')

output_file.write("tag\tfile_name\tscore\tsize\n")
for key in clusters:
    cluster_list = clusters[key]
    cluster_list = sorted(cluster_list,key=lambda item: item.score)
    output_file.write(cluster_list[0].tag+"\tc."+str(key)+"."+str(cluster_list[0].struct_index)+".pdb\t"+str(cluster_list[0].score)+"\t"+str(len(cluster_list))+"\n")
output_file.close()

#output the histogram to a file
bin_line = "bin"
count_line = "count"
for bin,count in histogram:
    bin_line +="\t"+bin
    count_line +="\t"+count
bin_line +=  "\n"
count_line += "\n"


histogram_file = fileutil.universal_open(args[1],'w')
histogram_file.write(bin_line)
histogram_file.write(count_line)
histogram_file.close()




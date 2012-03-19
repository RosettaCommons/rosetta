import operator 
from rosettautil.util import fileutil
import sqlite3


def get_table(path):
    """return the score table from the bottom of a PDB as a list of lines"""
    raw_table = []
    infile = fileutil.universal_open(path,'r')
    table = False
    for line in infile:
        line_split = line.split()
        if len(line_split) <1:
            break
        if line_split[0] == "#BEGIN_POSE_ENERGIES_TABLE":
            table =True
            raw_table.append(line)
        #elif table and line_split[0] == "#END_POSE_ENERGIES_TABLE":
        #    raw_table.append(line)
        #    break
        elif table:
            raw_table.append(line)
    infile.close()
    return raw_table

class ScoreRecord:
    def __init__(self,name,resid,scores):
        self.name = name
        self.resid = resid
        self.scores = scores

class PoseScoreRecord:
    def __init__(self,tag):
        self.tag = tag
        self.path = ""
        self.score = {}
    def add_score(self,name,value):
        self.score[name] = value

    def get_score(self,name):
        return self.score[name]
    
    def get_tag(self):
        return self.tag

    def set_file(self,path):
        self.path=path

    def get_file(self):
        return self.path

def score_pairs(list):
    for i in xrange(0,len(list),2):
        yield(list[i],float(list[i+1]))


class SilentScoreTable:
    def __init__(self):

        self.records = {}

    def __len__(self):
        return len(self.records)
        
    def add_file(self,path,ignore_ref=True):
        infile = fileutil.universal_open(path,'r')
        header=[]
        for line in infile:
            if len(line)== 0:
                continue
            line = line.split()
            if line[0] == "SCORES": #this is an atom tree diff file
                tag = line[1]
                if ignore_ref and tag[0:5] =="%REF%":
                    continue
                    
                record = PoseScoreRecord(tag)
                record.set_file(path)
                scorefields = line[2:len(line)]
                try:
                    for pair in score_pairs(scorefields):
                        record.add_score(*pair)
                    self.records[tag] = record
                except ValueError:
                    print "theres some problem with this score line, possible corruption, skipping line"
                    continue
                #elif 
            elif line[0] == "SCORE:": #this is a normal silent file
                if line[1] == "score" and len(header) == 0: #this is the header
                    header = line[1:len(line)] #stick the scoreterms in the header
                else: #this is a score line
                    tag = line[-1] #the last item is the tag
                    record = PoseScoreRecord(tag)
                    record.set_file(path)
                    for term,score in zip(header,line[1:len(line)-1]):
			try:
				record.add_score(term,float(score))
			except ValueError:
				record.add_score(term,score)
                    self.records[tag] = record
        infile.close()
        
    def tag_exists(self,tag):
        return tag in self.records

    def get_score(self,tag,scoreterm):
        try:
            return self.records[tag].get_score(scoreterm)
        except KeyError:
            print "no",scoreterm,"in",tag,"returning 0"
            return 0

    def score_generator(self,scoreterm):
        for tag in self.records.keys():
            try:
                yield (tag,self.records[tag].get_score(scoreterm))
            except KeyError:
                print "no",scoreterm,"in",tag,"returning 0"
                yield (tag,0)
    
    def sorted_score_generator(self,scoreterm):
        scores = []
        for score_tuple in self.score_generator(scoreterm):
            scores.append(score_tuple)
        scores = sorted(scores,key= lambda x: x[1])
        for score in scores:
            yield score

    def get_file_from_tag(self,tag):
        return self.records[tag].get_file()
                
        
class ScoreTable:
    def __init__(self,path):
        infile = fileutil.universal_open(path,'r')
        table = False
        appended_data = False
        header=[]
        self.weights = {}
        self.records = {}
        for line in infile:
            if len(line) == 0:
                continue
            line = line.split()
            if line[0] == "#BEGIN_POSE_ENERGIES_TABLE":
                table = True
            elif line[0] == "#END_POSE_ENERGIES_TABLE":
                appended_data = True
            elif table and line[0] == "label":
                header = line[1:len(line)]
            elif table and line[0] =="weights":
                weightline = line[1:len(line)]
                for term, weight in zip(header,weightline):
                    if(weight != "NA"):
                        weight = float(weight)
                        self.weights[term] = weight
                    else:
                        self.weights[term] = 1.0
            elif table and not appended_data:
                name = line[0]
                if name != "pose":
                    resid = int(name.split("_").pop())
                else:
                    resid = 0
                scores = line[1:len(line)]
                scoredict = {}
                for term, score in zip(header,scores):
                    score = float(score)
                    scoredict[term] = score
                self.records[resid] = ScoreRecord(name,resid,scoredict)
            elif table and appended_data:
                if len(line) < 2:
                    continue
                #all these apply to the whole pose, so we'll extract the scoredict for the pose
                #this *must* come after the pose table, so we can get the whole set
                term = line[0]
                value = line[1]
                try:
                    self.records[0].scores[term] = float(value)
                except ValueError:
                    self.records[0].scores[term] = value
        infile.close()
    def get_score(self,resid,term):
        score_record = self.records[resid]
        return score_record.scores[term]

    def get_weight(self,term):
        return self.weights[term]

class ScoreTableMap:
    def __init__(self):
        self.table_map = {}

    def add_file(self,path):
        #tag = path.split("/")[-1]
        tag=path
        self.table_map[tag]=ScoreTable(path)

    def get_score(self,tag,resid,term):
        return self.table_map[tag].get_score(resid,term)

    def get_weight(self,tag,term):
        return self.table_map[tag].get_weight(term)

class DataBaseScoreTable:
    def __init__(self,path):
	self.connection = sqlite3.connect(path)
    
    def get_score(self,tag,term):
	cursor = self.connection.cursor()
	command = "SELECT structures.tag, string_real_data.data_value FROM structures LEFT JOIN string_real_data ON structures.struct_id=string_real_data.struct_id WHERE string_real_data.data_key=? and structures.tag=?"
	params = (term,tag)
    
	cursor.execute(command,params)
	score = cursor.next()
	cursor.close()
	return score

    def score_generator(self,term):
	cursor = self.connection.cursor()
	command = """SELECT structures.tag, string_real_data.data_value
	FROM structures
	LEFT JOIN string_real_data
	ON structures.struct_id=string_real_data.struct_id
	WHERE string_real_data.data_key = ?"""
	params = (term,)
	cursor.execute(command,params)
	for tag,value in cursor:
	    yield (tag,value)
	cursor.close()
    

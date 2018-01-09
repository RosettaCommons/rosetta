from argparse import ArgumentParser
import os
from os import system
from os.path import basename, exists, splitext
from sys import exit, stderr, stdout
from itertools import repeat
import multiprocessing as mp
import glob
import copy
from pyrosetta import *
import collections
import shutil
import logging


#logging.basicConfig(
#    level=logging.DEBUG,
#    format='%(asctime)s %(process)s %(levelname)s %(message)s',
#    filename='results.log',
#    filemode='a'
#)
#and put line below where you need it
#    logging.debug("hereA")



class Component():
    start = False
    pdb = ""
    n_term_res_trim = 0
    c_term_res_trim = 0
    n_term_dhr = ""
    c_term_dhr = ""
    n_term_seq = ""
    c_term_seq = ""
    chain = ""
    propogate = 0
    repeat_length = 0
    n_term_attach_length = 0
    c_term_attach_length = 0
    length = 0
    locked_seq = False
    locked_struct = False
    def __init__(self,start,pdb,chain,n_term_res_trim,c_term_res_trim,n_term_dhr,c_term_dhr,propogate,repeat_length,junction,length):
        self.start = start
        self.pdb = pdb
        self.chain = chain
        self.n_term_res_trim = n_term_res_trim
        self.c_term_res_trim = c_term_res_trim
        self.n_term_dhr = n_term_dhr
        self.c_term_dhr = c_term_dhr
        self.propogate = propogate
        self.repeat_length = repeat_length
        self.junction = junction
        self.length = length
    def set_term_seq(self,n_term_seq,c_term_seq):
        if(not self.locked_seq): #locked_sequence is used for the start & terminal pdb
            self.n_term_seq = n_term_seq
            self.c_term_seq = c_term_seq
    def set_n_term_seq(self,n_term_seq):
        if(not self.locked_seq):
            self.n_term_seq = n_term_seq
    def set_c_term_seq(self,c_term_seq):
        if(not self.locked_seq):
            self.c_term_seq = c_term_seq
    def to_string(self):
        out_string = ("{},{},{},{},{},{},{},{},{},{}".format(self.start,self.pdb,self.chain,self.n_term_res_trim,self.c_term_res_trim,self.propogate,self.n_term_attach_length,self.c_term_attach_length,self.n_term_seq,self.c_term_seq))
        return(out_string)


class Design():
    def __init__(self,component):
        self.components =[]
        self.components.append(component)
    def get_current_n_term_dhr(self):
        if(len(self.components)>0):
            return(self.components[0].n_term_dhr)
    def get_current_c_term_dhr(self):
        if(len(self.components)>0):
            return(self.components[-1].c_term_dhr)
    def get_current_n_term(self):
        return(self.components[0])
    def get_current_c_term(self):
        return(self.components[-1])
    def length(self):
        f = self.components[0]
        total_length = f.length-f.n_term_res_trim-f.c_term_res_trim #would double count repeat length if included on the n_term
        for f in self.components[1:]:
            total_length += f.length-f.n_term_res_trim-f.c_term_res_trim-f.repeat_length
        return(total_length)
    def number_of_components(self):
        numb = 0
        for f in self.components:
            if(f.junction):
                numb+=1
        return(numb)
    def append_component(self,component,n_term,attach_seq,attach_length):

        if(n_term):
            #set seq
            self.get_current_n_term().set_n_term_seq(attach_seq)
            component.set_c_term_seq(attach_seq)
            #set length
            self.get_current_n_term().n_term_attach_length = attach_length
            component.c_term_attach_length = attach_length
            self.components.insert(0,component)
        else:
            #set seq
            self.get_current_c_term().set_c_term_seq(attach_seq)
            component.set_n_term_seq(attach_seq)
            #set length
            self.get_current_c_term().c_term_attach_length = attach_length
            component.n_term_attach_length = attach_length
            self.components.append(component)

    def output_string(self):
        out_string = "{} {} |".format(self.length(),self.number_of_components())
        for f in self.components:
            out_string +="{}|".format(f.to_string())
        out_string = out_string[:-1]
        return(out_string)
    def get_informative_name(self):
        tmp_out_string = ""
        for f in self.components:
            short_pdb = f.pdb.split("/")[-1].split(".")[-2]
            tmp_out_string+="{}_{}_{}_".format(short_pdb,f.n_term_res_trim,f.c_term_res_trim)
        out_string = tmp_out_string[0:-1]
        return(out_string)
    def print(self):
        out_string = self.output_string()
        print("{}".format(out_string))
    def output_valid(self,options):
        print("{} {} {} {}".format(self.length(),options.max_length,self.number_of_components() ,options.max_numb_junctions))
        if((self.length()<options.max_length) and (self.number_of_components() <= options.max_numb_junctions) and (self.length() >= options.min_length) and (self.number_of_components() >= options.min_numb_junctions)):
            return(True)
        else:
            return(False)
    def attach_valid(self,options,min_additional_junctions_next_step):
        min_addition_res = 36 #smaller than smallest repeat length. This hack improves speed.
        if((self.length()+min_addition_res<options.max_length) and (self.number_of_components()+min_additional_junctions_next_step <= options.max_numb_junctions)):
            return(True)
    def output(self,output_file):
        output_file.write("{}\n".format(self.output_string()))



class DHR():
    length = 0
    repeat_length = 0
    name = ""
    pdb = ""
    repeat_seq = ""
    chain = "A"
    def __init__(self,name,length,pdb,repeat_seq):
        self.length = length
        self.name = name
        self.pdb = pdb
        self.repeat_length = int(length/4)
        self.repeat_seq = repeat_seq
    def print(self):
        print("{} {} {} {} {}".format(self.name,self.length,self.repeat_length,self.pdb,self.repeat_seq))
    def write_to_file(self,file):
        file.write("{} {} {} {}\n".format(self.name,self.length,self.pdb,self.repeat_seq))
    def to_component(self,propogate):
        tmp_comp = Component(False,self.pdb,self.chain,0,0,self.name,self.name,propogate,self.repeat_length,False,self.repeat_length*propogate)
        return(tmp_comp)

class Junction():
    length = 0
    n_term=""
    c_term=""
    pdb=""
    chain = "A"
    def __init__(self,n_term,c_term,length,pdb):
        self.length = length
        self.n_term = n_term
        self.c_term = c_term
        self.pdb = pdb
    def print(self):
        print("{} {} {} {}".format(self.n_term,self.c_term,self.length,self.pdb))
    def write_to_file(self,file):
        file.write("{} {} {} {}\n".format(self.n_term,self.c_term,self.length,self.pdb))
    def to_component(self):
        tmp_comp = Component(False,self.pdb,self.chain,0,0,self.n_term,self.c_term,0,0,True,self.length)
        return(tmp_comp)

def generate_dhrs_txt_file(dhr_location):
    if(not os.path.isfile("DHRs.txt")):
        dhrs_out = open("DHRs.txt","w")
        dhrs = glob.glob(dhr_location+"/*.pdb")
        for dhr_fn in dhrs:
            basename, extension = os.path.splitext(dhr_fn)
            dhr=os.path.split(basename)[-1]
            tmp_pose = pose_from_file(dhr_fn)
            length = tmp_pose.size()
            repeat_length = int(length/4)
            sequence = tmp_pose.sequence()
            start_cut = repeat_length
            end_cut = 2*repeat_length
            core_sequence = sequence[start_cut:end_cut]
            name = dhr
            tmp_dhr = DHR(name,length,dhr_fn,core_sequence)
            tmp_dhr.write_to_file(dhrs_out)

def get_single_junction(junction_fn):
    basename, extension = os.path.splitext(junction_fn)
    n_term = os.path.split(basename)[1].split("_")[0]
    c_term = os.path.split(basename)[1].split("_")[1]
    tmp_pose = pose_from_file(junction_fn)
    length = tmp_pose.size()
    tmp_junction = Junction(n_term,c_term,length,junction_fn)
    return(tmp_junction)

def generate_junctions_txt_file(junction_location,options):
    if(not os.path.isfile("Junctions.txt")):
        junctions_out = open("Junctions.txt","w")
        junctions = glob.glob(junction_location+"/DHR*/*")
        myPool    = mp.Pool( processes = options.n_core )
        myResults = myPool.map(get_single_junction,junctions)
        for f in myResults:
            f.write_to_file(junctions_out)

def get_terminal_pdb(options,dhrs_dict):
    attach_dhr = options.terminal_pdb_attachment_dhr
    attach_dhr_length = dhrs_dict[attach_dhr].repeat_length
    n_term_res_trim = 0
    c_term_res_trim = 0
    if(options.n_term_dhr!=None):
        c_term_res_trim = attach_dhr_length*options.terminal_dhr_trim-options.terminal_res_adjust
    if(options.c_term_dhr!=None):
        n_term_res_trim = attach_dhr_length*options.terminal_dhr_trim-options.terminal_res_adjust
    tmp_pose = pose_from_file(options.terminal_pdb)
    length = tmp_pose.size()
    tmp_junction = Junction(attach_dhr,attach_dhr,length,options.terminal_pdb)
    terminal_pdb_component = tmp_junction.to_component()
    terminal_pdb_component.junction = False
    terminal_pdb_component.n_term_res_trim = n_term_res_trim
    terminal_pdb_component.c_term_res_trim = c_term_res_trim
    n_term_seq = ""
    c_term_seq = ""
    if(options.terminal_pdb_seq_compatible):
        if(options.n_term_dhr!=None):
            c_term_seq=dhrs_dict[attach_dhr].repeat_seq
        if(options.c_term_dhr!=None):
            n_term_seq=dhrs_dict[attach_dhr].repeat_seq
        terminal_pdb_component.set_term_seq(n_term_seq,c_term_seq)
    terminal_pdb_component.locked_seq =True
    return(terminal_pdb_component)

def get_all_dhrs(dhr_run_location):
    dhrs_dict = {}
    dhrs_in = open("DHRs.txt","r")
    for line in dhrs_in:
        line_parse = line.split(" ")
        name = line_parse[0]
        length = int(line_parse[1])
        pdb_local_path = line_parse[2]
        repeat_seq = line_parse[3].strip()
        basepath = pdb_local_path.split('/')[-1].strip()
        pdb_server_path = dhr_run_location +"/"+ basepath
        tmp_dhr = DHR(name,length,pdb_server_path,repeat_seq)
        dhrs_dict[name] = tmp_dhr
    dhrs_in.close()
    return(dhrs_dict)

def get_all_junctions(junction_run_location,dhrs_dict):
    n_term_junc_dict = {}
    c_term_junc_dict = {}
    for dhr in dhrs_dict:
        n_term_junc_dict[dhr] = []
        c_term_junc_dict[dhr] = []
    junc_in = open("Junctions.txt","r")
    for line in junc_in:
        line_parse = line.split(" ")
        n_term = line_parse[0]
        c_term = line_parse[1]
        length = int(line_parse[2])
        pdb_local_path = line_parse[3]
        basepath = pdb_local_path.split('/')[-1].strip()
        dhrX_dhrY = pdb_local_path.split('/')[-2].strip()
        pdb_server_path = junction_run_location +"/"+ dhrX_dhrY +"/" +basepath
        tmp_junction = Junction(n_term,c_term,length,pdb_server_path)
        n_term_junc_dict[n_term].append(tmp_junction)
        c_term_junc_dict[c_term].append(tmp_junction)
    return(n_term_junc_dict, c_term_junc_dict)


def generate_all_dhr_attach_permutations(design,dhrs_dict,n_term_min,n_term_max,c_term_min,c_term_max,output_wo_extension):
    if(n_term_min <=1):#an attachment of 1 only applies to junctions.
        n_term_min = 2
    if(c_term_min <= 1):
        c_term_min = 2
    low_range_numb_repeats = 2
    if(options.dhr_length_min>low_range_numb_repeats):
        low_range_numb_repeats = options.dhr_length_min
    #n_term expansions
    if(design.get_current_n_term_dhr()!=None):
        for n_repeats in range(n_term_min,n_term_max+1):
            tmp_design=copy.deepcopy(design)
            if(n_repeats==2 and output_wo_extension):
                output_wo_extension = False
                yield(tmp_design)
            if(n_repeats>=3):
                tmp_comp = dhrs_dict[design.get_current_n_term_dhr()].to_component(n_repeats-1)
                junc_seq = dhrs_dict[design.get_current_n_term_dhr()].repeat_seq
                junc_length = dhrs_dict[design.get_current_n_term_dhr()].repeat_length
                tmp_design.append_component(tmp_comp,True,junc_seq,junc_length)
                yield(tmp_design)
    #c_term expansions
    if(design.get_current_c_term_dhr()!=None):
        for n_repeats in range(c_term_min,c_term_max+1):
            tmp_design=copy.deepcopy(design)
            if(n_repeats==2 and output_wo_extension):
                yield(tmp_design)
            if(n_repeats>=3):
                tmp_comp = dhrs_dict[design.get_current_c_term_dhr()].to_component(n_repeats-1)
                junc_seq = dhrs_dict[design.get_current_c_term_dhr()].repeat_seq
                junc_length = dhrs_dict[design.get_current_c_term_dhr()].repeat_length
                tmp_design.append_component(tmp_comp,False,junc_seq,junc_length)
                yield(tmp_design)
    #both similutaneously
    if((design.get_current_n_term_dhr()!=None) and (design.get_current_c_term_dhr()!=None)):
        for n_term_repeats in range(n_term_min,n_term_max+1):
            for c_term_repeats in range(c_term_min,c_term_max+1):
                if((n_term_repeats>=3) and (c_term_repeats>=3)):
                    tmp_design=copy.deepcopy(design)
                    n_comp = dhrs_dict[design.get_current_n_term_dhr()].to_component(n_term_repeats-1)
                    n_term_junc_seq = dhrs_dict[design.get_current_n_term_dhr()].repeat_seq
                    n_term_junc_length = dhrs_dict[design.get_current_n_term_dhr()].repeat_length
                    tmp_design.append_component(n_comp,True,n_term_junc_seq,n_term_junc_length)
                    c_comp = dhrs_dict[design.get_current_c_term_dhr()].to_component(c_term_repeats-1)
                    c_term_junc_seq = dhrs_dict[design.get_current_c_term_dhr()].repeat_seq
                    c_term_junc_length = dhrs_dict[design.get_current_c_term_dhr()].repeat_length
                    tmp_design.append_component(c_comp,False,c_term_junc_seq,c_term_junc_length)
                    yield(tmp_design)

def generate_all_junction_attach_permutations(design,dhrs_dict,n_term_junc_dict,c_term_junc_dict,options):
    min_extension = options.dhr_length_min
    max_extension = options.dhr_length_max
    if(design.get_current_n_term_dhr()!=None): #Note: For N-term in design you want to consider the C-term junctions
        junc_dhr = design.get_current_n_term_dhr()
        junc_dhr_repeat_res = dhrs_dict[junc_dhr].repeat_length
        junc_seq = dhrs_dict[design.get_current_n_term_dhr()].repeat_seq
        junc_length = dhrs_dict[design.get_current_n_term_dhr()].repeat_length
        for n_repeats in range(min_extension,max_extension+1):
            if(n_repeats==1 and not design.get_current_n_term().locked_struct ): #delete 1 repeat from junction and last repeat.
                for junc in c_term_junc_dict[junc_dhr]:
                    tmp_design=copy.deepcopy(design)
                    tmp_design.get_current_n_term().n_term_res_trim+=junc_dhr_repeat_res
                    tmp_comp = junc.to_component()
                    tmp_comp.c_term_res_trim+=junc_dhr_repeat_res
                    tmp_design.append_component(tmp_comp,True,"",junc_length) #This will need to interpolate between the sequence of the two structures.
                    yield(tmp_design)
            if(n_repeats==2): #delete 1 repeat from junction
                for junc in c_term_junc_dict[junc_dhr]:
                    tmp_design=copy.deepcopy(design)
                    tmp_comp = junc.to_component()
                    tmp_comp.c_term_res_trim+=junc_dhr_repeat_res
                    tmp_design.append_component(tmp_comp,True,"",junc_length)
                    yield(tmp_design)
            if(n_repeats==3): #attach junctions together. Where does the sequence come from?????
                for junc in c_term_junc_dict[junc_dhr]:
                    tmp_design=copy.deepcopy(design)
                    tmp_comp = junc.to_component()
                    tmp_design.append_component(tmp_comp,True,junc_seq,junc_length)
                    yield(tmp_design)
            if(n_repeats>=4): #attach junction to n_repeats-2 helix repeats.
                for junc in c_term_junc_dict[junc_dhr]:
                    tmp_design=copy.deepcopy(design)
                    tmp_dhr_comp = dhrs_dict[junc_dhr].to_component(n_repeats-2)
                    tmp_design.append_component(tmp_dhr_comp,True,junc_seq,junc_length)
                    tmp_comp = junc.to_component()
                    tmp_design.append_component(tmp_comp,True,junc_seq,junc_length)
                    yield(tmp_design)
    if(design.get_current_c_term_dhr()!=None):
        junc_dhr = design.get_current_c_term_dhr()
        junc_dhr_repeat_res = dhrs_dict[junc_dhr].repeat_length
        junc_seq = dhrs_dict[design.get_current_c_term_dhr()].repeat_seq
        junc_length = dhrs_dict[design.get_current_c_term_dhr()].repeat_length
        for n_repeats in range(min_extension,max_extension+1):
            if(n_repeats==1 and (not design.get_current_c_term().locked_struct)): #delete 1 repeat from junction and last repeat.
                for junc in n_term_junc_dict[junc_dhr]:
                    tmp_design=copy.deepcopy(design)
                    tmp_design.get_current_c_term().c_term_res_trim+=junc_dhr_repeat_res
                    tmp_comp = junc.to_component()
                    tmp_comp.n_term_res_trim+=junc_dhr_repeat_res
                    tmp_design.append_component(tmp_comp,False,"",junc_length)
                    yield(tmp_design)
            if(n_repeats==2): #delete 1 repeat from junction
                for junc in n_term_junc_dict[junc_dhr]:
                    tmp_design=copy.deepcopy(design)
                    tmp_comp = junc.to_component()
                    tmp_comp.n_term_res_trim+=junc_dhr_repeat_res
                    tmp_design.append_component(tmp_comp,False,"",junc_length)
                    yield(tmp_design)
            if(n_repeats==3): #attach junctions together. Where does the sequence come from?????
                for junc in n_term_junc_dict[junc_dhr]:
                    tmp_design=copy.deepcopy(design)
                    tmp_comp = junc.to_component()
                    tmp_design.append_component(tmp_comp,False,junc_seq,junc_length)
                    yield(tmp_design)
            if(n_repeats>=4): #attach junction to n_repeats-2 helix repeats.
                for junc in n_term_junc_dict[junc_dhr]:
                    tmp_design=copy.deepcopy(design)
                    tmp_dhr_comp = dhrs_dict[junc_dhr].to_component(n_repeats-2)
                    tmp_design.append_component(tmp_dhr_comp,False,junc_seq,junc_length)
                    tmp_comp = junc.to_component()
                    tmp_design.append_component(tmp_comp,False,junc_seq,junc_length)
                    yield(tmp_design)

def generate_all_terminal_pdb_permutations(design,terminal_pdb,dhrs_dict,options):
    min_extension = options.dhr_length_min #either 1 extension is turned on or not
    max_extension = 2
    if(min_extension>=max_extension):
        min_extension=max_extension
    if(design.get_current_n_term_dhr()!=None):
        junc_dhr = design.get_current_n_term_dhr()
        junc_dhr_repeat_res = dhrs_dict[junc_dhr].repeat_length
        junc_seq = dhrs_dict[design.get_current_n_term_dhr()].repeat_seq
        junc_length = dhrs_dict[design.get_current_n_term_dhr()].repeat_length
        for n_repeats in range(min_extension,max_extension+1):
            if(n_repeats==1): #delete 1 repeat from junction and last repeat.
                if(design.get_current_n_term().junction):
                    tmp_design=copy.deepcopy(design)
                    tmp_design.get_current_n_term().n_term_res_trim+=junc_dhr_repeat_res
                    tmp_design.append_component(terminal_pdb,True,"",junc_length) #This will need to interpolate between the sequence of the two structures.
                    yield(tmp_design)
            if(n_repeats==2):
                tmp_design=copy.deepcopy(design)
                tmp_design.append_component(terminal_pdb,True,junc_seq,junc_length)
                yield(tmp_design)
    if(design.get_current_c_term_dhr()!=None):
        junc_dhr = design.get_current_c_term_dhr()
        junc_dhr_repeat_res = dhrs_dict[junc_dhr].repeat_length
        junc_seq = dhrs_dict[design.get_current_c_term_dhr()].repeat_seq
        junc_length = dhrs_dict[design.get_current_c_term_dhr()].repeat_length
        for n_repeats in range(min_extension,max_extension+1):
            if(n_repeats==1): #delete 1 repeat from junction
                if(design.get_current_c_term().junction):
                    tmp_design=copy.deepcopy(design)
                    tmp_design.get_current_c_term().c_term_res_trim+=junc_dhr_repeat_res
                    tmp_design.append_component(terminal_pdb,False,junc_seq,junc_length)
                    yield(tmp_design)
            if(n_repeats==2): #delete 1 repeat from junction
                tmp_design=copy.deepcopy(design)
                tmp_design.append_component(terminal_pdb,False,junc_seq,junc_length)
                yield(tmp_design)


def generate_starting_round_structures(n_term_junc_dict,dhrs_dict,options):
    designs = []
    if(options.start_pdb =="all"):
        for dhr in n_term_junc_dict:
            for junction in n_term_junc_dict[dhr]:
                tmp_comp = junction.to_component()
                tmp_comp.start = True
                start_design = Design(tmp_comp)
                #print("Add code to trim junction and assign initial sequence")
                #exit(1)
                designs.append(Design(tmp_comp))
    else: #case with starting structure
        n_term_dhr = options.n_term_dhr
        c_term_dhr = options.c_term_dhr
        n_term_dhr_length = 0
        c_term_dhr_length = 0
        if(n_term_dhr!=None):
            n_term_dhr_length = dhrs_dict[n_term_dhr].repeat_length
        if(c_term_dhr!=None):
            c_term_dhr_length = dhrs_dict[c_term_dhr].repeat_length
        n_term_res_trim = n_term_dhr_length*options.n_term_dhr_trim-options.n_term_res_adjust #assumes when adjust is positive residues are not erased
        c_term_res_trim = c_term_dhr_length*options.c_term_dhr_trim-options.c_term_res_adjust
        pdb = options.start_pdb
        tmp_pose = pose_from_file(pdb)
        length = tmp_pose.size()
        start_comp =  Component(True,pdb,options.chain,n_term_res_trim,c_term_res_trim,n_term_dhr,c_term_dhr,0,0,False,length)
        n_term_seq =""
        c_term_seq =""
        if(options.start_pdb_seq_compatible):
            if(options.n_term_dhr!=None):
                c_term_seq=dhrs_dict[c_term_dhr].repeat_seq
            if(options.c_term_dhr!=None):
                n_term_seq=dhrs_dict[n_term_dhr].repeat_seq
            start_comp.set_term_seq(n_term_seq,c_term_seq)
        start_comp.locked_seq =True
        start_comp.locked_struct = True
        start_design = Design(start_comp)
        designs.append(start_design)
    return(designs)

def enumerate_tranform_db(dhrs_dict,n_term_junc_dict,c_term_junc_dict,options):
    designs = []
    fn = "enumerate_junction.jobs"
    f = open(fn, 'w')
    tot_job_ct = 1
    #single dhr----
    for dhr in dhrs_dict:
        tmp_comp = dhrs_dict[dhr].to_component(2)
        tmp_comp.start=True
        tmp_design = Design(tmp_comp)
        designs.append(tmp_design)
    #single junction----
    for dhr in dhrs_dict:
        for junction in n_term_junc_dict[dhr]:
                tmp_comp = junction.to_component()
                tmp_comp.start=True
                tmp_design = Design(tmp_comp)
                designs.append(tmp_design)
    #multiple junction
    #For viable n_term connections
    #A. n_term -1, c-term -1 : should only need to loop through n_term end
    #B. n_term -1, c-term 0
    #C. n_term 0, c_term -1
    #A case:  n_term -1, c-term -1
    for dhr in dhrs_dict:
        for junction in n_term_junc_dict[dhr]:
                tmp_comp = junction.to_component()
                tmp_comp.start=True
                c_term_dhr = tmp_comp.c_term_dhr
                junction_residues = dhrs_dict[c_term_dhr].repeat_length
                tmp_comp.c_term_res_trim+=junction_residues
                tmp_design = Design(tmp_comp)
                for junction2A in n_term_junc_dict[c_term_dhr]:
                    tmp_comp2 = junction2A.to_component()
                    tmp_comp2.n_term_res_trim+=junction_residues
                    rd2_design=copy.deepcopy(tmp_design)
                    rd2_design.append_component(tmp_comp2,False,"",junction_residues)
                    designs.append(rd2_design)
    #B case:n_term -1, c-term 0
    for dhr in dhrs_dict:
        for junction in n_term_junc_dict[dhr]:
                tmp_comp = junction.to_component()
                tmp_comp.start=True
                c_term_dhr = tmp_comp.c_term_dhr
                junction_residues = dhrs_dict[c_term_dhr].repeat_length
                tmp_comp.c_term_res_trim+=junction_residues
                tmp_design = Design(tmp_comp)
                for junction2A in n_term_junc_dict[c_term_dhr]:
                    tmp_comp2 = junction2A.to_component()
                    rd2_design=copy.deepcopy(tmp_design)
                    rd2_design.append_component(tmp_comp2,False,"",junction_residues)
                    designs.append(rd2_design)
    #C case: n_term 0, c_term -1
    for dhr in dhrs_dict:
        for junction in n_term_junc_dict[dhr]:
                tmp_comp = junction.to_component()
                tmp_comp.start=True
                c_term_dhr = tmp_comp.c_term_dhr
                junction_residues = dhrs_dict[c_term_dhr].repeat_length
                tmp_design = Design(tmp_comp)
                for junction2A in n_term_junc_dict[c_term_dhr]:
                    tmp_comp2 = junction2A.to_component()
                    tmp_comp2.n_term_res_trim+=junction_residues
                    rd2_design=copy.deepcopy(tmp_design)
                    rd2_design.append_component(tmp_comp2,False,"",junction_residues)
                    designs.append(rd2_design)

    for design in designs:
        #name = "{}_{}".format(options.output_structure_name,tot_job_ct)
        name = design.get_informative_name()
        m = design.output_string()
        f.write("{} {}\n".format(name,m))
        tot_job_ct+=1
        f.flush()
    f.close()



def generate_structures(design,n_term_junc_dict,c_term_junc_dict,dhrs_dict,options,queue,terminal_pdb,output_wo_extension,recursive):
    #step1 Generate terminal structures & output
    for output_design in generate_all_dhr_attach_permutations(design,dhrs_dict,options.dhr_length_min,options.dhr_length_max,options.dhr_length_min,options.dhr_length_max,output_wo_extension):
        if(terminal_pdb==None):
            if(output_design.output_valid(options)):
                queue.put(output_design.output_string() +"\n")
        else:
            if((output_design.get_current_n_term_dhr() == terminal_pdb.c_term_dhr) or (output_design.get_current_c_term_dhr() == terminal_pdb.n_term_dhr)):
                for terminal_pdb_output in generate_all_terminal_pdb_permutations(output_design,terminal_pdb,dhrs_dict,options):
                    if(terminal_pdb_output.output_valid(options)):
                        queue.put(terminal_pdb_output.output_string() +"\n")
    #step2 Add non-terminal extension & recursively call if appropriate. (First iteration not recursive so I can parallelize that stage)
    if(design.attach_valid(options,1)):
        for next_round_design in generate_all_junction_attach_permutations(design,dhrs_dict,n_term_junc_dict,c_term_junc_dict,options):
            if(recursive):
                if(next_round_design.attach_valid(options,0)):
                    for next_rd in generate_structures(next_round_design,n_term_junc_dict,c_term_junc_dict,dhrs_dict,options,queue,terminal_pdb,output_wo_extension,recursive):
                        yield(next_rd)
            else:
                yield(next_round_design)

Pool_data=collections.namedtuple('Pool_data', 'design n_term_junc_dict c_term_junc_dict dhrs_dict options terminal_pdb output_wo_extension recursive')
Listener_data=collections.namedtuple('Pool_data','q options')

def listener(q,options):
#def listener(q):
    legos = set()
    print("listener launched")
    #Listens to messages on the q, writes to files
    tot_job_ct = 1
    job_ct = 0
    fn_ct = 1
    max_job_ct = options.output_jobs_per_file
    fn_start = ""
    if(options.output_tmp_directory == None):
        fn_start = options.output_file_name
    else:
        fn_start = "{}/{}".format(options.output_tmp_directory,options.output_file_name)
    fn = "{}_{}.jobs".format(fn_start,fn_ct)
    f = open(fn, 'w')
    while 1:
        m = q.get()
        if m == 'kill':
            print("listener killed")
            break
        if(job_ct >= max_job_ct):
            fn_ct+=1
            job_ct=0
            fn = "{}_{}.jobs".format(fn_start,fn_ct)
            f.close()
            f = open(fn, 'w')
        tmp_m = m.replace("True","False")
        if(tmp_m not in legos):
            legos.add(tmp_m)
            name = "{}_{}".format(options.output_structure_name,tot_job_ct)
            f.write("{} {}".format(name,m))
            job_ct+=1
            tot_job_ct +=1
            f.flush()
    f.close()

def generate_structures_wrapper(data,q):
    print("process id = {}".format(os.getpid()))
    for tmp_design in generate_structures(data.design,data.n_term_junc_dict,data.c_term_junc_dict,data.dhrs_dict,data.options,q,data.terminal_pdb,data.output_wo_extension,data.recursive):
        pass

if __name__=='__main__':
    parser = ArgumentParser( description="" )
    parser.add_argument("-m", "--machine",type=str,default="jojo",help="computer to run on hyak,jojo or digs")
    parser.add_argument("-s", "--start_pdb",type=str,default="all",help="pdb or all")
    parser.add_argument("--terminal_pdb",type=str,default="-",help="last connection")
    parser.add_argument("--terminal_pdb_attachment_dhr",type=str,default="-",help="DHR of N or C term. WARNING:Could fail if start structure allowed to expand in both directions" )
    parser.add_argument("--chain",type=str,default="A",help="chain of start pdb")
    parser.add_argument("--start_pdb_seq_compatible",type=bool,default=False,help="does the starting pdb have a sequence compatible junction")
    parser.add_argument("--terminal_pdb_seq_compatible",type=bool,default=False,help="does the terminal pdb have a sequence compatible junction")
    parser.add_argument("-n", "--n_term_dhr",type=str,default=None, help="n_term_dhr")
    parser.add_argument("-c", "--c_term_dhr",type=str,default=None, help="c_term_dhr")
    parser.add_argument("--terminal_dhr_trim",type=int,default = 0, help="how many repeats to trim of the terminal pdb")
    parser.add_argument("--terminal_res_adjust",type=int,default = 0, help="how many residues to trim of the terminal pdb")
    parser.add_argument("--n_term_dhr_trim",type=int,default = 0, help="how many repeats to trim the start pdb on n_term")
    parser.add_argument("--c_term_dhr_trim",type=int,default = 0, help="how many repeats to trim the start pdb on c_term")
    parser.add_argument("--n_term_res_adjust",type=int,default = 0, help="how many residues to trim the start pdb on n_term")
    parser.add_argument("--c_term_res_adjust",type=int,default = 0, help="how many residues to trim the start pdb on c_term")
    parser.add_argument("--min_numb_junctions",type=int,default = 0, help="min number of junctions to output structure")
    parser.add_argument("--max_numb_junctions",type=int,default = 5, help="how many junctions total")
    parser.add_argument("--max_length",type=int,default=275,help="max number of residues allowed")
    parser.add_argument("--min_length",type=int,default=0,help="min number of residues allowed")
    parser.add_argument("--dhr_length_min",type=int,default = 1, help="min number of repeats between junctions.")
    parser.add_argument("--dhr_length_max",type=int,default = 5, help="max number of repeats between junctions.")
    #parser.add_argument("--max_jobs_per_file",type=int,default = 10000, help="max number of jobs in each output file")
    parser.add_argument("--output_structure_name",type=str,default="lego",help="short name used to identify a set of attachments")
    parser.add_argument("--output_file_name",type=str,default ="repeat_blueprint", help="output filename")
    parser.add_argument("--output_tmp_directory",type=str,default=None,help ="repeat blueprints stored here during parallel stage")
    parser.add_argument("--output_jobs_per_file",type=int,default=10000,help="number of jobs per output file")
    parser.add_argument("--avoid_duplicates",type=bool,default=True,help="avoid duplicates when the only different is start location. Bad memory implications")
    parser.add_argument("--n_core",type=int,default=20, help="number of cores to use during parallel stage")
    parser.add_argument("--enumerate_tranform_db",type=bool,default=False,help="output two junction possibiities for building transform possibiities")
    options = parser.parse_args()
    current_machine = "jojo" #expects this to be set before running script
    DHR_location_jojo = "/work/brunette/DBs/repeats/DHRs"
    #Junction_location_jojo = "/home/brunette/experiments/junction_assemblies/one_junction/test"
    Junction_location_jojo="/home/brunette/experiments/junction_assemblies/oligomers/pentamer/bex2C5_G2/rd6_selected_blueprints/selected_junctions"
    #Junction_location_jojo = "/home/brunette/DBs/repeats/junctions_verified/db_format"
    #Junction_location_jojo = "/work/brunette/DBs/repeats/junctions"
    DHR_location_hyak = "/gscratch/baker/brunette/DBs/repeats/DHRs"
    Junction_location_hyak = "/gscratch/baker/brunette/DBs/repeats/junctions"
    DHR_location_digs = "/work/brunette/DBs/repeats/DHRs"
    Junction_location_digs = "/work/brunette/DBs/repeats/junctions"
    DHR_location_short = ""
    Junction_location_short = ""
    dhr_location = DHR_location_jojo
    junction_location = Junction_location_jojo
    if(current_machine == "digs"):
        dhr_location = DHR_location_digs
        junction_location = Junction_location_digs
    if(current_machine == "hyak"):
        dhr_location = DHR_location_hyak
        junction_location = Junction_location_hyak
    dhr_run_location = DHR_location_jojo
    junction_run_location =  Junction_location_jojo
    if(options.machine == "digs"):
        dhr_run_location = DHR_location_digs
        junction_run_location = Junction_location_digs
    if(options.machine == "hyak"):
        dhr_run_location = DHR_location_hyak
        junction_run_location = Junction_location_hyak
    if(options.machine == "short"):
        dhr_run_location = DHR_location_short
        junction_run_location = Junction_location_short
    #prep output directory
    if(options.output_tmp_directory != None):
        if(os.path.exists(options.output_tmp_directory)):
            shutil.rmtree(options.output_tmp_directory)
        os.mkdir(options.output_tmp_directory)
    manager = mp.Manager()
    q = manager.Queue()
    pool = mp.Pool( processes = options.n_core )
    #put listener to work first
    watcher = pool.apply_async(listener, (q,options,) )
    #test for logical errors in options
    if(options.n_term_res_adjust>0 and options.n_term_dhr_trim == 0):
        sys.exit("if you are deleting {} n_term residues you must also be trimming at least 1 n_term_dhr".format(options.n_term_res_adjust))
    if(options.c_term_res_adjust>0 and options.c_term_dhr_trim == 0):
        sys.exit("if you are deleting {} c_term residues you must also be trimming at least 1 n_term_dhr".format(options.c_term_res_adjust))
    init()
    generate_dhrs_txt_file(dhr_location)
    generate_junctions_txt_file(junction_location,options)
    dhrs_dict = get_all_dhrs(dhr_run_location)
    n_term_junc_dict,c_term_junc_dict = get_all_junctions(junction_run_location,dhrs_dict)
    if(options.enumerate_tranform_db):
        enumerate_tranform_db(dhrs_dict,n_term_junc_dict,c_term_junc_dict,options)
        quit()
    terminal_pdb = None
    if(not options.terminal_pdb == "-"):
        terminal_pdb = get_terminal_pdb(options,dhrs_dict)
    designs_rd1 = generate_starting_round_structures(n_term_junc_dict,dhrs_dict,options)
    designs_rd2 = []
    for design in designs_rd1:
        for tmp_design in  generate_structures(design,n_term_junc_dict,c_term_junc_dict,dhrs_dict,options,q,terminal_pdb,False,False):
            tmp_data = Pool_data(tmp_design,n_term_junc_dict,c_term_junc_dict,dhrs_dict,options,terminal_pdb,True,True)
            designs_rd2.append(tmp_data)
    pool.starmap(generate_structures_wrapper,zip(designs_rd2,repeat(q)))
    q.put("kill")
    pool.close()
    pool.join()


#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   /GUIs/pyrosetta_toolkit/window_modules/ligand_ncaa_ptm_manager/ligand_ncaa_ptm_manager.py
## @brief  Window for exploring + enabling Ligands/NCAA/PTMs.
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Rosetta Imports
from rosetta import *
from rosetta.core.pose import add_variant_type_to_pose_residue

#Python Imports
import os
import glob
import re

#Tkinter Imports
from Tkinter import *
import tkFont
import tkSimpleDialog

#Toolkit Imports
from app.pyrosetta_toolkit.window_main import global_variables
from app.pyrosetta_toolkit.modules.tools import output as output_tools
from app.pyrosetta_toolkit.modules.protocols.DesignProtocols import DesignProtocols
from app.pyrosetta_toolkit.modules.tools import input as input_tools
from app.pyrosetta_toolkit.modules.tools import general_tools as gen_tools
from app.pyrosetta_toolkit.window_modules.scorefunction.ScoreFxnControl import ScoreFxnControl
#from window_main.IO.GUIInput import GUIInput

class ligand_ncaa_ptm_manager:
    """
    Class for exploring Rosetta's variants, NCAA's, and ligands.  Enable ligands/NCAA + mutate to them.
    Can use outside of windows.  Just don't use the setTk and shoTk functions.  Everything is setup in init.  To use mutate function, just use getters and setters to set a particular prop first!
    """
    def __init__(self, input_class, score_class, pose):
        self.pose = pose
        self.input_class = input_class
        self.score_class = score_class

        #Ignore this.  It is for Komodo autocomplete
        if 0:
            #self.input_class = GUIInput() Cannot check import - some cyclical dependancy issue.
            self.score_class = ScoreFxnControl()

        self.param_paths = []; #[array path string]
        self.fa_std_dir = rosetta.rosetta_database_from_env()+'/chemical/residue_type_sets/fa_standard'
        self.patch_directory = self.fa_std_dir+ '/patches'
        self.patches_on = self.fa_std_dir+"/patches.txt"

        self.params_directory = self.fa_std_dir+'/residue_types'
        self.params_on = self.fa_std_dir+"/residue_types.txt"



        self.patch_type_map = dict(); #[string "type_selection"]:[array Patch_name string]
        self.patch_name_map = dict(); #[string "name_selection]:[instance PatchProperties]

        self.param_main_map = dict(); #[string "main_selection"]:[array type_selection string]
        self.param_type_map = dict(); #[string "type_selection"]:[array param_name string]
        self.param_name_map = dict(); #[string "name_selection]:[instance ParamProperties]

        self.three_letter_to_prop_map = dict(); #[string three_letter_code]:[PropertyClass instance]
        self.variant_map = dict(); #[string variant]:[array param_name string]

        self.__read_patches__()
        self.__read_params__()


        self.selections = ["patch", "ligand", "polymer"]

        #ListBox Selections
        self.current_main_selection = StringVar()
        self.current_type_selection = StringVar()
        self.current_name_selection = StringVar()

        #Variables set when name is selected
        self.current_variant = StringVar()
        self.current_three_letter = StringVar()
        self.current_default = StringVar()
        self.current_name = StringVar()

        ###Callbacks###
        self.current_main_selection.trace_variable('w', self.__main_selection_callback__)

        #self.type_selection.trace_variable('w', self.__type_selection_callback__)
        #self.name_selection.trace_variable('w', self.__name_selection_callback__)

        self.please_delete = []; #Array of paths that are deleted when the module is closed

    def __exit__(self):
        for path in self.please_delete:
            os.remove(path)

    def setTk(self, main):
        """
        Set TK objects
        """
        main.title("Ligand/NCAA/PTM Manager")
        print "\n"
        print "To permenantly enable patches and parameters for amino acids, polymers, and ligands uncomment the file in: "
        print self.patches_on
        print self.params_on

        print "It is recommended to switch at least the statistical residue-based pair potential to the coulumbic atom-based fa_elec potential (Ligands/PTM)"
        print "Consider using the mm_std scorefunction for NCAA (No DNA)"
        print "Consider using the orbitals scorefunction for DNA/RNA/NCAA"
        print "Note that this may increase run time, especially for larger selections"
        #print "In preliminary results, it has also been shown to increase rotamer recovery for cannonical AA's (unpublished/Dunbrack Lab)"
        #print "To add a variant type to a pose, use the Per Residue Control in the advanced menu."
        self.main = main
        self.main_selection = OptionMenu(self.main, self.current_main_selection, *self.selections)
        self.type_listbox = Listbox(self.main)
        self.type_listbox.bind("<ButtonRelease-1>", lambda event:self.__type_selection_callback__())
        self.name_listbox = Listbox(self.main)
        self.name_listbox.bind("<ButtonRelease-1>", lambda event:self.__name_selection_callback__())
        self.name_listbox.bind("<Double-Button-1>", lambda event :self.open_param_or_patch())
        self.label_variant = Label(self.main, text="Variant", font=tkFont.Font(weight='bold'))
        self.label_three_letter=Label(self.main, text="Code Read", font=tkFont.Font(weight='bold'))
        self.label_default = Label(self.main, text="On by Default", font=tkFont.Font(weight='bold'))

        self.show_variant = Label(self.main, textvariable = self.current_variant, justify=CENTER)
        self.show_three_letter = Label(self.main, textvariable = self.current_three_letter, justify=CENTER)
        self.show_default = Label(self.main, textvariable = self.current_default, justify=CENTER)
        self.show_name = Label(self.main, textvariable = self.current_name)
        self.add_param_to_session_button = Button(self.main, text="Enable", command = lambda: self.enable())
        self.mutate_button = Button(self.main, text="Mutate", command = lambda: self.mutate())
        self.reload_pdb_button = Button(self.main, text = "Reload PDB", command = lambda:self.input_class.choose_load_pose())


        self.e_function_label = Label(self.main, text = "Energy Function Optimization", font=tkFont.Font(weight='bold', size=13))
        self.electrostatics_label = Label(self.main, text = "Electrostatic")
        self.switch_pair_to_elec_button=Button(self.main, text = "statistical -> coulombic", command = lambda: self.set_fa_elec(True))
        self.switch_elec_to_pair_button=Button(self.main, text = "statistical <- coulombic", command = lambda: self.set_fa_elec(False))

        self.scorefunction_label = Label(self.main, text = "Score Function")
        self.switch_to_orbitals_button = Button(self.main, text = "Orbital Based: orbitals", command = lambda: self.set_orbitals())
        self.switch_to_mm_std_button = Button(self.main, text = "Molecular Mechanics Based: mm_std", command = lambda: self.set_mm_std())

        self.current_main_selection.set(self.selections[0])

        #Photo.  This is for Jason Labonte, for his help in understanding how Rosetta works with NCAA/Patches
        DesignPhoto =PhotoImage(file = os.path.dirname(os.path.abspath(__file__))+"/media/glcnac_clean.gif")
        self.Photo = Label(master=self.main, image=DesignPhoto)
        self.Photo.image = DesignPhoto

    def shoTk(self, r, c):
        """
        Show TK objects
        """

        self.main_selection.grid(row=r, column=c+1, columnspan=2); #Cannot decide to columnspan or not!
        self.type_listbox.grid(row=r+1, column=c+1, padx=15, rowspan=6); self.name_listbox.grid(row=r+1, column=c+2, rowspan=6)
        self.add_param_to_session_button.grid(row=r+2, column=c+3)
        self.mutate_button.grid(row=r+3, column=c+3)
        self.reload_pdb_button.grid(row=r+4, column=c+3);

        self.Photo.grid(row=10, column=c+3, rowspan=9, sticky=W+E, padx=10)
        self.show_name.grid(row=r+8, column=c+1, columnspan=2)
        self.label_variant.grid(row=r+9, column=c+1); self.label_three_letter.grid(row=r+9, column=c+2)

        self.show_variant.grid(row=r+10, column=c+1); self.show_three_letter.grid(row=r+10, column=c+2)

        self.label_default.grid(row=r+11, column=c+1, columnspan=2)
        self.show_default.grid(row=r+12, column=c+1, columnspan=2)




        self.e_function_label.grid(row=r+13, column=c+1, columnspan=2, sticky=W+E, pady=7)
        self.electrostatics_label.grid(row=r+14, column=c+1, columnspan=2)
        self.switch_pair_to_elec_button.grid(row=r+15, column=c+1, columnspan=2)
        self.switch_elec_to_pair_button.grid(row=r+16, column=c+1, columnspan=2)

        self.scorefunction_label.grid(row=r+17, column=c+1, columnspan=2, pady=5)
        self.switch_to_orbitals_button.grid(row=r+18, column=c+1, columnspan=2)
        self.switch_to_mm_std_button.grid(row=r+19, column=c+1, columnspan=2)


### 'Public' Getter and Setters ###

    def get_prop(self):
        """
        Gets current property objects
        """
        return self.prop

    def set_prop(self, prop):
        """
        Sets prop.  Use this before mutating.
        """
        self.prop = prop

    def get_prop_from_three_letter_code(self, code):
        return self.three_letter_to_prop_map[code]

    def get_param_main_map(self):return self.param_main_map

    def get_param_type_map(self):return self.param_type_map

    def get_param_name_map(self):return self.param_name_map

    def get_patch_type_map(self):return self.patch_type_map

    def get_patch_name_map(self):return self.patch_name_map


    def get_patch_property_class(self, name):
        """
        Returns a patch property class from name.
        """
        return self.patch_name_map[name]

    def get_patch_property(self, name, type):
        """
        Returns a string of the patch property.  IO string, variant, etc.  See below for options.
        """
        return self.patch_name_map[name].get_property(type)

    def get_param_property_class(self, name):
        """
        Returns a param property class from name.
        """
        return self.param_name_map[name]

    def get_param_property(self, name, type):
        """
        Returns a string of the param property.  IO_String, etc.  See below for options.
        """
        return self.param_name_map[name].get_property(type)

    ### Main Functions ###
    def mutate(self, resnum=False):
        if self.current_main_selection.get()=="patch":
            variant = self.current_variant.get()
            print "Mutating to "+variant
            if not resnum:
                residue = tkSimpleDialog.askstring(title = "residue", prompt="Enter the residue you wish to mutate (resNum chain)")
                if not residue: return
                resnum = self.pose.pdb_info().pdb2pose(residue.split()[1].upper(), int(residue.split()[0]))

            add_variant_type_to_pose_residue(self.pose, variant, resnum)

            task = TaskFactory.create_packer_task(self.pose)
            task.temporarily_fix_everything()
            task.restrict_to_repacking()
            task.temporarily_set_pack_residue(resnum, True)
            pack_mover = PackRotamersMover(self.score_class.score, task)
            print self.score_class.score(self.pose)
            pack_mover.apply(self.pose)
            print self.score_class.score(self.pose)

            print "Variant type added + Packed"

        else:
            #Since It seems I cannot pass a ResidueTypeSet to packer task.  If you know how to do this, please let me know or add the code here.  Thanks.
            if not self.prop.rosetta_read_state:
                print "Cannot mutate.  Please enable ncaa through removing # in the database and reload GUI.  See documentation for instructions."
                print "Make sure to have the appropriate NCAA rotamer library downloaded or created and copied to database/rotamer/ncaa_rotlibs/"
                return

            mutant = self.current_three_letter.get()
            print "Mutating to "+mutant
            residue = tkSimpleDialog.askstring(title = "residue", prompt="Enter the residue you wish to mutate (resNum:chain)")
            if not residue:return

            #Check to make sure that particular residue can be mutated to the desired mutation.
            resnum = self.pose.pdb_info().pdb2pose(residue.split(":")[1], int(residue.split(":")[0]))
            residue_object = self.pose.residue(resnum)
            set = residue_object.residue_type_set()
            if not set.has_name(self.prop.found_name.get()):
                print "Residue does not have ResidueType in it's ResidueTypeSet."
                print "Unable to mutate"
                return
            ResDic = dict()
            ResDic[residue]=["NC:"+mutant,]
            output_tools.save_resfile_w_designdic(self.pose, ResDic, global_variables.current_directory+"/temp_resfile.txt")

            task = TaskFactory.create_packer_task(self.pose)
            parse_resfile(self.pose, task, global_variables.current_directory+"/temp_resfile.txt")
            design_mover = PackRotamersMover(self.score_class.score, task)
            design_mover.apply(self.pose)
            os.remove(global_variables.current_directory+"/temp_resfile.txt")
            print self.score_class.score(self.pose)
            print "Position mutated + Packed"

    def set_mm_std(self):
        """
        Sets the mm_std scorefunction.
        """
        self.score_class.set_scorefunction('mm_std')
        print "Please cite:  P. Douglas Renfrew, Eun Jung Choi, Brian Kuhlman, 'Using Noncanonical Amino Acids in Computational Protein-Peptide Interface Design' (2011) PLoS One."

    def set_orbitals(self):
        """
        Sets the orbitals scorefunction.
        """
        self.score_class.set_scorefunction('orbitals_talaris2013')
        print "Orbitals Scorefunction created by Steven Combs - Meiler Lab, Vanderbilt University"
        print "No paper to cite as of yet."

    def set_prop(self, prop):
        """
        Sets the current propertyclass, whether it be param or patch
        """
        self.prop = prop

    def get_prop(self):
        """
        Gets the current propertyclass, whether it be param or patch
        """
        return self.prop

    def enable(self):
        """
        The original reason the window was created.
        To use outside of gui window, just use set_prop before running.
        """

        #self.prop is set from name callback.  Now we get the path.  Add it.
        #We regenerate the nonstondard_residueset every time the user adds another.
        #Since something about the residuetypeset is global in rosetta, we check that the residue has not already been loaded.
        path = self.prop.path
        self.input_class.param_paths.append(path)
        self.input_class.nonstandard_ResidueTypeSet, self.input_class.loaded_paths = input_tools.get_residuetypeset_from_path_array(self.input_class.param_paths, self.input_class.loaded_paths)
        if not self.input_class.nonstandard_ResidueTypeSet.has_name(self.prop.found_name.get()):
            print "Residue type not loaded! "


        #Just to make sure it was loaded properly.
        check = dict()
        for path in self.input_class.param_paths:
            check[path]=0

        for path in check:
            #print path
            if not self.input_class.nonstandard_ResidueTypeSet.has_name(self.prop.found_name.get()):

                print "Residue type not loaded! "+self.prop.found_name.get()


    def set_fa_elec(self, bool):
        if bool:
            fa_pair_weight = self.score_class.score.get_weight(fa_pair)
            if fa_pair_weight!=0:
                self.score_class.score.set_weight(fa_pair, 0)
                self.score_class.score.set_weight(fa_elec, weight)
                print "Setting fa_elec at previous fa_pair weight."
            elif self.score_class.score.get_weight(fa_elec) != 0:
                print "fa_elec already set in scorefunction."
                return
            else:
                self.score_class.score.set_weight(fa_elec, .49)
                print "Setting fa_elec at weight of .49  This may need manual adjustment depending on the scorefunction."
        else:
            weight = self.score_class.score.get_weight(fa_elec)
            self.score_class.score.set_weight(fa_pair, weight)
            print "Setting fa_pair at same weight as previous fa_elec weight"
        print "E function changed."

    def open_param_or_patch(self):
        """
        Opens a text editor from system for the param or patch.
        This is the double button function of the name listbox.
        """
        if gen_tools.getOS()=="Linux":
            open_command = "gedit "
            #open_command = "xdg-open "
        else:
            open_command = "open "
        if re.search(".txt", os.path.basename(self.prop.path)):
            os.system(open_command+self.prop.path+" &")
        else:

            base = os.path.basename(self.prop.path).split(".")[0]
            new_base = base+".txt"
            os.system("cp "+self.prop.path.strip("\n")+" "+os.path.dirname(self.prop.path)+"/"+new_base)
            os.system(open_command+os.path.dirname(self.prop.path)+"/"+new_base)
            #print "open "+os.path.dirname(self.prop.path)+"/"+new_base+" &"
            #os.system("rm "+os.path.dirname(self.prop.path)+"/"+new_base)
            self.please_delete.append(os.path.dirname(self.prop.path)+"/"+new_base)

    ###Callbacks###

    def __main_selection_callback__(self, name, index, mode):
        """
        Callback for when you select a main (patch, ligand, polymer)
        """

        varValue = self.current_main_selection.get()
        #Exception handling for using the module without windows (Full Control)
        try:
            self.name_listbox.delete(0, END)
            self.mutate_button.config(state=DISABLED)
            if varValue=="patch":
                self.add_param_to_session_button.config(state=DISABLED)

                self.type_listbox.delete(0, END)
                for key in sorted(self.patch_type_map):
                    self.type_listbox.insert(END, key)
            else:
                try:
                    #self.add_param_to_session_button.config(state=NORMAL)

                    self.type_listbox.delete(0, END)
                    for key in sorted(self.param_main_map[varValue]):
                        self.type_listbox.insert(END, key)
                    self.current_variant.set("")
                except AttributeError:
                    return
        except AttributeError:
            pass

    def __type_selection_callback__(self):
        """
        Callback for selecting a type in the type_listbox
        """
        self.current_type_selection.set(self.type_listbox.get(self.type_listbox.curselection()))
        if self.current_main_selection.get()=="patch":
            self.name_listbox.delete(0, END)
            for value in self.patch_type_map[self.current_type_selection.get()]:
                self.name_listbox.insert(END, value)
        else:
            self.name_listbox.delete(0, END)
            for value in self.param_type_map[self.current_type_selection.get()]:
                self.name_listbox.insert(END, value)

    def __name_selection_callback__(self):
        """
        Callback for selecting a name in the name_listbox
        """
        self.current_name.set(self.name_listbox.get(self.name_listbox.curselection()))
        self.mutate_button.config(state=NORMAL)
        if self.current_main_selection.get()=="patch":
            self.prop = self.patch_name_map[self.current_name.get()]
        else:
            self.prop = self.param_name_map[self.current_name.get()]

        if not self.prop.rosetta_read_state:
            self.add_param_to_session_button.config(state=NORMAL)
        else:
            self.add_param_to_session_button.config(state=DISABLED)

        if self.current_main_selection.get()=="patch":
            self.current_variant.set(self.prop.molecule_type.get())
        else:
            self.current_variant.set(self.prop.variant.get())
        self.current_three_letter.set(self.prop.three_letter_name.get())
        self.current_default.set(str(self.prop.rosetta_read_state))

    ###Setup###

    def __read_patches__(self):
        """
        Read patches from directory
        """

        files = glob.glob(self.patch_directory+"/*.txt")
        for path in files:
            p = PatchProperties(path, self.patches_on)
            self.__populate_type_map__(self.patch_type_map, p)
            self.__populate_name_map__(self.patch_name_map, p)
            #self.populate_variant_map(self.variant_map, p)

        #Grab for carbohydrates.  Maybe all of them will be organized soon...
        self.patch_types = os.walk(self.patch_directory).next()[1]
        for type in self.patch_types:
            dir_path = self.patch_directory+"/"+type
            files = glob.glob(dir_path+"/*.txt")
            for path in files:
                p = PatchProperties(path, self.patches_on)

                self.__populate_type_map__(self.patch_type_map, p)
                self.__populate_name_map__(self.patch_name_map, p)
                self.__populate_three_letter_map__(self.three_letter_to_prop_map, p)

    def __read_params__(self):
        """
        Read params from directory
        """
        self.param_types = os.walk(self.params_directory).next()[1]
        for type in self.param_types:
            dir_path = self.params_directory+"/"+type
            files = glob.glob(dir_path+"/*.params")
            for path in files:
                p = ParamsProperties(path, type, self.params_on)
                self.__populate_main_map__(self.param_main_map, p)
                self.__populate_type_map__(self.param_type_map, p)
                self.__populate_name_map__(self.param_name_map, p)
                self.__populate_three_letter_map__(self.three_letter_to_prop_map, p)

    def __populate_type_map__(self, dictionary, property_class):
        """
        Populate type map for use in type_listbox
        """
        if property_class.molecule_type.get()=="":
            property_class.molecule_type.set("unknown")
        if dictionary.has_key(property_class.molecule_type.get().lower()):
                dictionary[property_class.molecule_type.get().lower()].append(property_class.name.get())
        else:
            dictionary[property_class.molecule_type.get().lower()]=[]
            self.__populate_type_map__(dictionary, property_class)

    def __populate_name_map__(self, dictionary, property_class):
        """
        Populate name map for use in name_listbox
        """
        if dictionary.has_key(property_class.name.get()):
            pass
        else:
            dictionary[property_class.name.get()]=property_class

    def __populate_main_map__(self, dictionary, property_class):
        """
        Populate main_map for use in main callback - differentiate between ligand and polymer for params.
        """
        if property_class.main_type.get()=="":
            property_class.main_type.set('unknown')

        if dictionary.has_key(property_class.main_type.get().lower()):
            dictionary[property_class.main_type.get().lower()][(property_class.molecule_type.get().lower())]=0
        else:
            dictionary[property_class.main_type.get().lower()]=dict()
            dictionary[property_class.main_type.get().lower()][(property_class.molecule_type.get().lower())]=0

    def __populate_three_letter_map__(self, dictionary, property_class):
        if not property_class.three_letter_name.get():
            dictionary[property_class.three_letter_name.get()] = property_class
            #property_class.three_letter_name.set('unknown')


class Properties:
    """
    Properties Base class. Lightweight.  Could maybe use decorator instead once I understand it.
    """

    def __init__(self, path, patch_or_param_on_file):
        self.path = path
        self.basename = os.path.basename(self.path).split(".")[0]
        self.properties=dict()
        self.main_file = patch_or_param_on_file

    def get_property(self, type):
        """
        Returns specific property from self.properties dictionary.
        """
        if not self.properties:
            print "Properties not set."
            return

        return self.properties[type].get()

    def read_properties(self):
        """
        Sets properties specified in self.properties dictionary.
        """

        FILE = open(self.path)
        #Doing it manually.  Want control of variable names.
        for line in FILE:
            line = line.strip()
            lineSP = line.split()
            try:
                if self.properties.has_key(lineSP[0]):
                    self.properties[lineSP[0]].set(lineSP[1])
            except IndexError:
                continue
        FILE.close()

    def check_if_on_by_default(self):
        """
        checks if patch/param is on by default.
        """

        FILE = open(self.main_file)
        self.rosetta_read_state = False
        for line in FILE:
            line = line.strip()
            if re.search("#", line):
                pass

            if re.search(self.basename, line) and not re.search("#", line):
                self.rosetta_read_state = True
                break
        FILE.close()

class PatchProperties(Properties):
    """
    Reads patch files.  Grabs properties to show/manipulate in GUI.
    Container for properties, which are stored in StringVar objects and can be accessed through properties dictionary or get_property() function
    """

    def __init__(self, path, patch_on_path):
        Properties.__init__(self, path, patch_on_path)
        #Initialize since setattr did not work.
        self.name = StringVar()
        self.name.set(os.path.basename(path).split(".")[0])
        self.molecule_type=StringVar()
        self.three_letter_name=StringVar()
        self.protein_or_dna=StringVar()
        self.found_name = StringVar(); #Name in file
        #So that we organize by variant type.

        self.properties = {
            "NAME":self.found_name,
            "PROPERTY":self.protein_or_dna,
            "SET_IO_STRING":self.three_letter_name,
            "TYPES":self.molecule_type}

        self.read_properties()
        self.check_if_on_by_default()


class ParamsProperties(Properties):
    """
    Reads param files.  Grabs properties to show/manipulate in GUI.
    Container for properties, which are stored in StringVar objects and can be accessed through properties dictionary or get_property() function
    """
    def __init__(self, path, base_type, params_on_path):
        Properties.__init__(self, path, params_on_path)
        self.molecule_type = StringVar()
        self.molecule_type.set(base_type); #Type of the param.  Since this information is not sepecifically in the params files, we grab from the directory name.

        #Initialize
        self.name = StringVar()
        self.name.set(os.path.basename(path).split(".")[0])
        self.found_name = StringVar(); #Name in file
        self.three_letter_name=StringVar()
        self.main_type = StringVar()
        self.variant = StringVar()
        self.properties={
            "NAME":self.found_name,
            "IO_STRING":self.three_letter_name,
            #This looks like either Ligand or Polymer.
            "VARIANT":self.variant,
            "TYPE":self.main_type
        }
        self.read_properties()
        self.check_if_on_by_default()



if __name__ == '__main__':
    """
    For Testing UI.
    """
    main = Tk()
    input_class = ""; score_class=""; pose=""
    manager = ligand_ncaa_ptm_manager(input_class, score_class, pose)
    manager.setTk(main)
    manager.shoTk(0, 0)
    main.mainloop()



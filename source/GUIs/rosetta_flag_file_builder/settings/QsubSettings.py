#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   /GUIs/pyrosetta_toolkit/window_modules/rosetta_tools/settings/QsubSettings.py
## @brief  Used to set default paths and queue settings for qsub setup
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

from Tkinter import *
import os


class QsubSettings():
    """
    Window for setting up Qsub Settings.  Holds Qsub Settings. Can run by itself.
    """
    
    def __init__(self):
        self.pwd = (self.location())[0]
        self.queue = StringVar(); #Queue in Entry box
        self.current_queue = StringVar(); #Queue in 
        self.qsub_path = StringVar()
        self.maui_showq_path = StringVar()
        self.qsub_temp = StringVar()
        self.tempscripts = StringVar()
        self.QUEUE_LIST = []
        #File names
        self.path_file_name = "/QSUBPATHS.txt"
        self.queue_file_name = "/QSUBQUEUES.txt"
        
        self.load_settings()
        self.current_queue.set(self.QUEUE_LIST[0])
        

        
    def setTk(self, main):
        """
        Sets Tk window variables
        """
        self.main = main
        self.queue_label = Label(self.main, text = "Queue")
        self.queue_entry = Entry(self.main, textvariable= self.queue)

        self.queue_options = OptionMenu(self.main, self.current_queue, *self.QUEUE_LIST)
        
        self.qsub_label = Label(self.main, text = "Qsub exec  Path")
        self.qsub_entry = Entry(self.main, textvariable = self.qsub_path)
        
        self.maui_label = Label(self.main, text = "Maui showq Path")
        self.maui_entry = Entry(self.main, textvariable = self.maui_showq_path)
        
        self.qsub_temp_label = Label(self.main, text="Qsub Output Dir")
        self.qsub_temp_entry = Entry(self.main, textvariable = self.qsub_temp)
        
        self.tempscripts_label = Label(self.main, text = "Temp Shell Script Dir")
        self.tempscripts_entry = Entry(self.main, textvariable = self.tempscripts)
        
        self.queue_add_button = Button(self.main, text = "Add", command = lambda: self.add_to_queue_list())
        self.queue_rm_button =  Button(self.main, text = "Remove", command = lambda: self.rm_from_queue_list())
        
        self.save_button = Button(self.main, text = "Save Settings", command = lambda: self.save_settings())
        
        
        
    def shoTk(self, main, r, c):
        """
        Shows Tk window variables
        """
        self.main = main
        self.queue_label.grid(row = r, column = c, rowspan = 2)
        self.queue_entry.grid(row = r, column = c+1, sticky = E+W)
        self.queue_add_button.grid(row = r, column = c+2, sticky = E+W)
        self.queue_options.grid(row = r+1, column = c+1, sticky = E+W)
        self.queue_rm_button.grid(row = r+1, column=c+2, sticky = W+E)
        
        self.qsub_label.grid(row = r+2, column=c)
        self.qsub_entry.grid(row = r+2, column = c+1, sticky=W+E)
        
        self.maui_label.grid(row = r+3, column=c)
        self.maui_entry.grid(row = r+3, column = c+1, sticky=W+E)
        
        self.qsub_temp_label.grid(row=r+4, column=c)
        self.qsub_temp_entry.grid(row=r+4, column=c+1)
        
        self.tempscripts_label.grid(row=r+5, column=c)
        self.tempscripts_entry.grid(row=r+5, column=c+1)
        
        self.save_button.grid(row = r+6, column = c, columnspan = 3, sticky = W+E)
        
    def load_settings(self):
        """
        Loads settings when GUI is run.  Sets settings in self.
        """
        
        if not os.path.exists(self.pwd+self.path_file_name):
            self.qsub_path.set('/usr/local/bin/qsub')
            self.maui_showq_path.set('/usr/local/maui/bin/showq')
            
        else:
            FILE = open(self.pwd+self.path_file_name)
            for line in FILE:
                line = line.strip()
                lineSP = line.split()
                if lineSP[0]=="QSUB":
                    self.qsub_path.set(lineSP[1])
                if lineSP[0]=="MAUI_SHOWQ":
                    self.maui_showq_path.set(lineSP[1])
                if lineSP[0]=="QSUBOUTPUT":
                    self.qsub_temp.set(lineSP[1])
                if lineSP[0]=="TEMPSCRIPT":
                    self.tempscripts.set(lineSP[1])
            FILE.close()
        if not os.path.exists(self.pwd+self.queue_file_name):
            
            print "No queue list to laod.  Please set queue."
            self.QUEUE_LIST.append("")
        else:
            FILE = open(self.pwd+self.queue_file_name)
            for line in FILE:
                line = line.strip()
                self.QUEUE_LIST.append(line)
            FILE.close()
    
    def add_to_queue_list(self):
        self.QUEUE_LIST.append(self.queue.get())
        self.queue_options["menu"].add_command(label=self.queue.get(), command=lambda temp = self.queue.get(): self.queue_options.setvar(self.queue_options.cget("textvariable"), value = temp))
        self.current_queue.set(self.queue.get())
        
    def rm_from_queue_list(self):
        ind = self.QUEUE_LIST.index(self.current_queue.get())
        self.QUEUE_LIST.pop(ind)
        self.repopulate_queue_menu()
        self.current_queue.set(self.QUEUE_LIST[0])
        
    def repopulate_queue_menu(self):
        self.queue_options["menu"].delete(0, END)
        for i in self.QUEUE_LIST:
            self.queue_options["menu"].add_command(label=i, command=lambda temp = i: self.queue_options.setvar(self.queue_options.cget("textvariable"), value = temp))
    
    def save_settings(self):
        FILE = open(self.pwd+self.path_file_name, 'w')
        FILE.write("QSUB\t"+self.qsub_path.get()+"\n")
        FILE.write("MAUI_SHOWQ\t"+self.maui_showq_path.get()+"\n")
        FILE.write("QSUBOUTPUT\t"+self.qsub_temp.get()+"\n")
	FILE.write("TEMPSCRIPT\t"+self.tempscripts.get()+"\n")
	FILE.close()
        
        FILE = open(self.pwd+self.queue_file_name, "w")
        for queue in self.QUEUE_LIST:
            if queue:
		FILE.write(queue+"\n")
        FILE.close()
        print "Settings Saved.."
        
    def location(self):
        """
        Allows the script to be self-aware of it's path.
        So that it can be imported from anywhere.
        """
        
        p = os.path.abspath(__file__)
        pathSP = os.path.split(p)
        return pathSP
    
if __name__ == '__main__':
    main = Tk()
    setup = QsubSettings()
    setup.setTk(main)
    setup.shoTk(main, 0, 0)
    main.mainloop()
    

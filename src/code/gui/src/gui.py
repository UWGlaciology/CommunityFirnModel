'''
Created on Oct 6, 2014

@author: huongnvo
'''

import Tkinter
import tkFileDialog
from Tkinter import Menu
from Tkinter import *
import os
import csv
import json

class frame(Tkinter.Tk):
    
    def __init__(self,parent):
        Tkinter.Tk.__init__(self,parent)
        self.grid()
        self.parent = parent
        self.initialize()     
        
    def initialize(self):
        menubar = Menu(self)
        self.config(menu = menubar)
        filemenu = Menu(menubar, tearoff = 0)
        menubar.add_cascade(label = 'File', menu = filemenu)
        filemenu.add_command(label = "Exit", command = quit)
        filemenu.add_command(label = 'Manual', command = self.manual)
        
        openButton   = Tkinter.Button(text = "Select a data directory", height = 1, width = 20, command = self.openData)
        jsonButton   = Tkinter.Button(text = "Select a config file", height = 1, width = 20, command = self.openConfig)    
        configButton = Tkinter.Button(text = "Make a config file", height = 1, width = 20, command = self.configJson)
        saveData     = Tkinter.Button(text = "Output directory", height = 1, width = 20, command = self.saveData)
        runButton    = Tkinter.Button(text = "Run model", height = 1, width = 45, command = self.run)
        saveButton   = Tkinter.Button(text = "Save config file", height = 1, width = 45, command = self.save)
        clearButton  = Tkinter.Button(text = "Clear console", height = 1, width = 45, command = self.clear)
        
        self.text = Text(self.parent)
        self.text.insert(Tkinter.END, "Welcome to the Community Firn Model\n")
        self.text.insert(Tkinter.END, "To get started, load a data directory and a .json file\n")
        self.text.insert(Tkinter.END, "or configure your own .json file\n")
        self.text.insert(Tkinter.END, "\n")
        
        openButton.grid(row = 0, column = 0)
        jsonButton.grid(row = 1, column = 0)
        configButton.grid(row = 0, column = 1)   
        saveData.grid(row = 1, column = 1)
        
        runButton.grid(row = 2, column = 0, columnspan = 2)
        saveButton.grid(row = 3, column = 0, columnspan = 2)
        self.text.grid(row = 0, column = 2, rowspan = 35, pady = 5, padx = 5)
        clearButton.grid(row = 36, column = 2, pady = 5)
    
    def openData(self):
        data = tkFileDialog.askdirectory()
        files = os.listdir(data)
        self.write(data + " has been opened\n")
        self.write("The following files have been loaded:\n")
        for dataFile in files:
            self.write("    " + dataFile + "\n")
        self.write("\n")
            
    def openConfig(self):
        file_opt = options = {} 
        options['filetypes'] = [('all files', '.*'), ('text files', '.txt'), ('json files', '.json')]
        options['title'] = 'Choose a config file'
        config = tkFileDialog.askopenfilename(**file_opt)
        global cc
        with open(config) as f:
            json_data = f.read()
            cc = json.loads(json_data)
            self.write(config + " has been opened\n")
            self.write("Fields in .json file\n")
            self.write(cc)
            self.write("\n")
            self.write("\n")

    def write(self, txt):
        self.text.insert(Tkinter.END, str(txt))
         
    def configJson(self):
        
        self.top = Toplevel(self)
        depthEntry      = Entry(self.top)
        adChoiceEntry   = Entry(self.top)
        siteChoiceEntry = Entry(self.top)
        z_resEntry      = Entry(self.top)
        stepsEntry      = Entry(self.top)
        userDataEntry   = Entry(self.top)
        gravityEntry    = Entry(self.top)
        thermalEntry    = Entry(self.top)
        pressureEntry   = Entry(self.top)
        conzoneEntry    = Entry(self.top)
        rho0Entry       = Entry(self.top)
        tGradEntry      = Entry(self.top)
        diffuEntry      = Entry(self.top)
        gasChoiceEntry  = Entry(self.top)
        runTypeEntry    = Entry(self.top)
          
        errorBox = Text(self.top, height = 30, width = 40)
        generateJson = Tkinter.Button(self.top, text = "Select a data directory", height = 1, width = 20, command = generate)
        
        Label(self.top, text="Depth:").grid(row = 0, sticky = W)
        Label(self.top, text="Advection method: ").grid(row = 1, sticky = W)
        Label(self.top, text="Sitechoice: ").grid(row = 2, sticky = W)
        Label(self.top, text="Z resolution: ").grid(row = 3, sticky = W)
        Label(self.top, text="Steps size: ").grid(row = 4, sticky = W)
        Label(self.top, text="Users data: ").grid(row = 5, sticky = W)
        Label(self.top, text="Gravity: ").grid(row = 6, sticky = W)
        Label(self.top, text="Thermal: ").grid(row = 7, sticky = W)
        Label(self.top, text="Pressure: ").grid(row = 8, sticky = W)
        Label(self.top, text="ConZone Depth: ").grid(row = 9, sticky = W)
        Label(self.top, text="Rho_0: ").grid(row = 10, sticky = W)
        Label(self.top, text="T Grad: ").grid(row = 11, sticky = W)
        Label(self.top, text="Diffusion: ").grid(row = 12, sticky = W)
        Label(self.top, text="Gas choices: ").grid(row = 13, sticky = W)
        Label(self.top, text="Run type: ").grid(row = 14, sticky = W)
    
        depthEntry.grid(row = 0, column = 1)
        adChoiceEntry.grid(row = 1, column = 1)
        siteChoiceEntry.grid(row = 2, column =1)
        z_resEntry.grid(row = 3, column = 1)
        stepsEntry.grid(row = 4, column = 1)
        userDataEntry.grid(row = 5, column = 1)
        gravityEntry.grid(row = 6, column = 1)
        thermalEntry.grid(row = 7, column = 1)
        pressureEntry.grid(row = 8, column = 1)
        conzoneEntry.grid(row = 9, column = 1)
        rho0Entry.grid(row = 10, column = 1)
        tGradEntry.grid(row = 11, column = 1)
        diffuEntry.grid(row = 12, column = 1)
        gasChoiceEntry.grid(row = 13, column = 1)
        runTypeEntry.grid(row = 14, column = 1)
        errorBox.grid(row = 0, column = 2, rowspan = 15, pady = 5, padx = 5)
        generateJson.grid(row = 15, column = 0, columnspan = 13 )
        
        try:
            errorBox.insert(Tkinter.END, cc)      
        except NameError:
            errorBox.insert(Tkinter.END, "Config file is currently empty\n")
        
        def generate(self):
            self.top.errorBox.insert(Tkinter.END, "Config file is currently being generated")
            
            depth       = self.depthEntry.get()
            adchoice    = self.adChoiceEntry.get()
            sitechoice  = self.siteChoiceEntry.get()
            z_res       = self.z_resEntry.get()
            steps       = self.stepsEntry.get()
            userdata    = self.userDataEntry.get()
            gravity     = self.gravityEntry.get()
            thermal     = self.thermalEntry.get()
            pressure    = self.pressureEntry.get()
            conzone     = self.conzoneEntry.get()
            rho0        = self.rho0Entry.get()
            tgrad       = self.tGradEntry.get()
            diffu       = self.diffuEntry.get()
            gaschoice   = self.gasChoiceEntry.get()
            runtype     = self.runTypeEntry.get()
    
            cc["depth"] = depth
            cc["ad_method"] = adchoice
            cc["sitechoice"] = sitechoice
            cc["z_resolution"] = z_res
            cc["StepsPerYear"] = steps
            cc["UserData"] = userdata
            cc["gravity"] = gravity
            cc["thermal"] = thermal
            cc["Pressure"] = pressure
            cc["ConZoneDepth"] = conzone
            cc["rho0"] = rho0
            cc["Tgrad"] = tgrad
            cc["diffu"] = diffu
            cc["gaschoice"] = gaschoice
            cc["runtype"] = runtype
            
            self.top.destroy()
                
    def run(self):
        print "running"
        
    def save(self):
        print "saving"
#         
    def clear(self):
        self.text.delete(1.0, Tkinter.END)
    
    def saveData(self):
        print "saving data"
    
    def manual(self):
        print "print manual"
        
if __name__ == "__main__":
    app = frame(None)
    app.title("Community Firn Model")
    app.mainloop()
    
    
    

        

'''
Created on Oct 6, 2014

@author: huongnvo
'''

import Tkinter
import tkFileDialog
from tkFileDialog import askopenfile
from Tkinter import Menu

class frame(Tkinter.Tk):
    def __init__(self,parent):
        Tkinter.Tk.__init__(self,parent)
        self.parent = parent
        self.initialize()
        
    def initialize(self):
        self.grid()
        menubar = Menu(self)
        self.config(menu = menubar)
        filemenu = Menu(menubar, tearoff = 0)
        menubar.add_cascade(label='File', menu = filemenu)
        filemenu.add_command(label ="Open", command = self.openfile)
        filemenu.add_separator()
        filemenu.add_command(label = "Exit", command = quit)
        
    
    def openfile(self):
        return askopenfile(mode = 'r')

        
if __name__ == "__main__":
    app = frame(None)
    app.title("Community Firn Model")
    app.mainloop()
    

        

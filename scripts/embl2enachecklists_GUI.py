import sys, os

# IMPORTS FOR GUI
import Tkinter as tk
import tkFileDialog
from tkFileDialog import askopenfilename

# GET DIR
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'embl2enachecklists'))

# IMPORT OWN CODE
import Embl2enachecklistsMain as EMBL2ENAclMain
import MyExceptions as ME

__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>'
__copyright__ = 'Copyright (C) 2016-2018 Michael Gruenstaeudl'
__info__ = 'embl2enachecklists'
__version__ = '2018.06.27.2030'

# WINDOW
window = tk.Tk()
window.title("embl2enachecklist")
window.geometry("360x150")

# VARIABLES
outFile = tk.StringVar()
inFile = tk.StringVar()
outFile.set("")
inFile.set("")
data = ['ITS','rRNA','trnK_matK','IGS','genomic CDS']
clType = tk.StringVar()
clType.set('ITS')

# FUNCTIONS
def submit():
    try:
        EMBL2ENAclMain.embl2enachecklists(inFile.get(), outFile.get(), clType.get())
        window.destroy()
    except Exception as error:
        errorMessage(error)


def chooseFile():
    fileChooser = tk.Tk()
    fileChooser.withdraw()
    filename = askopenfilename()

    fileChooser.destroy()
    inFile.set(filename)

def errorMessage(error):
    errorWindow = tk.Tk()
    errorWindow.title("Error " + error.getErrorNumber() + ": " + error.getErrorName())
    errorWindow.geometry("500x100")

    #LABEL
    errorLabel = tk.Label(master = errorWindow, text = error.value)
    errorLabel.grid(row = 0, column = 0)

#LABEL
chooseEMBL = tk.Label(text = "EMBL File:")
chooseOutput = tk.Label(text = "Output:")
chooseCLType = tk.Label(text = "Checklist Type:")

#BUTTON
submitButton = tk.Button(master = window, text = "Submit", command = submit)
chooseFileButton = tk.Button(master = window, text = "Choose File", command = chooseFile)
closeButton = tk.Button(master = window, text="Close", command=window.quit)

#ENTRY
outputEntry = tk.Entry(master = window, textvariable = outFile)
inputEntry = tk.Entry(master = window, textvariable = inFile)

#Dropdown
dropdown = tk.OptionMenu(window, clType, *data)

# PRINTER
chooseEMBL.grid(row=0, column=0, sticky=tk.W+tk.E)
inputEntry.grid(row=0, column=1,  sticky=tk.W+tk.E)
chooseFileButton.grid(row=0, column=2,  sticky=tk.W+tk.E)
chooseOutput.grid(row=1, column=0, sticky=tk.W+tk.E)
outputEntry.grid(row=1, column=1,  sticky=tk.W+tk.E)
chooseCLType.grid(row=2, column=0, sticky=tk.W+tk.E)
dropdown.grid(row=2, column=1,  sticky=tk.W+tk.E)
submitButton.grid(row=3, column=2,  sticky=tk.W+tk.E)
closeButton.grid(row=3, column=0, sticky=tk.W+tk.E)


window.mainloop()

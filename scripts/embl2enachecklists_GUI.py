###################
# IMPORTS FOR GUI #
###################

import sys, os
import Tkinter as tk
import tkFileDialog
from tkFileDialog import askopenfilename

###########
# GET DIR #
###########

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'embl2enachecklists'))

###################
# IMPORT OWN CODE #
###################

import Embl2enachecklistsMain as EMBL2ENAclMain
import globalVariables as GlobVars
import MyExceptions as ME

__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>'
__copyright__ = 'Copyright (C) 2016-2018 Michael Gruenstaeudl'
__info__ = 'embl2enachecklists'
__version__ = '2018.06.27.2030'

##########
# WINDOW #
##########

window = tk.Tk()
window.title("embl2enachecklist")
window.geometry("360x110")

#############
# VARIABLES #
#############

outFile = tk.StringVar()
inFile = tk.StringVar()
outFile.set("sdsd")
inFile.set("/home/yanjo/Desktop/Yannick/Programs/embl2enachecklist/examples/input/matK2.embl")
data = ['ITS','ETS','rRNA','trnK_matK','IGS','genomic_CDS']
clType = tk.StringVar()
clType.set('trnK_matK')

#############
# FUNCTIONS #
#############

def submit():
    try:
        if inFile.get().split('.')[-1] == 'embl':
            EMBL2ENAclMain.embl2enachecklists(inFile.get(), "output/" + outFile.get()+ ".enachecklist", clType.get())
            if len(GlobVars.warnings) != 0:
                warningMessage(GlobVars.warnings)
        else:
            raise ME.WrongInputFile('Be shure your inputfile is in embl format')
        window.destroy()
        doneMessage()
    except Exception as error:
        try:
            errorMessage(error)
        except:
            errorMessage(ME.ErrorNotFound(error.value))


def chooseFile():
    fileChooser = tk.Tk()
    fileChooser.withdraw()
    filename = askopenfilename()

    fileChooser.destroy()
    inFile.set(filename)

#########
# LABEL #
#########

chooseEMBL = tk.Label(text = "EMBL File:")
chooseOutput = tk.Label(text = "Output:")
chooseCLType = tk.Label(text = "Checklist Type:")

##########
# BUTTON #
##########

submitButton = tk.Button(master = window, text = "Submit", command = submit)
chooseFileButton = tk.Button(master = window, text = "Choose File", command = chooseFile)
closeButton = tk.Button(master = window, text="Close", command=window.quit)

#########
# ENTRY #
#########

outputEntry = tk.Entry(master = window, textvariable = outFile)
inputEntry = tk.Entry(master = window, textvariable = inFile)

############
# Dropdown #
############

dropdown = tk.OptionMenu(window, clType, *data)

###########
# PRINTER #
###########

chooseEMBL.grid(row=0, column=0, sticky=tk.W+tk.E)
inputEntry.grid(row=0, column=1,  sticky=tk.W+tk.E)
chooseFileButton.grid(row=0, column=2,  sticky=tk.W+tk.E)
chooseOutput.grid(row=1, column=0, sticky=tk.W+tk.E)
outputEntry.grid(row=1, column=1,  sticky=tk.W+tk.E)
chooseCLType.grid(row=2, column=0, sticky=tk.W+tk.E)
dropdown.grid(row=2, column=1,  sticky=tk.W+tk.E)
submitButton.grid(row=3, column=2,  sticky=tk.W+tk.E)
closeButton.grid(row=3, column=0, sticky=tk.W+tk.E)

######################
# Additional windows #
######################

def errorMessage(error):
    errorWindow = tk.Tk()
    errorWindow.title("Error " + error.getErrorNumber() + ": " + error.getErrorName())
    errorWindow.geometry("500x100")

    #LABEL
    errorMessage= tk.Message(master = errorWindow, text = error.value, width='495')
    errorMessage.grid(row = 0, column = 0)

    #Button
    errorWindow = tk.Button(master = errorWindow, text="OK", command=errorWindow.quit)
    errorWindow.grid(row = 1, column = 0)

def warningMessage(warningList):
    warningWindow = tk.Tk()
    warningWindow.title("Warnings")
    warningWindow.geometry("500x100")

    warningString = ''
    for warning in warningList:
        warningString = warningString + warning + '\n'

    #LABEL
    warningMessage = tk.Message(master = warningWindow, text = warningString, width='495')
    warningMessage.grid(row = 0, column = 0)

    #Button
    warningButton = tk.Button(master = warningWindow, text="OK", command=warningWindow.quit)
    warningButton.grid(row = 1, column = 0)

def doneMessage():
    doneWindow = tk.Tk()
    doneWindow.title("embl2enachecklist")
    doneWindow.geometry("350x70")

    #LABEL
    doneLabel = tk.Label(master = doneWindow, text = "Done! :D\n You can find your output here: " + "output/" + outFile.get()+ ".enachecklist")
    doneLabel.grid(row = 0, column = 0)

    #Button
    doneButton = tk.Button(master = doneWindow, text="Nice", command=doneWindow.quit)
    doneButton.grid(row = 1, column = 0, sticky=tk.W+tk.E)

############
# Mainloop #
############

window.mainloop()

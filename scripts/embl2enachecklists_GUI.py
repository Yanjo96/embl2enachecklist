 ###
  #  #    # #####   ####  #####  #####  ####
  #  ##  ## #    # #    # #    #   #   #
  #  # ## # #    # #    # #    #   #    ####
  #  #    # #####  #    # #####    #        #
  #  #    # #      #    # #   #    #   #    #
 ### #    # #       ####  #    #   #    ####

import sys, os
import Tkinter as tk
import tkFileDialog
from tkFileDialog import askopenfilename, asksaveasfilename

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'embl2enachecklists'))

import Embl2enachecklistsMain as E2E
import globalVariables as GlobVars
import MyExceptions as ME

__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>'
__copyright__ = 'Copyright (C) 2016-2018 Michael Gruenstaeudl'
__info__ = 'embl2enachecklists'
__version__ = '2018.10.04, Version 0.1.0'

class GUI():

    def __init__(self):

         #     # ### #     # ######  ####### #     #
         #  #  #  #  ##    # #     # #     # #  #  #
         #  #  #  #  # #   # #     # #     # #  #  #
         #  #  #  #  #  #  # #     # #     # #  #  #
         #  #  #  #  #   # # #     # #     # #  #  #
         #  #  #  #  #    ## #     # #     # #  #  #
          ## ##  ### #     # ######  #######  ## ##


        self.window = tk.Tk()
        self.window.title("embl2enachecklist")
        self.window.geometry("460x400")
        self.window.configure(background=GlobVars.bg)

         #     #    #    ######  ###    #    ######  #       #######  #####
         #     #   # #   #     #  #    # #   #     # #       #       #     #
         #     #  #   #  #     #  #   #   #  #     # #       #       #
         #     # #     # ######   #  #     # ######  #       #####    #####
          #   #  ####### #   #    #  ####### #     # #       #             #
           # #   #     # #    #   #  #     # #     # #       #       #     #
            #    #     # #     # ### #     # ######  ####### #######  #####

        self.outFile = tk.StringVar()
        self.inFile = tk.StringVar()
        self.outFile.set("")
        self.isEnv = 'no'
        self.inFile.set("")
        self.data = GlobVars.allowed_checklists

        self.checklistHelper = {'ITS':'For ITS rDNA region. This checklist allows generic annotation of the ITS components (18S rRNA, ITS1, 5.8S rRNA, ITS2 and 28S rRNA). For annotation of the rRNA component only, please use the rRNA gene checklist.',
                                'IGS':'For intergenic spacer (IGS) sequences'
                                      + ' between neighbouring genes (e.g.,'
                                      + ' psbA-trnH IGS, 16S-23S rRNA IGS).'
                                      + ' Inclusion of the flanking genes is'
                                      + ' allowed.',
                          'trnK_matK':'For complete or partial matK gene within'
                                      + ' the chloroplast trnK gene intron.',
                                'ETS':'For submission of External Transcribed'
                                      + ' Spacer (ETS) regions of the'
                                      + ' eukaryotic rDNA transcript; a region'
                                      + ' often used to study intrageneric'
                                      + ' relationships.',
                               'rRNA':'For ribosomal RNA genes from'
                                      + ' prokaryotic, nuclear or organellar'
                                      + ' DNA. All rRNAs are considered'
                                      + ' partial.',
                        'gene_intron':'For complete or partial single gene'
                                      + ' intron.',
                        'genomic_CDS':'Not implemented yet'
                                    }

        self.clType = tk.StringVar()
        self.clType.set('ITS')
        self.checklistLabel = tk.Message(text = self.checklistHelper[self.clType.get()], width=320, background = GlobVars.bg, foreground = GlobVars.fg, anchor="nw",bd=2)
        self.clType.trace("w",self.callback)

        self.gui()

    def gui(self):

         #          #    ######  ####### #
         #         # #   #     # #       #
         #        #   #  #     # #       #
         #       #     # ######  #####   #
         #       ####### #     # #       #
         #       #     # #     # #       #
         ####### #     # ######  ####### #######


        chooseEMBL = tk.Label(text = "Input:", background = GlobVars.bg, foreground = GlobVars.fg)
        chooseOutput = tk.Label(text = "Output:", background = GlobVars.bg, foreground = GlobVars.fg)
        chooseCLType = tk.Label(text = "Checklist Type:", background = GlobVars.bg, foreground = GlobVars.fg)
        isEnv = tk.Message(text="Is your organism from an environmental/uncultured sample?", background = GlobVars.bg, foreground = GlobVars.fg, anchor="nw", width=300)
        self.helperLabel = tk.Message(master = self.window, text="", width=240, background = GlobVars.bg, foreground = GlobVars.fg)


         ######  #     # ####### ####### ####### #     #
         #     # #     #    #       #    #     # ##    #
         #     # #     #    #       #    #     # # #   #
         ######  #     #    #       #    #     # #  #  #
         #     # #     #    #       #    #     # #   # #
         #     # #     #    #       #    #     # #    ##
         ######   #####     #       #    ####### #     #

        submitButton =     tk.Button(master = self.window, text = "Submit", command = self.submit, background = GlobVars.buttonbg, foreground = GlobVars.buttonfg)
        chooseFileButton = tk.Button(master = self.window, text = "Choose File", command = self.chooseFile, background = GlobVars.buttonbg, foreground = GlobVars.buttonfg)
        chooseOutFileButton = tk.Button(master = self.window, text = "Choose File", command = self.saveFile, background = GlobVars.buttonbg, foreground = GlobVars.buttonfg)
        closeButton =      tk.Button(master = self.window, text = "Close", command = self.window.quit, background = GlobVars.buttonbg, foreground = GlobVars.buttonfg)

        self.yesCheckbutton = tk.Checkbutton(master = self.window, text = "yes", background = GlobVars.buttonbg, command=self.isEnvTrue)
        self.noCheckbutton = tk.Checkbutton(master = self.window, text = "no", background = GlobVars.buttonbg, command=self.isEnvFalse)
        self.noCheckbutton.select()

         ####### #     # ####### ######  #     #
         #       ##    #    #    #     #  #   #
         #       # #   #    #    #     #   # #
         #####   #  #  #    #    ######     #
         #       #   # #    #    #   #      #
         #       #    ##    #    #    #     #
         ####### #     #    #    #     #    #

        outputEntry = tk.Entry(master = self.window, textvariable = self.outFile, background = GlobVars.entrybg, foreground = GlobVars.entryfg)
        inputEntry = tk.Entry(master = self.window, textvariable = self.inFile, background = GlobVars.entrybg, foreground = GlobVars.entryfg)


         ######
         #     # #####   ####  #####  #####   ####  #    # #    #
         #     # #    # #    # #    # #    # #    # #    # ##   #
         #     # #    # #    # #    # #    # #    # #    # # #  #
         #     # #####  #    # #####  #    # #    # # ## # #  # #
         #     # #   #  #    # #      #    # #    # ##  ## #   ##
         ######  #    #  ####  #      #####   ####  #    # #    #

        dropdown = tk.OptionMenu(self.window, self.clType, *self.data)
         ######  ######  ### #     # ####### ####### ######
         #     # #     #  #  ##    #    #    #       #     #
         #     # #     #  #  # #   #    #    #       #     #
         ######  ######   #  #  #  #    #    #####   ######
         #       #   #    #  #   # #    #    #       #   #
         #       #    #   #  #    ##    #    #       #    #
         #       #     # ### #     #    #    ####### #     #

        chooseEMBL.place(x=10,y=10, width=110, height=40)
        inputEntry.place(x=130,y=10, width=210, height=40)
        chooseFileButton.place(x=350,y=10, width=100, height=40)

        chooseOutput.place(x=10,y=70, width=110, height=40)
        outputEntry.place(x=130,y=70, width=210, height=40)
        chooseOutFileButton.place(x=350,y=70, width=100, height=40)

        chooseCLType.place(x=10,y=130, width=110, height=40)
        dropdown.place(x=130,y=130, width=210, height=40)

        self.checklistLabel.place(x=130,y=170, width=320, height=200)

        isEnv.place(x=10,y=280,width=300,height=40)
        self.yesCheckbutton.place(x=270,y=290)
        self.noCheckbutton.place(x=330,y=290)

        closeButton.place(x=10,y=350, width=100, height=40)
        submitButton.place(x=350,y=350, width=100, height=40)
        self.helperLabel.place(x=110,y=340, width=240, height=50)

        submitButton.bind("<Enter>", self.on_enterSubmitButton)
        submitButton.bind("<Leave>", self.on_leave)

        chooseFileButton.bind("<Enter>", self.on_enterInput)
        chooseFileButton.bind("<Leave>", self.on_leave)

        chooseOutFileButton.bind("<Enter>", self.on_enterOutput)
        chooseOutFileButton.bind("<Leave>", self.on_leave)

        dropdown.bind("<Enter>", self.on_enterChecklisttype)
        dropdown.bind("<Leave>", self.on_leave)

        self.yesCheckbutton.bind("<Enter>", self.on_enterEnvSample)
        self.yesCheckbutton.bind("<Leave>", self.on_leave)

        self.noCheckbutton.bind("<Enter>", self.on_enterEnvSample)
        self.noCheckbutton.bind("<Leave>", self.on_leave)

        closeButton.bind("<Enter>", self.on_enterClose)
        closeButton.bind("<Leave>", self.on_leave)

         #     #
         ##   ##   ##   # #    # #       ####   ####  #####
         # # # #  #  #  # ##   # #      #    # #    # #    #
         #  #  # #    # # # #  # #      #    # #    # #    #
         #     # ###### # #  # # #      #    # #    # #####
         #     # #    # # #   ## #      #    # #    # #
         #     # #    # # #    # ######  ####   ####  #

        self.window.mainloop()

      #######
      #       #    # #    #  ####  ##### #  ####  #    #  ####
      #       #    # ##   # #    #   #   # #    # ##   # #
      #####   #    # # #  # #        #   # #    # # #  #  ####
      #       #    # #  # # #        #   # #    # #  # #      #
      #       #    # #   ## #    #   #   # #    # #   ## #    #
      #        ####  #    #  ####    #   #  ####  #    #  ####

    def submit(self):
        GlobVars.warnings = []
        try:
            if self.inFile.get()[-4:] == 'embl':
                E2E.embl2enachecklists(self.inFile.get(), self.outFile.get(), 'embl', self.clType.get(), self.isEnv)
                self.ready()
            elif self.inFile.get()[-2:] == 'gb':
                E2E.embl2enachecklists(self.inFile.get(), self.outFile.get(), 'gb', self.clType.get(), self.isEnv)
                self.ready()
            else:
                raise ME.WrongInputFile('Be sure your inputfile is in embl format')

        except Exception as error:
            try:
                self.errorMessage(error)
            except:
                print error
                #self.errorMessage(ME.ErrorNotFound(error.value))

    def ready(self):
        if len(GlobVars.warnings) != 0:
            self.warningMessage(GlobVars.warnings)

        self.doneMessage()

    def callback(self,*args):
        self.checklistLabel.configure(text=self.checklistHelper[self.clType.get()])

    def isEnvTrue(self):
        if self.isEnv == 'no':
            self.noCheckbutton.deselect()
            self.isEnv = 'yes'
        else:
            self.yesCheckbutton.select()

    def isEnvFalse(self):
        if self.isEnv == 'yes':
            self.yesCheckbutton.deselect()
            self.isEnv = 'no'
        else:
            self.noCheckbutton.select()

    def on_enterSubmitButton(self, event):
        self.helperLabel.configure(text="Click on submit and start the program")

    def on_enterInput(self, event):
        self.helperLabel.configure(text="Click to choose an inputfile in genbank or embl format")

    def on_enterOutput(self, event):
        self.helperLabel.configure(text="Click create an outputfile if the choosen on do not match with your wishes")

    def on_enterChecklisttype(self, event):
        self.helperLabel.configure(text="Choose a checklist type which match with your inputfile")

    def on_enterEnvSample(self, event):
        self.helperLabel.configure(text="Is your organism from an environmental/uncultured sample? yes or no")

    def on_enterClose(self, event):
        self.helperLabel.configure(text="Click to close the program")

    def on_leave(self, enter):
        self.helperLabel.configure(text="")

    def chooseFile(self):
        fileChooser = tk.Tk()
        fileChooser.withdraw()
        filename = askopenfilename()

        fileChooser.destroy()
        self.inFile.set(filename)
        self.outFile.set(filename.split("/")[-1].split(".")[0] + ".tsv")

    def saveFile(self):
        fileChooser = tk.Tk()
        fileChooser.withdraw()
        filename = asksaveasfilename(initialdir = ".",title = "Select file",filetypes = (("LibreOffice",".tsv"),("ENA checklist","*.checklist"),("Excel",".xlsx"),("all files","*.*")))

        fileChooser.destroy()
        self.outFile.set(filename)



     ###### #####  #####   ####  #####     #    # # #    # #####   ####  #    #
     #      #    # #    # #    # #    #    #    # # ##   # #    # #    # #    #
     #####  #    # #    # #    # #    #    #    # # # #  # #    # #    # #    #
     #      #####  #####  #    # #####     # ## # # #  # # #    # #    # # ## #
     #      #   #  #   #  #    # #   #     ##  ## # #   ## #    # #    # ##  ##
     ###### #    # #    #  ####  #    #    #    # # #    # #####   ####  #    #

    def errorMessage(self, error):
        errorWindow = tk.Tk()
        errorWindow.title("Error " + error.getErrorNumber() + ": " + error.getErrorName())
        errorWindow.geometry("420x100")
        errorWindow.configure(background=GlobVars.bg)

        #Label
        errorMessage= tk.Message(master = errorWindow, text = error.value, width='400',background = GlobVars.bg, foreground = GlobVars.fg)
        errorMessage.place(x=10,y=10, width=400, height=40)

        #Button
        errorWindow = tk.Button(master = errorWindow, text="OK", command=errorWindow.destroy)
        errorWindow.place(x=160, y=50, width=110, height=40)


     #    #   ##   #####  #    # # #    #  ####     #    # # #    # #####   ####  #    #
     #    #  #  #  #    # ##   # # ##   # #    #    #    # # ##   # #    # #    # #    #
     #    # #    # #    # # #  # # # #  # #         #    # # # #  # #    # #    # #    #
     # ## # ###### #####  #  # # # #  # # #  ###    # ## # # #  # # #    # #    # # ## #
     ##  ## #    # #   #  #   ## # #   ## #    #    ##  ## # #   ## #    # #    # ##  ##
     #    # #    # #    # #    # # #    #  ####     #    # # #    # #####   ####  #    #

    def warningMessage(self, warningList):
        warningWindow = tk.Tk()
        warningWindow.title("Warnings")
        height = len(warningList)*30+10+50
        warningWindow.geometry("520x" + str(height))
        warningWindow.configure(background=GlobVars.bg)

        warningString = ''
        for counter, warning in enumerate(warningList):

            #Label
            warningLabel = tk.Label(master = warningWindow, text = warning, background = GlobVars.bg, foreground = GlobVars.fg)
            warningLabel.place(x=10,y=10+counter*30, width=500, height=20)

        #Button
        warningButton = tk.Button(master = warningWindow, text="OK",command=warningWindow.destroy, background = GlobVars.buttonbg, foreground = GlobVars.buttonfg)
        warningButton.place(x=205, y=10+len(warningList)*30, width=110, height=40)

     #####   ####  #    # ######    #    # # #    # #####   ####  #    #
     #    # #    # ##   # #         #    # # ##   # #    # #    # #    #
     #    # #    # # #  # #####     #    # # # #  # #    # #    # #    #
     #    # #    # #  # # #         # ## # # #  # # #    # #    # # ## #
     #    # #    # #   ## #         ##  ## # #   ## #    # #    # ##  ##
     #####   ####  #    # ######    #    # # #    # #####   ####  #    #

    def doneMessage(self):
        doneWindow = tk.Tk()
        doneWindow.title("embl2enachecklist")
        doneWindow.geometry("420x160")
        doneWindow.configure(background=GlobVars.bg)

        #Label
        doneLabel = tk.Message(master = doneWindow, text = "Done! :D\nYou can find your output here: " + self.outFile.get(), width = 400, background = GlobVars.bg, foreground = GlobVars.fg)
        doneLabel.place(x=10, y=10, width=400, height=100)

        #Button
        doneButton = tk.Button(master = doneWindow, text="Nice", command=doneWindow.destroy, background = GlobVars.buttonbg, foreground = GlobVars.buttonfg)
        doneButton.place(x=160, y=110, width=110, height=40)



 #     #
 ##   ##   ##   # #    #
 # # # #  #  #  # ##   #
 #  #  # #    # # # #  #
 #     # ###### # #  # #
 #     # #    # # #   ##
 #     # #    # # #    #

if __name__ == '__main__':
    GUI()

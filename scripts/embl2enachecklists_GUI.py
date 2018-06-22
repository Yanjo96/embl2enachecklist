import Tkinter as tk


window = tk.Tk()
window.title("embl2enachecklist")
window.geometry("400x400")

#LABEL
title = tk.Label(text = "embl2enachecklist")
title.grid(column=0, row=0)

#BUTTON
button1 = tk.Button(text = "Submit")
button1.grid()

#ENTRY
entry1 = tk.Entry()
entry1.grid()

#Dropdown
data = ['ITS','rRNA','trnk_matK','IGS','genomic CDS']
var = tk.StringVar()
var.set('ITS')
dropdown = tk.OptionMenu(window, var, *data)
dropdown.grid()

#

window.mainloop()

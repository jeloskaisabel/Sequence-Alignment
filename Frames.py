from tkinter import *
from tkinter import ttk
from Bio import pairwise2
from GlobalAlignment import GlobalAlignment
from Alignment import Alignment
from pandastable import Table, TableModel

class ControlFrame(object):

    def __init__(self, master, result, scrollFrame):
        self.master = master
        self.result = result
        self.scrollFrame =scrollFrame

        self.first_seq = StringVar()
        self.second_seq = StringVar()
        self.match_score  = DoubleVar()
        self.mismatch_score = DoubleVar()
        self.gap_penalty = DoubleVar()
        self.score = StringVar()


        # Default scoring scheme for sequence alignment
        self.match_score.set(3.0)
        self.mismatch_score.set(1.0)
        self.gap_penalty.set(-1.0)
        self.score.set("Alignment score:")

        self.label1 = ttk.Label(self.master, text = 'Especifique dos secuencias de ADN con bases A, T, G y C')
        self.label2 = ttk.Label(self.master, text = 'Primera secuencia')
        self.label3 = ttk.Label(self.master, text = 'Segunda secuencia')
        self.label4 = ttk.Label(self.master, text = 'Match score')
        self.label5 = ttk.Label(self.master, text = 'Mismatch penalty')
        self.label6 = ttk.Label(self.master, text = 'Gap penalty')

        self.button1 = ttk.Button(self.master, text = 'Alineamiento Global', command = self.global_alignment)
        self.button2 = ttk.Button(self.master, text = 'Alineamiento Local', command = self.local_alignment)

        self.entry1 = ttk.Entry(self.master, textvariable = self.first_seq)
        self.entry2 = ttk.Entry(self.master, textvariable = self.second_seq)
        self.entry3 = ttk.Entry(self.master, textvariable = self.match_score)
        self.entry4 = ttk.Entry(self.master, textvariable = self.mismatch_score)
        self.entry5 = ttk.Entry(self.master, textvariable = self.gap_penalty)

        self.label1.grid(column = 0, row = 0, columnspan = 2)
        self.label2.grid(column = 0, row = 1, sticky = W)
        self.label3.grid(column = 0, row = 2, sticky = W)
        self.label4.grid(column = 0, row = 3, sticky = W)
        self.label5.grid(column = 0, row = 4, sticky = W)
        self.label6.grid(column = 0, row = 5, sticky = W)


        self.button1.grid(column = 2, row = 5, sticky = (N, S, E, W))
        self.button2.grid(column = 2, row = 6, sticky = (N, S, E, W))

        self.entry1.grid(column = 1, row = 1)
        self.entry2.grid(column = 1, row = 2)
        self.entry3.grid(column = 1, row = 3)
        self.entry4.grid(column = 1, row = 4)
        self.entry5.grid(column = 1, row = 5)


    def createTable(self, df, name):
        newWindow = Toplevel(self.master)
        newWindow.title(name)
        newWindow.geometry("400x400")
        newWindow.table = Table(
                    newWindow, dataframe=df,
                    showtoolbar=False,
                    showstatusbar=True,
                    editable=False)


        newWindow.table.show()

    def global_alignment(self):
        """Align two sequences globally"""
        first = self.first_seq.get()
        second = self.second_seq.get()
        ga = GlobalAlignment(first, second, self.match_score.get(), self.mismatch_score.get(), \
                        self.gap_penalty.get())
        ans = ga.globalAlignment()
        print(ans.seqA, ans.seqB, ans.score, ans.start, ans.end)
        text = pairwise2.format_alignment(ans.seqA, ans.seqB, ans.score, ans.start, ans.end)
        '''alignments = pairwise2.align.globalms(first, second, \
            self.match_score.get(), self.mismatch_score.get(), self.gap_penalty.get(), \
            self.gap_extension.get(), one_alignment_only = True)'''

        #print(alignments[0])
        #alignments[0]
        #text = pairwise2.format_alignment(*alignments[0])
        self.result.text.delete(1.0, END)
        self.result.text.insert(1.0, text)
        print(ga.dfScoringArray)
        print(ga.dfTracebackArray)
        self.createTable(ga.dfScoringArray, "Score table")
        self.createTable(ga.dfTracebackArray, "Traceback table")

        ansT = ga.ans
        self.scrollFrame.add(ansT)

        # sets the title of the
        # Toplevel widget



    def local_alignment(self):
        """Align two sequences locally"""
        '''first = self.first_seq.get()
        second = self.second_seq.get()
        alignments = pairwise2.align.localms(first, second, \
            self.match_score.get(), self.mismatch_score.get(), self.gap_opening.get(), \
            self.gap_extension.get(), one_alignment_only = True)
        text = pairwise2.format_alignment(*alignments[0])
        self.result.text.delete(1.0, END)
        self.result.text.insert(1.0, text)'''

class ResultFrame(object):

    def __init__(self, master):
        self.master = master

        self.label = ttk.Label(master, text = 'Result of sequence alignment:')
        self.text = Text(master, width = 65, height = 10)

        self.label.grid(column = 0, row = 0, sticky = W)
        self.text.grid(column = 0, row = 1, columnspan = 2)


class ScrollBar:

    # constructor
    def __init__(self, master):
        self.master = master
        self.label = ttk.Label(master, text='Result of sequence alignment:')
        v = Scrollbar(master, orient='vertical')
        v.pack(side=RIGHT, fill=Y)
        self.t = Text(master, width=65, height=15, wrap=NONE,
                 yscrollcommand=v.set)

        self.t.pack(side=TOP, fill=X)
        v.config(command=self.t.yview)

    def add(self, ans):
        for i in ans:
            self.t.insert(END, i)

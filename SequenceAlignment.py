from tkinter import *
from tkinter import ttk
from Frames import ControlFrame, ResultFrame, ScrollBar


class AppFrame(object):

    def __init__(self, master):
        self.master = master
        self.master.title("Sequence Alignment")
        self.controlFrame = ttk.Frame(self.master, padding=5, relief='solid')
        self.resultFrame = ttk.Frame(self.master, padding=5, relief='solid')
        self.scrollFrame = ttk.Frame(self.master, padding=5, relief='solid')

        self.result_frame = ResultFrame(self.resultFrame)
        self.third_frame = ScrollBar(self.scrollFrame)
        self.first_frame = ControlFrame(self.controlFrame, self.result_frame, self.third_frame)


        self.controlFrame.grid(column=0, row=0)
        self.resultFrame.grid(column=0, row=1)
        self.scrollFrame.grid(column=0, row=2)

        # Add padding to all widgets in the control frame.
        for child in self.controlFrame.winfo_children():
            child.grid_configure(padx=3, pady=3)

        # Add padding to all widgets in the result frame.
        for child in self.resultFrame.winfo_children():
            child.grid_configure(padx=3, pady=3)


def main():
    root = Tk()
    app = AppFrame(root)

    # Add padding to all widgets in the root.
    for child in root.winfo_children():
        child.grid_configure(padx=3, pady=3)
    root.mainloop()


if __name__ == '__main__':
    main()
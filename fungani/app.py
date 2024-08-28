import tkinter as tk
import tkinter.ttk as ttk
from tkinter.filedialog import askopenfile, askdirectory
import time

ws = tk.Tk()
ws.title("FungANI")
ws.geometry("400x200")


state = {}


def select_outdir():
    pathname = askdirectory()
    if pathname is not None:
        return pathname


def select_genome(target):
    filepath = askopenfile(
        mode="r", filetypes=[("Fasta Files", ["*.fasta", "*.fna", "*.fas"])]
    )
    if filepath is not None:
        state[target] = filepath


label_reference = ttk.Label(ws, text="Reference Genome", anchor="e", width=20)
label_reference.grid(row=0, column=0, padx=10)
button_reference = ttk.Button(
    ws, text="Open", command=lambda: select_genome("reference")
)
button_reference.grid(row=0, column=1)

label_test = ttk.Label(ws, text="Test Genome", anchor="e", width=20)
label_test.grid(row=1, column=0, padx=10)
button_test = ttk.Button(ws, text="Open", command=lambda: select_genome("test"))
button_test.grid(row=1, column=1)

label_percent = ttk.Label(ws, text="Genome sample size (0-100%)", anchor="e", width=20)
label_percent.grid(row=2, column=0, padx=10)
value_percent = tk.IntVar()
edit_percent = ttk.Entry(ws, textvariable=value_percent, width=8, justify="right")
value_percent.set(100)
edit_percent.grid(row=2, column=1)

label_size = ttk.Label(ws, text="Window size (bp)", anchor="e", width=20)
label_size.grid(row=3, column=0, padx=10)
value_size = tk.IntVar()
edit_size = ttk.Entry(ws, textvariable=value_size, width=8, justify="right")
value_size.set(1000)
edit_size.grid(row=3, column=1)

label_overlap = ttk.Label(ws, text="Window overlap (bp)", anchor="e", width=20)
label_overlap.grid(row=4, column=0, padx=10)
value_overlap = tk.IntVar()
edit_overlap = ttk.Entry(ws, textvariable=value_overlap, width=8, justify="right")
value_overlap.set(500)
edit_overlap.grid(row=4, column=1)

label_cpus = ttk.Label(ws, text="Number of cores", anchor="e", width=20)
label_cpus.grid(row=4, column=0, padx=10)
value_cpus = tk.IntVar()
edit_cpus = ttk.Entry(ws, textvariable=value_cpus, width=8, justify="right")
value_cpus.set(4)
edit_cpus.grid(row=4, column=1)

label_outdir = ttk.Label(ws, text="Output directory", anchor="e", width=20)
label_outdir.grid(row=5, column=0, padx=10)
button_outdir = ttk.Button(ws, text="Open", command=select_outdir)
button_outdir.grid(row=5, column=1)

label_clean = ttk.Label(ws, text="Clean intermediate files", anchor="e", width=20)
label_clean.grid(row=6, column=0, padx=10)
value_clean = tk.IntVar()
button_clean = ttk.Checkbutton(ws, variable=value_clean)
button_clean.grid(row=6, column=1)

button_quit = ttk.Button(ws, text="Quit", command=ws.destroy)
button_quit.grid(row=7, column=2)

ws.mainloop()

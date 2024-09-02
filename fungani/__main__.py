import os
import time
import tkinter as tk
from tkinter import font
import tkinter.ttk as ttk
from tkinter.filedialog import askdirectory, askopenfile
from types import SimpleNamespace

from fungani.core import main

FONT_SIZE = 12

app = tk.Tk()
app.title("FungANI")
app.geometry("480x440")
default_font = font.nametofont("TkDefaultFont")
default_font.configure(size=FONT_SIZE)
label_default_font = font.Font(size=FONT_SIZE)
label_optional_font = font.Font(size=FONT_SIZE, slant="italic")

state = {}
state["outdir"] = None


def handle_click(event):
    args = collect_values()
    tic = time.time()
    # forward mode: test -> reference
    args.mode = "fwd"
    main(args)
    # reverse mode: reference -> test
    args.mode = "rev"
    args.test, args.reference = args.reference, args.test
    success = main(args, tic)
    if success:
        value_run.set("Done")


def collect_values():
    values = dict(
        reference=state["reference"],
        test=state["test"],
        outdir=state["outdir"],
        threshold=int(edit_threshold.get()),
        percent=int(edit_percent.get()),
        size=int(edit_size.get()),
        overlap=int(edit_overlap.get()),
        cpus=int(edit_cpus.get()),
        clean="selected" in checkbox_clean.state(),
    )
    args = SimpleNamespace(**values)
    return args


def select_directory():
    global state
    pathname = askdirectory()
    if pathname is not None:
        state["outdir"] = pathname
        value_outdir.set(pathname)


def select_file(target):
    global state
    filepath = askopenfile(
        mode="r", filetypes=[("Fasta Files", ["*.fasta", "*.fna", "*.fas"])]
    )
    if filepath is not None:
        state[target] = filepath.name
        if target == "reference":
            value_reference.set(os.path.basename(filepath.name))
        if target == "test":
            value_test.set(os.path.basename(filepath.name))


label_reference = ttk.Label(app, text="Reference genome", anchor="e", width=20)
label_reference.grid(row=0, column=0, padx=10, pady=20)
button_reference = ttk.Button(
    app, text="Select", command=lambda: select_file("reference")
)
button_reference.grid(row=0, column=1)
value_reference = tk.StringVar()
edit_reference = ttk.Entry(
    app, textvariable=value_reference, state="readonly", justify="left"
)
edit_reference.grid(row=0, column=2, padx=10)

label_test = ttk.Label(app, text="Test genome", anchor="e", width=20)
label_test.grid(row=1, column=0, pady=0)
button_test = ttk.Button(app, text="Select", command=lambda: select_file("test"))
button_test.grid(row=1, column=1)
value_test = tk.StringVar()
edit_test = ttk.Entry(app, textvariable=value_test, state="readonly", justify="left")
edit_test.grid(row=1, column=2)

label_threshold = ttk.Label(app, text="ANI threshold (%)", anchor="e", width=20)
label_threshold.grid(row=2, column=0, pady=10)
value_threshold = tk.IntVar()
edit_threshold = ttk.Entry(
    app, textvariable=value_threshold, font=("", FONT_SIZE), width=8, justify="right"
)
value_threshold.set(80)
edit_threshold.grid(row=2, column=1)

label_percent = ttk.Label(app, text="Genome size (%)", anchor="e", width=20)
label_percent.grid(row=3, column=0, pady=10)
value_percent = tk.IntVar()
edit_percent = ttk.Entry(
    app, textvariable=value_percent, font=("", FONT_SIZE), width=8, justify="right"
)
value_percent.set(10)
edit_percent.grid(row=3, column=1)

label_size = ttk.Label(app, text="Window size (bp)", anchor="e", width=20)
label_size.grid(row=4, column=0, pady=10)
value_size = tk.IntVar()
edit_size = ttk.Entry(
    app, textvariable=value_size, font=("", FONT_SIZE), width=8, justify="right"
)
value_size.set(1000)
edit_size.grid(row=4, column=1)

label_overlap = ttk.Label(app, text="Window overlap (bp)", anchor="e", width=20)
label_overlap.grid(row=5, column=0, pady=10)
value_overlap = tk.IntVar()
edit_overlap = ttk.Entry(
    app, textvariable=value_overlap, font=("", FONT_SIZE), width=8, justify="right"
)
value_overlap.set(500)
edit_overlap.grid(row=5, column=1)

label_cpus = ttk.Label(app, text="Number of cores", anchor="e", width=20)
label_cpus.grid(row=6, column=0, pady=10)
value_cpus = tk.IntVar()
edit_cpus = ttk.Entry(
    app, textvariable=value_cpus, font=("", FONT_SIZE), width=8, justify="right"
)
value_cpus.set(4)
edit_cpus.grid(row=6, column=1)

label_outdir = ttk.Label(
    app, text="Output directory", font=label_optional_font, anchor="e", width=20
)
label_outdir.grid(row=7, column=0, pady=10)
button_outdir = ttk.Button(app, text="Select", command=select_directory)
button_outdir.grid(row=7, column=1)
value_outdir = tk.StringVar()
edit_outdir = ttk.Entry(
    app, textvariable=value_outdir, state="readonly", justify="left"
)
edit_outdir.grid(row=7, column=2)

label_clean = ttk.Label(app, text="Clean intermediate files", anchor="e", width=20)
label_clean.grid(row=8, column=0, padx=10)
value_clean = tk.IntVar()
checkbox_clean = ttk.Checkbutton(app, variable=value_clean)
checkbox_clean.grid(row=8, column=1)

button_run = ttk.Button(app, text="Run", command=collect_values)
button_run.grid(row=9, column=0, pady=30)
value_run = tk.StringVar()
edit_run = ttk.Entry(
    app,
    textvariable=value_run,
    width=8,
    state="readonly",
    justify="center",
)
edit_run.grid(row=9, column=1)

button_quit = ttk.Button(app, text="Quit", command=app.destroy)
button_quit.grid(row=9, column=2)

button_run.bind("<Button-1>", handle_click)

[app.grid_columnconfigure(x, weight=1) for x in range(3)]
[app.grid_rowconfigure(x, weight=1) for x in range(9)]

app.mainloop()

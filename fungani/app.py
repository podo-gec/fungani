import tkinter as tk
import tkinter.ttk as ttk
from tkinter.filedialog import askopenfile
import time

ws = tk.Tk()
ws.title("FungANI")
ws.geometry("400x200")


def open_file():
    file_path = askopenfile(
        mode="r", filetypes=[("Fasta Files", ["*.fasta", "*.fna", "*.fas"])]
    )
    if file_path is not None:
        pass


def uploadFiles():
    pb1 = ttk.Progressbar(ws, orient=ttk.HORIZONTAL, length=300, mode="determinate")
    pb1.grid(row=4, columnspan=3, pady=20)
    for i in range(5):
        ws.update_idletasks()
        pb1["value"] += 20
        time.sleep(1)
    pb1.destroy()
    ttk.Label(ws, text="File Uploaded Successfully!", foreground="green").grid(
        row=4, columnspan=3, pady=10
    )


adhar = ttk.Label(ws, text="Reference Genome")
adhar.grid(row=0, column=0, padx=10)

adharbtn = ttk.Button(ws, text="Choose File", command=lambda: open_file())
adharbtn.grid(row=0, column=1)

dl = ttk.Label(ws, text="Test Genome")
dl.grid(row=1, column=0, padx=10)

dlbtn = ttk.Button(ws, text="Choose File", command=lambda: open_file())
dlbtn.grid(row=1, column=1)

upld = ttk.Button(ws, text="Upload Files", command=uploadFiles)
upld.grid(row=2, columnspan=3, pady=10)


ws.mainloop()

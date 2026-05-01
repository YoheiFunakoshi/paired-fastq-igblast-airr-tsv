from __future__ import annotations

from pathlib import Path
import queue
import threading
import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext, ttk

from .igblast import IgBlastConfig
from .pipeline import run_paired_igblast
from .prepare import ReadTransform


class App(ttk.Frame):
    def __init__(self, master: tk.Tk) -> None:
        super().__init__(master, padding=12)
        self.master = master
        self.messages: queue.Queue[tuple[str, str]] = queue.Queue()
        self.vars: dict[str, tk.Variable] = {}
        self._build_variables()
        self._build_layout()
        self._poll_messages()

    def _build_variables(self) -> None:
        defaults: dict[str, str | int | float | bool] = {
            "r1": "",
            "r2": "",
            "out": "",
            "query_fasta": "",
            "igblastn": "igblastn",
            "germline_db_v": "",
            "germline_db_d": "",
            "germline_db_j": "",
            "auxiliary_data": "",
            "organism": "human",
            "domain_system": "imgt",
            "ig_seqtype": "Ig",
            "num_threads": 1,
            "read_selection": "both",
            "r1_orientation": "forward",
            "r2_orientation": "reverse-complement",
            "trim_left_r1": 0,
            "trim_right_r1": 0,
            "trim_left_r2": 0,
            "trim_right_r2": 0,
            "min_length": 0,
            "max_n_rate": 1.0,
            "query_name_template": "{read_id}|{read}",
            "strict_ids": True,
        }
        for key, value in defaults.items():
            if isinstance(value, bool):
                self.vars[key] = tk.BooleanVar(value=value)
            elif isinstance(value, int):
                self.vars[key] = tk.IntVar(value=value)
            elif isinstance(value, float):
                self.vars[key] = tk.DoubleVar(value=value)
            else:
                self.vars[key] = tk.StringVar(value=value)

    def _build_layout(self) -> None:
        self.grid(sticky="nsew")
        self.master.columnconfigure(0, weight=1)
        self.master.rowconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)

        row = 0
        row = self._add_file_row(row, "R1 FASTQ", "r1", "open")
        row = self._add_file_row(row, "R2 FASTQ", "r2", "open")
        row = self._add_file_row(row, "Output TSV", "out", "save")
        row = self._add_file_row(row, "Keep query FASTA", "query_fasta", "save")

        separator = ttk.Separator(self)
        separator.grid(row=row, column=0, columnspan=3, sticky="ew", pady=8)
        row += 1

        row = self._add_file_row(row, "igblastn", "igblastn", "open")
        row = self._add_file_row(row, "V DB prefix", "germline_db_v", "open")
        row = self._add_file_row(row, "D DB prefix", "germline_db_d", "open")
        row = self._add_file_row(row, "J DB prefix", "germline_db_j", "open")
        row = self._add_file_row(row, "Aux file", "auxiliary_data", "open")

        row = self._add_entry_row(row, "Organism", "organism")
        row = self._add_entry_row(row, "Domain system", "domain_system")
        row = self._add_entry_row(row, "Seq type", "ig_seqtype")
        row = self._add_spin_row(row, "Threads", "num_threads", 1, 128)

        separator = ttk.Separator(self)
        separator.grid(row=row, column=0, columnspan=3, sticky="ew", pady=8)
        row += 1

        row = self._add_combo_row(row, "Reads", "read_selection", ("both", "r1", "r2"))
        row = self._add_combo_row(row, "R1 orientation", "r1_orientation", ("forward", "reverse-complement"))
        row = self._add_combo_row(row, "R2 orientation", "r2_orientation", ("forward", "reverse-complement"))
        row = self._add_spin_row(row, "Trim left R1", "trim_left_r1", 0, 10000)
        row = self._add_spin_row(row, "Trim right R1", "trim_right_r1", 0, 10000)
        row = self._add_spin_row(row, "Trim left R2", "trim_left_r2", 0, 10000)
        row = self._add_spin_row(row, "Trim right R2", "trim_right_r2", 0, 10000)
        row = self._add_spin_row(row, "Min length", "min_length", 0, 100000)
        row = self._add_entry_row(row, "Max N rate", "max_n_rate")
        row = self._add_entry_row(row, "Query name", "query_name_template")

        strict_ids = ttk.Checkbutton(self, text="Require matching R1/R2 IDs", variable=self.vars["strict_ids"])
        strict_ids.grid(row=row, column=1, columnspan=2, sticky="w", pady=3)
        row += 1

        self.run_button = ttk.Button(self, text="Run", command=self._start_run)
        self.run_button.grid(row=row, column=1, sticky="w", pady=8)
        row += 1

        self.log = scrolledtext.ScrolledText(self, height=10, width=80)
        self.log.grid(row=row, column=0, columnspan=3, sticky="nsew")
        self.rowconfigure(row, weight=1)

    def _add_file_row(self, row: int, label: str, key: str, mode: str) -> int:
        ttk.Label(self, text=label).grid(row=row, column=0, sticky="w", padx=(0, 8), pady=3)
        ttk.Entry(self, textvariable=self.vars[key]).grid(row=row, column=1, sticky="ew", pady=3)
        ttk.Button(self, text="Browse", command=lambda: self._browse(key, mode)).grid(row=row, column=2, padx=(8, 0))
        return row + 1

    def _add_entry_row(self, row: int, label: str, key: str) -> int:
        ttk.Label(self, text=label).grid(row=row, column=0, sticky="w", padx=(0, 8), pady=3)
        ttk.Entry(self, textvariable=self.vars[key]).grid(row=row, column=1, columnspan=2, sticky="ew", pady=3)
        return row + 1

    def _add_spin_row(self, row: int, label: str, key: str, start: int, end: int) -> int:
        ttk.Label(self, text=label).grid(row=row, column=0, sticky="w", padx=(0, 8), pady=3)
        ttk.Spinbox(self, from_=start, to=end, textvariable=self.vars[key], width=12).grid(
            row=row,
            column=1,
            sticky="w",
            pady=3,
        )
        return row + 1

    def _add_combo_row(self, row: int, label: str, key: str, values: tuple[str, ...]) -> int:
        ttk.Label(self, text=label).grid(row=row, column=0, sticky="w", padx=(0, 8), pady=3)
        combo = ttk.Combobox(self, textvariable=self.vars[key], values=values, state="readonly")
        combo.grid(row=row, column=1, sticky="w", pady=3)
        return row + 1

    def _browse(self, key: str, mode: str) -> None:
        if mode == "save":
            path = filedialog.asksaveasfilename()
        else:
            path = filedialog.askopenfilename()
        if path:
            self.vars[key].set(path)

    def _start_run(self) -> None:
        self.run_button.configure(state="disabled")
        self._log("Starting IgBLAST run...")
        thread = threading.Thread(target=self._run_pipeline, daemon=True)
        thread.start()

    def _run_pipeline(self) -> None:
        try:
            extra_query_fasta = self.vars["query_fasta"].get().strip()
            result = run_paired_igblast(
                r1_path=self.vars["r1"].get(),
                r2_path=self.vars["r2"].get(),
                output_tsv=self.vars["out"].get(),
                query_fasta=extra_query_fasta or None,
                igblast_config=IgBlastConfig(
                    igblastn=self.vars["igblastn"].get(),
                    germline_db_v=self.vars["germline_db_v"].get(),
                    germline_db_d=self.vars["germline_db_d"].get() or None,
                    germline_db_j=self.vars["germline_db_j"].get(),
                    auxiliary_data=self.vars["auxiliary_data"].get() or None,
                    organism=self.vars["organism"].get(),
                    domain_system=self.vars["domain_system"].get(),
                    ig_seqtype=self.vars["ig_seqtype"].get(),
                    num_threads=int(self.vars["num_threads"].get()),
                ),
                read_selection=self.vars["read_selection"].get(),
                r1_transform=ReadTransform(
                    self.vars["r1_orientation"].get(),
                    int(self.vars["trim_left_r1"].get()),
                    int(self.vars["trim_right_r1"].get()),
                ),
                r2_transform=ReadTransform(
                    self.vars["r2_orientation"].get(),
                    int(self.vars["trim_left_r2"].get()),
                    int(self.vars["trim_right_r2"].get()),
                ),
                min_length=int(self.vars["min_length"].get()),
                max_n_rate=float(self.vars["max_n_rate"].get()),
                query_name_template=self.vars["query_name_template"].get(),
                strict_ids=bool(self.vars["strict_ids"].get()),
            )
        except Exception as exc:
            self.messages.put(("error", str(exc)))
            return

        stats = result.stats
        message = (
            f"Done: {Path(result.output_tsv)}\n"
            f"pairs={stats.total_pairs}, records={stats.records_written}, "
            f"R1={stats.r1_written}, R2={stats.r2_written}, "
            f"skipped_short={stats.skipped_too_short}, skipped_N={stats.skipped_n_rate}"
        )
        self.messages.put(("done", message))

    def _poll_messages(self) -> None:
        try:
            while True:
                kind, message = self.messages.get_nowait()
                if kind == "error":
                    self._log("ERROR: " + message)
                    messagebox.showerror("AIRR IgBLAST", message)
                    self.run_button.configure(state="normal")
                elif kind == "done":
                    self._log(message)
                    messagebox.showinfo("AIRR IgBLAST", message)
                    self.run_button.configure(state="normal")
        except queue.Empty:
            pass
        self.after(100, self._poll_messages)

    def _log(self, message: str) -> None:
        self.log.insert("end", message + "\n")
        self.log.see("end")


def main() -> None:
    root = tk.Tk()
    root.title("AIRR IgBLAST GUI")
    root.geometry("900x760")
    App(root)
    root.mainloop()

from __future__ import annotations

from pathlib import Path
import queue
import shutil
import subprocess
import threading
import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext, ttk

from .igblast import IgBlastConfig
from .naming import (
    default_data_folder,
    default_output_tsv_path,
    default_query_fasta_path,
    default_results_folder,
)
from .pipeline import default_work_dir, run_paired_igblast
from .prepare import ReadTransform


def _has_blast_db(prefix: Path) -> bool:
    return Path(str(prefix) + ".nsq").exists()


def _find_preferred_refdata_root() -> Path | None:
    desktop = Path.home() / "Desktop"
    data_folder = desktop / "RG Paired Fastq IgBLAST AIRR tsv"
    legacy_data_folder = desktop / "Paired Fastq IgBLAST AIRR tsv"
    local_refdata = data_folder / "refdata" / "IgBlast_refdata_edit_imgt"
    legacy_local_refdata = legacy_data_folder / "refdata" / "IgBlast_refdata_edit_imgt"
    desktop_refdata = desktop / "IgBlast_refdata_edit_imgt"
    nested_refdata = desktop / "大切なフォルダ レパトア解析" / "IgBlast_refdata_edit_imgt"

    if not desktop_refdata.exists() and nested_refdata.exists():
        subprocess.run(
            ["cmd", "/c", "mklink", "/J", str(desktop_refdata), str(nested_refdata)],
            capture_output=True,
            text=True,
            check=False,
        )

    candidates = [
        local_refdata,
        legacy_local_refdata,
        desktop_refdata,
        nested_refdata,
    ]
    for root in candidates:
        if (
            _has_blast_db(root / "db" / "IMGT_IGHV.imgt")
            and _has_blast_db(root / "db" / "IMGT_IGHD.imgt")
            and _has_blast_db(root / "db" / "IMGT_IGHJ.imgt")
            and (root / "optional_file" / "human_gl.aux").exists()
        ):
            return root
    return None


class App(ttk.Frame):
    def __init__(self, master: tk.Tk) -> None:
        super().__init__(master, padding=12)
        self.master = master
        self.messages: queue.Queue[tuple[str, str]] = queue.Queue()
        self.vars: dict[str, tk.Variable] = {}
        self.data_folder = default_data_folder()
        self.results_folder = default_results_folder(self.data_folder)
        self._auto_output_path = ""
        self._auto_query_fasta_path = ""
        self._output_overridden = False
        self._query_fasta_overridden = False
        self._setting_auto_path = False
        self._build_variables()
        self._build_layout()
        self._attach_path_traces()
        self._update_default_result_paths()
        self._poll_messages()

    def _build_variables(self) -> None:
        data_folder = self.data_folder
        imgt_vdj = Path.home() / "Desktop" / "IgWork" / "IMGT_VDJ"
        refdata_root = _find_preferred_refdata_root()
        igblast_root = Path("C:/Program Files/NCBI/igblast-1.21.0")
        igblastn = shutil.which("igblastn") or str(igblast_root / "bin" / "igblastn.exe")
        aux_file = (
            refdata_root / "optional_file" / "human_gl.aux"
            if refdata_root
            else igblast_root / "optional_file" / "human_gl.aux"
        )

        if refdata_root:
            germline_db_v = refdata_root / "db" / "IMGT_IGHV.imgt"
            germline_db_d = refdata_root / "db" / "IMGT_IGHD.imgt"
            germline_db_j = refdata_root / "db" / "IMGT_IGHJ.imgt"
        else:
            germline_db_v = imgt_vdj / "human_IGHV_IMGT"
            germline_db_d = imgt_vdj / "human_IGHD_IMGT"
            germline_db_j = imgt_vdj / "human_IGHJ_IMGT"

        defaults: dict[str, str | int | float | bool] = {
            "r1": "",
            "r2": "",
            "out": "",
            "query_fasta": "",
            "igblastn": igblastn,
            "germline_db_v": str(germline_db_v) if _has_blast_db(germline_db_v) else "",
            "germline_db_d": str(germline_db_d) if _has_blast_db(germline_db_d) else "",
            "germline_db_j": str(germline_db_j) if _has_blast_db(germline_db_j) else "",
            "auxiliary_data": str(aux_file) if aux_file.exists() else "",
            "organism": "human",
            "domain_system": "imgt",
            "ig_seqtype": "Ig",
            "num_threads": 4,
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
        row = self._add_db_row(row, "V DB prefix", "germline_db_v")
        row = self._add_db_row(row, "D DB prefix", "germline_db_d")
        row = self._add_db_row(row, "J DB prefix", "germline_db_j")
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

    def _add_db_row(self, row: int, label: str, key: str) -> int:
        ttk.Label(self, text=label).grid(row=row, column=0, sticky="w", padx=(0, 8), pady=3)
        ttk.Entry(self, textvariable=self.vars[key]).grid(row=row, column=1, sticky="ew", pady=3)
        ttk.Button(self, text="Browse", command=lambda: self._browse_db_prefix(key)).grid(
            row=row,
            column=2,
            padx=(8, 0),
        )
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
            path = filedialog.asksaveasfilename(**self._save_dialog_options(key))
        else:
            path = filedialog.askopenfilename(**self._open_dialog_options(key))
        if path:
            self.vars[key].set(path)
            if key in {"r1", "r2"}:
                self._update_default_result_paths()

    def _open_dialog_options(self, key: str) -> dict[str, object]:
        options: dict[str, object] = {}
        if key in {"r1", "r2"} and self.data_folder.exists():
            options["initialdir"] = str(self.data_folder)
            options["filetypes"] = [
                ("FASTQ files", "*.fastq *.fastq.gz *.fq *.fq.gz"),
                ("All files", "*.*"),
            ]
        return options

    def _save_dialog_options(self, key: str) -> dict[str, object]:
        options: dict[str, object] = {}
        if key == "out":
            suggested = self._suggested_output_path()
            options["defaultextension"] = ".tsv"
            options["filetypes"] = [("TSV files", "*.tsv"), ("All files", "*.*")]
        elif key == "query_fasta":
            suggested = self._suggested_query_fasta_path()
            options["defaultextension"] = ".fasta"
            options["filetypes"] = [("FASTA files", "*.fasta *.fa"), ("All files", "*.*")]
        else:
            suggested = None

        initialdir = suggested.parent if suggested else self.results_folder
        try:
            initialdir.mkdir(parents=True, exist_ok=True)
        except OSError:
            pass
        options["initialdir"] = str(initialdir)
        if suggested:
            options["initialfile"] = suggested.name
        return options

    def _attach_path_traces(self) -> None:
        self.vars["r1"].trace_add("write", lambda *_: self._update_default_result_paths())
        self.vars["r2"].trace_add("write", lambda *_: self._update_default_result_paths())
        self.vars["out"].trace_add("write", lambda *_: self._mark_output_overridden())
        self.vars["query_fasta"].trace_add("write", lambda *_: self._mark_query_fasta_overridden())

    def _mark_output_overridden(self) -> None:
        if self._setting_auto_path:
            return
        value = str(self.vars["out"].get()).strip()
        self._output_overridden = bool(value and value != self._auto_output_path)

    def _mark_query_fasta_overridden(self) -> None:
        if self._setting_auto_path:
            return
        value = str(self.vars["query_fasta"].get()).strip()
        self._query_fasta_overridden = bool(value and value != self._auto_query_fasta_path)

    def _suggested_output_path(self) -> Path | None:
        r1 = str(self.vars["r1"].get()).strip()
        r2 = str(self.vars["r2"].get()).strip()
        if not r1:
            return None
        return default_output_tsv_path(r1, r2 or None, self.data_folder)

    def _suggested_query_fasta_path(self) -> Path | None:
        r1 = str(self.vars["r1"].get()).strip()
        r2 = str(self.vars["r2"].get()).strip()
        if not r1:
            return None
        return default_query_fasta_path(r1, r2 or None, self.data_folder)

    def _update_default_result_paths(self) -> None:
        suggested_output = self._suggested_output_path()
        suggested_query = self._suggested_query_fasta_path()
        if not suggested_output:
            return

        current_output = str(self.vars["out"].get()).strip()
        current_query = str(self.vars["query_fasta"].get()).strip()
        if (
            self._output_overridden
            and current_output != self._auto_output_path
            and Path(current_output).name.lower() != "result.airr.tsv"
        ):
            update_output = False
        else:
            update_output = True

        if self._query_fasta_overridden and current_query != self._auto_query_fasta_path:
            update_query = False
        else:
            update_query = suggested_query is not None

        self._setting_auto_path = True
        try:
            if update_output:
                self._auto_output_path = str(suggested_output)
                self.vars["out"].set(self._auto_output_path)
                self._output_overridden = False
            if update_query and suggested_query is not None:
                self._auto_query_fasta_path = str(suggested_query)
                self.vars["query_fasta"].set(self._auto_query_fasta_path)
                self._query_fasta_overridden = False
        finally:
            self._setting_auto_path = False

    def _browse_db_prefix(self, key: str) -> None:
        path = filedialog.askopenfilename()
        if not path:
            return
        selected = Path(path)
        db_suffixes = {
            ".ndb",
            ".nhr",
            ".nin",
            ".nog",
            ".nos",
            ".not",
            ".nsq",
            ".ntf",
            ".nto",
            ".phr",
            ".pin",
            ".pog",
            ".psd",
            ".psi",
            ".psq",
        }
        if selected.suffix.lower() in db_suffixes:
            selected = selected.with_suffix("")
        self.vars[key].set(str(selected))

    def _start_run(self) -> None:
        missing = self._missing_required_fields()
        if missing:
            messagebox.showerror(
                "AIRR IgBLAST",
                "次の項目を入力してください:\n\n" + "\n".join(f"- {name}" for name in missing),
            )
            return
        self.run_button.configure(state="disabled")
        self._log("Starting IgBLAST run...")
        thread = threading.Thread(target=self._run_pipeline, daemon=True)
        thread.start()

    def _missing_required_fields(self) -> list[str]:
        required = [
            ("R1 FASTQ", "r1"),
            ("R2 FASTQ", "r2"),
            ("Output TSV", "out"),
            ("igblastn", "igblastn"),
            ("V DB prefix", "germline_db_v"),
            ("J DB prefix", "germline_db_j"),
            ("Aux file", "auxiliary_data"),
        ]
        return [label for label, key in required if not str(self.vars[key].get()).strip()]

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
                work_dir=default_work_dir(),
            )
        except Exception as exc:
            self.messages.put(("error", str(exc)))
            return

        stats = result.stats
        pair_stats = result.pair_summary_stats
        message = (
            f"Done: {Path(result.output_tsv)}\n"
            f"pairs={stats.total_pairs}, records={stats.records_written}, "
            f"R1={stats.r1_written}, R2={stats.r2_written}, "
            f"skipped_short={stats.skipped_too_short}, skipped_N={stats.skipped_n_rate}"
        )
        if result.r1_tsv and result.r2_tsv and result.integrated_tsv:
            message += (
                f"\nR1 TSV: {result.r1_tsv}"
                f"\nR2 TSV: {result.r2_tsv}"
                f"\nIntegrated TSV: {result.integrated_tsv}"
            )
        if result.counts_tsv:
            message += f"\nCounts TSV: {result.counts_tsv}"
        if result.counts_xlsx:
            message += f"\nCounts Excel: {result.counts_xlsx}"
        if pair_stats:
            message += (
                f"\nintegrated_pairs={pair_stats.total_pairs}, "
                f"junction_aa_conflicts={pair_stats.junction_aa_conflicts}, "
                f"included_in_counts={pair_stats.included_in_counts}, "
                f"unique_final_clonotypes={pair_stats.unique_final_clonotypes}"
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
    root.title("RG Paired Fastq IgBLAST AIRR tsv")
    root.geometry("900x760")
    App(root)
    root.mainloop()

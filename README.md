# AIRR IgBLAST Paired

R1/R2 の paired-end FASTQ を **マージせず**、R1 と R2 を別々の IgBLAST query として `igblastn -outfmt 19` に渡し、AIRR Rearrangement TSV を出力するツールです。

GUI から実行することも、CLI から実行することもできます。

## 前提

- Python 3.10 以上
- `igblastn` が実行できること
- IgBLAST 用 germline DB が用意済みであること
  - 例: `human_gl_V`, `human_gl_D`, `human_gl_J`
  - DB 名は `makeblastdb` 後のファイル接尾辞を除いた prefix を指定します。
- CDR/FWR 情報が必要な場合は、IgBLAST 付属の `.aux` ファイルも指定します。

## インストール

```powershell
python -m pip install -e .
```

インストールせずに実行する場合:

```powershell
$env:PYTHONPATH = "src"
python -m airr_igblast_paired --help
```

## GUI

```powershell
airr-igblast-paired gui
```

GUI では次を指定できます。

- R1 FASTQ
- R2 FASTQ
- 出力 TSV
- IgBLAST executable
- V/D/J germline DB
- auxiliary data file
- R1/R2 の向き
- R1/R2 の左右 trim
- 最小配列長
- 最大 N 率
- FASTA query 名テンプレート

## CLI

R1/R2 から AIRR TSV を作る例:

```powershell
airr-igblast-paired run `
  --r1 sample_R1.fastq.gz `
  --r2 sample_R2.fastq.gz `
  --out sample.airr.tsv `
  --germline-db-v "C:\Program Files\NCBI\igblast-1.21.0\database\human_gl_V" `
  --germline-db-d "C:\Program Files\NCBI\igblast-1.21.0\database\human_gl_D" `
  --germline-db-j "C:\Program Files\NCBI\igblast-1.21.0\database\human_gl_J" `
  --auxiliary-data "C:\Program Files\NCBI\igblast-1.21.0\optional_file\human_gl.aux" `
  --organism human `
  --domain-system imgt `
  --num-threads 4
```

IgBLAST に渡す FASTA だけを作る例:

```powershell
airr-igblast-paired prepare `
  --r1 sample_R1.fastq.gz `
  --r2 sample_R2.fastq.gz `
  --out-fasta sample.queries.fasta
```

`run` でも `--query-fasta sample.queries.fasta` を付けると、中間 FASTA を保存できます。

## R1/R2 の扱い

- R1 と R2 はマージしません。
- FASTA query はデフォルトで `read_id|R1` と `read_id|R2` の2本になります。
- R1 はデフォルトで forward のまま使います。
- R2 はデフォルトで reverse-complement して使います。
- 出力 TSV は IgBLAST の `-outfmt 19` そのものです。

主な調整項目:

```powershell
--read-selection both
--read-selection r1
--read-selection r2
--r1-orientation forward
--r2-orientation reverse-complement
--trim-left-r1 0
--trim-right-r1 0
--trim-left-r2 0
--trim-right-r2 0
--min-length 80
--max-n-rate 0.10
--query-name-template "{read_id}|{read}"
```

標準エラーには簡単な処理統計が出ます。

```text
prepare stats: total_pairs=1000, records_written=2000, r1_written=1000, r2_written=1000, skipped_too_short=0, skipped_n_rate=0
```

## テスト

インストール後:

```powershell
python -m unittest discover -s tests
```

インストールしない場合:

```powershell
$env:PYTHONPATH = "src"
python -m unittest discover -s tests
```

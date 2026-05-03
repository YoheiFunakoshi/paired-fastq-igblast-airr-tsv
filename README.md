# RG Paired Fastq IgBLAST AIRR tsv

RG社で実施したPCRシーケンス由来のレパトア解析データを、ペアの FASTQ ファイル、つまり R1 と R2 から IgBLAST を実行して AIRR 形式の TSV ファイルにするツールです。

このツールでは R1 と R2 をマージしません。IgBLAST には `read_id|R1` と `read_id|R2` のように、R1 と R2 を別々の query として渡します。

RGデータの前提と、この解析システムの方針は次の通りです。

- RG社データはUMIなしのPCRシーケンスデータとして扱います。
- FASTQはトリムしません。標準設定では `Trim left/right R1/R2 = 0` です。
- R1/R2はマージしません。
- R1/R2を別々にIgBLASTへ渡し、pair単位の比較やCDR3/junction_aa rescueに使えるAIRR TSVを作成します。
- complete VDJやSHM解析を主目的とするものではなく、RG社レパトア解析結果との比較、特にQASAS用のV/J/junctionAA候補確認を主目的にしています。

## 基本情報

- 目的: RG社PCRシーケンス由来のpaired-end FASTQからAIRR TSVを作成すること
- 入力: R1 FASTQ と R2 FASTQ
- 中間処理: R1/R2 を別々の FASTA query に変換
- 解析: NCBI IgBLAST の `igblastn`
- 参照データ: IMGT由来のIGHV/IGHD/IGHJ配列をIgBLAST用に整形したBLAST DB
- 出力: IgBLAST `-outfmt 19` による AIRR Rearrangement TSV
- GUI: あり
- CLI: あり
- 確認済み IgBLAST: `igblastn 1.21.0`
- GUIの標準スレッド数: `4`

このリポジトリには IgBLAST 本体、IgBLAST germline database、大きな FASTQ データ、RG社Excelなどの企業解析結果は含めません。利用者のPCに IgBLAST と database を用意してから使います。

## KKF103での確認

このシステムは、KKF103データで動作と出力内容の妥当性を確認しました。

KKF103に関するFASTQ、RG社作成Excel、ChatGPTによる確認Word、比較用FASTA、AIRR TSVなどの実データはGitHubには置きません。共同研究で必要な場合は、船越がローカルに保有しているデータを別途共有します。

## 出力ファイルと統合方針

このツールの基本コンセプトは、正確さを厳しく追い込むよりも、まず情報を落とさずに残すことです。R1とR2はマージせず、それぞれを独立したIgBLAST queryとして解析します。解析後には、追跡しやすいように次のTSVを同じResultsフォルダに作成します。

```text
<sample>.airr.tsv
<sample>.R1.airr.tsv
<sample>.R2.airr.tsv
<sample>.integrated.tsv
<sample>.integrated_counts.tsv
<sample>.integrated_counts.xlsx
```

- `<sample>.airr.tsv`: IgBLAST `-outfmt 19` の生AIRR TSVです。R1とR2の行が両方入ります。
- `<sample>.R1.airr.tsv`: 生AIRR TSVからR1行だけを抜き出したAIRR TSVです。
- `<sample>.R2.airr.tsv`: 生AIRR TSVからR2行だけを抜き出したAIRR TSVです。
- `<sample>.integrated.tsv`: R1/R2をread pair単位で並べ、暫定final値を付けた追跡用TSVです。これはAIRR標準そのものではなく、共同解析で確認しやすくするためのサマリーです。全pairを残し、最終集計に入れるかどうかを `include_in_counts` と `exclude_reason` に記録します。
- `<sample>.integrated_counts.tsv`: `include_in_counts=true` のpairだけを、アリルなしV候補セット、アリルなしJ候補セット、`final_junction_aa` ごとに集計したTSVです。解析やグラフ作成ではこちらを使うと便利です。
- `<sample>.integrated_counts.xlsx`: `integrated_counts.tsv` と同じ内容のExcelファイルです。共同研究者がExcelで確認しやすいように作成します。

統合TSVでは、R1/R2の元データを消さず、`r1_*`、`r2_*`、`final_*`、`*_source`、`*_decision_reason` の列で判断過程を残します。

暫定ルール:

- `junction_aa` がR1/R2で一致する場合はその値を採用します。
- 片方だけに `junction_aa` がある場合は、その値を採用します。
- `junction_aa` が不一致の場合も落としません。長い方を `final_junction_aa` とし、同じ長さならR2を優先します。
- `v_call` はR2を優先します。R2はV領域側から読むためです。
- `j_call`、`d_call`、`productive`、`junction` は、原則として `final_junction_aa` に採用したread側を優先し、空欄の場合はもう一方を使います。
- QASASなどで使いやすいよう、`final_v_call`、`final_j_call`、`final_junction_aa` がそろう場合に `usable_for_qasas=true` とします。

この統合ルールは「落とさない」暫定版です。将来的に、conflictをより厳しく扱う正確性重視版を別に作る余地を残しています。

`integrated.tsv` では全pairを残します。ただし、次の条件をすべて満たすものだけを `include_in_counts=true` とし、`integrated_counts.tsv` と `integrated_counts.xlsx` の集計対象にします。

- `productive` がtrue相当であること
- V候補セットがあること
- J候補セットがあること
- `final_junction_aa` があること
- `final_junction_aa` にstop `*` がないこと
- `final_junction_aa` が `C` で始まること
- `final_junction_aa` が `W` または `F` で終わること
- `final_junction_aa` の長さが5-40 amino acidsであること

`integrated_counts.tsv` は、統合TSVのうち `include_in_counts=true` の行だけを集計した表です。主な列は次の通りです。

- `unique_v_gene_set`: アリルを外したV候補セット。例: `IGHV4-59,IGHV4-61`
- `unique_j_gene_set`: アリルを外したJ候補セット。例: `IGHJ4`
- `final_junction_aa`: 暫定採用したCDR3アミノ酸配列
- `read_pair_count`: 同じ `unique_v_gene_set` + `unique_j_gene_set` + `final_junction_aa` を持つread pair数
- `match_count`: R1/R2の `junction_aa` が一致したpair数
- `conflict_count`: R1/R2の `junction_aa` が不一致だったpair数
- `r1_only_count`: R1だけで `junction_aa` が得られたpair数
- `r2_only_count`: R2だけで `junction_aa` が得られたpair数
- `productive_true_count`: IgBLASTの `productive` がtrue相当だったread pair数
- `canonical_junction_aa_count`: canonical junctionAA条件を満たしたread pair数

`integrated.tsv` には `final_v_call`、`final_d_call`、`final_j_call` としてアリル付きのIgBLAST出力を残します。例えば `IGHV1-69*04` のような表記です。一方、`integrated_counts.tsv` ではRG社Excelや他社解析と比較しやすいよう、アリルを外した候補セットで集計します。D遺伝子はコールされないことが多く不安定なため、ユニーク判定には使いません。

`productive_true_count` は、RG社Excelの `In frame` と完全に同じ判定ではありません。IgBLASTの `productive` 列を使った参考指標です。RG社納品Excelに近い見方をする場合、`include_in_counts=true` の候補を主に確認します。

つまり、`integrated.tsv` は追跡用、`integrated_counts.tsv` と `integrated_counts.xlsx` はリード数集計・解析用です。

## IMGT参照データについて

IgBLASTでAIRR TSVを作成するには、V/D/Jのgermline databaseが必要です。このプロジェクトでは、IMGTから取得したヒトIgHのIGHV/IGHD/IGHJ FASTAをIgBLAST用に整形し、`makeblastdb` で作成したBLAST DBを使う想定です。

このPCでは、作業フォルダ内の次の参照データフォルダを標準として使います。

```text
C:\Users\Yohei Funakoshi\Desktop\RG Paired Fastq IgBLAST AIRR tsv\refdata\IgBlast_refdata_edit_imgt
```

このフォルダがない場合は、Desktop直下の `IgBlast_refdata_edit_imgt`、または既存の参照データフォルダも探索します。これはIgBLAST/BLASTが日本語を含む参照データパスで失敗する場合を避けるためです。

想定フォルダ構成:

```text
IgBlast_refdata_edit_imgt
├─ db
│  ├─ IMGT_IGHV.imgt.*
│  ├─ IMGT_IGHD.imgt.*
│  └─ IMGT_IGHJ.imgt.*
├─ internal_data
└─ optional_file
   └─ human_gl.aux
```

GUIはこのフォルダが見つかる場合、次の値を自動入力します。

```text
V DB prefix: ...\db\IMGT_IGHV.imgt
D DB prefix: ...\db\IMGT_IGHD.imgt
J DB prefix: ...\db\IMGT_IGHJ.imgt
Aux file:    ...\optional_file\human_gl.aux
```

IMGT FASTAのヘッダは、そのままだとIgBLASTで扱いにくい場合があります。推奨はIgBLAST付属の `edit_imgt_file.pl` でヘッダを整形してから `makeblastdb` する方法です。

例:

```powershell
perl "C:\Program Files\NCBI\igblast-1.21.0\bin\edit_imgt_file.pl" IMGT_IGHV.fasta > IMGT_IGHV.imgt.fasta
perl "C:\Program Files\NCBI\igblast-1.21.0\bin\edit_imgt_file.pl" IMGT_IGHD.fasta > IMGT_IGHD.imgt.fasta
perl "C:\Program Files\NCBI\igblast-1.21.0\bin\edit_imgt_file.pl" IMGT_IGHJ.fasta > IMGT_IGHJ.imgt.fasta
```

DB作成例:

```powershell
makeblastdb -parse_seqids -dbtype nucl -in IMGT_IGHV.imgt.fasta -out db\IMGT_IGHV.imgt
makeblastdb -parse_seqids -dbtype nucl -in IMGT_IGHD.imgt.fasta -out db\IMGT_IGHD.imgt
makeblastdb -parse_seqids -dbtype nucl -in IMGT_IGHJ.imgt.fasta -out db\IMGT_IGHJ.imgt
```

研究者間で結果を比較する場合は、IgBLASTのバージョン、IMGT参照データの取得日、ヘッダ整形方法、`makeblastdb` の条件を記録してください。

## IgBLAST について

このツールは内部で次の形式のコマンドを実行します。

```powershell
igblastn -query queries.fasta -out result.airr.tsv -outfmt 19 ...
```

`-outfmt 19` は IgBLAST が AIRR Rearrangement 形式の TSV を出力するための指定です。つまり、最終的な TSV はこのツール独自の表ではなく、IgBLAST が出力する AIRR TSV です。

開発環境で確認した IgBLAST バージョン:

```text
igblastn: 1.21.0
Package: igblast 1.21.0, build Apr 12 2023 18:52:21
```

利用者のPCでバージョンを確認するには、PowerShell で次を実行します。

```powershell
igblastn -version
```

## データ置き場

サンプル FASTQ や結果 TSV のやり取りには、Desktop に作成した次のフォルダを使っています。

```text
C:\Users\Yohei Funakoshi\Desktop\RG Paired Fastq IgBLAST AIRR tsv
```

例:

```text
sample_R1.fastq.gz
sample_R2.fastq.gz
Results of RG Paired Fastq IgBLAST AIRR tsv\sample.airr.tsv
```

非圧縮 FASTQ も使えます。

```text
sample_R1.fastq
sample_R2.fastq
```

大きい FASTQ は GitHub に載せない方針が安全です。GitHub にはソフト本体、説明文、再現手順だけを置きます。

## ローカル作業フォルダの構成

このPCでは、Desktop上の次のフォルダを「後で使うための作業フォルダ」として使います。

```text
C:\Users\Yohei Funakoshi\Desktop\RG Paired Fastq IgBLAST AIRR tsv
```

推奨構成:

```text
RG Paired Fastq IgBLAST AIRR tsv
├─ Open RG Paired Fastq IgBLAST AIRR tsv.lnk
├─ Launch RG Paired Fastq IgBLAST AIRR tsv.ps1
├─ RG Paired Fastq IgBLAST AIRR tsv.ico
├─ app
│  └─ ソフト本体のコピー
├─ refdata
│  └─ IgBlast_refdata_edit_imgt
│     ├─ db
│     ├─ internal_data
│     └─ optional_file
├─ Results of RG Paired Fastq IgBLAST AIRR tsv
│  └─ sample.airr.tsv
├─ sample_R1.fastq.gz
└─ sample_R2.fastq.gz
```

`Open RG Paired Fastq IgBLAST AIRR tsv.lnk` をダブルクリックすると、作業フォルダ内の `app` を使ってGUIを起動します。

`refdata/IgBlast_refdata_edit_imgt` にはIMGT由来のIgBLAST参照データを置きます。この参照データは解析に必要ですが、GitHubには含めません。

## GitHubに置くものと置かないもの

GitHubに置くもの:

- Pythonソースコード
- GUI/CLIの使い方
- IgBLASTとIMGT参照データの準備方法
- フォルダ構成
- テストコード

GitHubに置かないもの:

- 研究用FASTQデータ
- AIRR TSVなどの実解析結果
- RG社作成Excelなどの企業解析結果
- KKF103確認用の実データ一式
- IMGT/IgBLAST参照DBの実ファイル
- PC固有の個人フォルダにしか意味がない設定ファイル

参照データ本体はGitHubには置かず、READMEの手順に従って各利用者のPCで準備します。

KKF103でのデータ確認は実施済みです。確認に使ったFASTQ、RG社Excel、比較用FASTA、AIRR TSV、解析メモなどは船越が保有しています。

## 前提

- Python 3.10 以上
- `igblastn` が実行できること
- IgBLAST 用 germline DB が用意済みであること
  - 例: `IMGT_IGHV.imgt`, `IMGT_IGHD.imgt`, `IMGT_IGHJ.imgt`
  - DB 名は `makeblastdb` 後のファイル接尾辞を除いた prefix を指定します。
- CDR/FWR 情報が必要な場合は、IgBLAST 付属の `.aux` ファイルも指定します。

## インストール

このリポジトリを取得したフォルダで実行します。

```powershell
python -m pip install -e .
```

インストールせずに動作確認する場合:

```powershell
$env:PYTHONPATH = "src"
python -m airr_igblast_paired --help
```

## GUI の起動

```powershell
rg-paired-fastq-igblast-airr-tsv gui
```

GUI が開いたら、R1 FASTQ、R2 FASTQ、出力 TSV、IgBLAST database などを指定して `Run` を押します。

## GUI の入力欄とボタン

### FASTQ と出力ファイル

- `R1 FASTQ`
  - R1 側の FASTQ ファイルを指定します。
  - `.fastq` と `.fastq.gz` の両方に対応します。
  - `Browse` ボタンでファイルを選びます。

- `R2 FASTQ`
  - R2 側の FASTQ ファイルを指定します。
  - `.fastq` と `.fastq.gz` の両方に対応します。
  - `Browse` ボタンでファイルを選びます。

- `Output TSV`
  - IgBLAST の生AIRR TSV出力先です。
  - R1/R2 FASTQ を選ぶと、標準では作業フォルダ内の `Results of RG Paired Fastq IgBLAST AIRR tsv` に自動設定されます。
  - ファイル名にはFASTQの共通サンプル名を使います。例: `KKF103hG_S57_L001_R1_001.fastq` と `KKF103hG_S57_L001_R2_001.fastq` から `KKF103hG_S57_L001.airr.tsv` を作ります。
  - 実行後、同じフォルダに `KKF103hG_S57_L001.R1.airr.tsv`、`KKF103hG_S57_L001.R2.airr.tsv`、`KKF103hG_S57_L001.integrated.tsv`、`KKF103hG_S57_L001.integrated_counts.tsv`、`KKF103hG_S57_L001.integrated_counts.xlsx` も自動作成します。
  - `Browse` ボタンで保存先を指定します。

- `Keep query FASTA`
  - IgBLAST に渡す中間 FASTA の保存先です。
  - R1/R2 FASTQ を選ぶと、標準では `Results of RG Paired Fastq IgBLAST AIRR tsv` 内の `<サンプル名>.queries.fasta` に自動設定されます。
  - 空欄にした場合は、一時ファイルとして作成し、解析後に削除します。
  - トラブルシューティングやquery内容の確認をしたい場合は指定してください。

### IgBLAST 設定

- `igblastn`
  - `igblastn` の実行ファイルを指定します。
  - PATH が通っている場合は `igblastn` のままで使えます。
  - PATH が通っていない場合は、`igblastn.exe` を `Browse` で指定します。

- `V DB prefix`
  - IgBLAST の V germline database prefix を指定します。
  - 標準例: `...\IgBlast_refdata_edit_imgt\db\IMGT_IGHV.imgt`
  - `.nhr`, `.nin`, `.nsq` などの拡張子を除いた名前を指定します。

- `D DB prefix`
  - IgBLAST の D germline database prefix を指定します。
  - heavy chain 解析では通常指定します。
  - 標準例: `...\IgBlast_refdata_edit_imgt\db\IMGT_IGHD.imgt`
  - light chain だけを対象にする場合など、不要な条件では空欄でも使えます。

- `J DB prefix`
  - IgBLAST の J germline database prefix を指定します。
  - 標準例: `...\IgBlast_refdata_edit_imgt\db\IMGT_IGHJ.imgt`

- `Aux file`
  - IgBLAST の auxiliary data file を指定します。
  - 標準例: `...\IgBlast_refdata_edit_imgt\optional_file\human_gl.aux`
  - CDR/FWR 領域などの情報を出したい場合に重要です。

- `Organism`
  - IgBLAST の `-organism` に渡す値です。
  - human サンプルなら通常 `human` です。

- `Domain system`
  - IgBLAST の `-domain_system` に渡す値です。
  - 通常は `imgt` を使います。

- `Seq type`
  - IgBLAST の `-ig_seqtype` に渡す値です。
  - 免疫グロブリンなら通常 `Ig` です。

- `Threads`
  - IgBLAST の実行スレッド数です。
  - GUIの標準値は `4` です。
  - 大きな FASTQ ではCPUに余裕があれば `8` などに増やすと速くなる場合があります。

### R1/R2 の処理条件

- `Reads`
  - `both`: R1 と R2 の両方を IgBLAST に渡します。
  - `r1`: R1 だけを使います。
  - `r2`: R2 だけを使います。
  - デフォルトは `both` です。

- `R1 orientation`
  - R1 の向きを指定します。
  - デフォルトは `forward` です。

- `R2 orientation`
  - R2 の向きを指定します。
  - デフォルトは `reverse-complement` です。
  - paired-end sequencing では R2 が逆向きに読まれていることが多いため、この設定を標準にしています。

- `Trim left R1`
  - R1 の左端から削る塩基数です。
  - primer や低品質領域を除く場合に使います。

- `Trim right R1`
  - R1 の右端から削る塩基数です。

- `Trim left R2`
  - R2 の左端から削る塩基数です。

- `Trim right R2`
  - R2 の右端から削る塩基数です。

- `Min length`
  - IgBLAST に渡す最小配列長です。
  - trim 後にこの長さ未満になった read は除外します。

- `Max N rate`
  - 配列中の `N` の割合の上限です。
  - 例: `0.10` なら、N が10%を超える read を除外します。

- `Query name`
  - FASTA query 名の形式です。
  - デフォルトは `{read_id}|{read}` です。
  - 例: `M03603:...|R1`, `M03603:...|R2`

- `Require matching R1/R2 IDs`
  - チェックあり: R1 と R2 の read ID が対応していない場合に停止します。
  - チェックなし: read ID が違っても処理を続けます。
  - 通常はチェックありを推奨します。

### 実行ボタン

- `Browse`
  - ファイルの場所や保存先を選ぶためのボタンです。

- `Run`
  - FASTQ から中間 FASTA を作成し、IgBLAST を実行して AIRR TSV を作成します。
  - IgBLAST完了後、生AIRR TSVからR1 TSV、R2 TSV、統合TSV、集計TSV、集計Excelも作成します。
  - 実行中は処理が完了するまで待ちます。
  - 成功すると、出力TSV、R1 TSV、R2 TSV、統合TSV、集計TSV、集計Excelの場所と処理したread数が表示されます。
  - GUIではIgBLASTの作業ファイルをPC内のローカル一時フォルダで作り、完了後にResultsフォルダへコピーします。Desktop/OneDrive配下へ巨大TSVを長時間直接書き続けることを避けるためです。

### ログ表示欄

GUI 下部のログ欄には、処理開始、完了、エラー内容が表示されます。エラーが出た場合は、まず FASTQ のパス、IgBLAST のパス、database prefix、aux file の指定を確認してください。

## CLI

R1/R2 から AIRR TSV を作る例:

```powershell
rg-paired-fastq-igblast-airr-tsv run `
  --r1 "C:\Users\Yohei Funakoshi\Desktop\RG Paired Fastq IgBLAST AIRR tsv\sample_R1.fastq.gz" `
  --r2 "C:\Users\Yohei Funakoshi\Desktop\RG Paired Fastq IgBLAST AIRR tsv\sample_R2.fastq.gz" `
  --out "C:\Users\Yohei Funakoshi\Desktop\RG Paired Fastq IgBLAST AIRR tsv\Results of RG Paired Fastq IgBLAST AIRR tsv\sample.airr.tsv" `
  --germline-db-v "C:\Users\Yohei Funakoshi\Desktop\RG Paired Fastq IgBLAST AIRR tsv\refdata\IgBlast_refdata_edit_imgt\db\IMGT_IGHV.imgt" `
  --germline-db-d "C:\Users\Yohei Funakoshi\Desktop\RG Paired Fastq IgBLAST AIRR tsv\refdata\IgBlast_refdata_edit_imgt\db\IMGT_IGHD.imgt" `
  --germline-db-j "C:\Users\Yohei Funakoshi\Desktop\RG Paired Fastq IgBLAST AIRR tsv\refdata\IgBlast_refdata_edit_imgt\db\IMGT_IGHJ.imgt" `
  --auxiliary-data "C:\Users\Yohei Funakoshi\Desktop\RG Paired Fastq IgBLAST AIRR tsv\refdata\IgBlast_refdata_edit_imgt\optional_file\human_gl.aux" `
  --organism human `
  --domain-system imgt `
  --num-threads 4 `
  --work-dir "$env:LOCALAPPDATA\PairedFastqIgblastAirrTsv\work"
```

IgBLAST に渡す FASTA だけを作る例:

```powershell
rg-paired-fastq-igblast-airr-tsv prepare `
  --r1 "C:\Users\Yohei Funakoshi\Desktop\RG Paired Fastq IgBLAST AIRR tsv\sample_R1.fastq.gz" `
  --r2 "C:\Users\Yohei Funakoshi\Desktop\RG Paired Fastq IgBLAST AIRR tsv\sample_R2.fastq.gz" `
  --out-fasta "C:\Users\Yohei Funakoshi\Desktop\RG Paired Fastq IgBLAST AIRR tsv\Results of RG Paired Fastq IgBLAST AIRR tsv\sample.queries.fasta"
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

## 出力されるファイル

- `*.airr.tsv`
  - IgBLAST `-outfmt 19` による AIRR TSV です。

- `*.queries.fasta`
  - `--query-fasta` または GUI の `Keep query FASTA` を指定した場合だけ保存されます。
  - IgBLAST に渡した query 配列の確認に使えます。

## トラブルシューティング

- `igblastn` が見つからない
  - GUI の `igblastn` 欄に `igblastn.exe` のフルパスを指定してください。

- database が見つからない
  - `V DB prefix`, `D DB prefix`, `J DB prefix` が正しいか確認してください。
  - `.nhr`, `.nin`, `.nsq` などの拡張子を付けず、prefix だけを指定します。

- AIRR TSV が空に近い
  - `R1 orientation`, `R2 orientation` を確認してください。
  - `Min length` や `Max N rate` で read が除外されすぎていないか確認してください。
  - germline DB と organism がサンプルに合っているか確認してください。

- R1/R2 ID mismatch と出る
  - R1 と R2 が同じサンプル由来か確認してください。
  - read ID の形式が特殊な場合は `Require matching R1/R2 IDs` を外すか、CLI では `--allow-id-mismatch` を使います。

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

# RG Paired Fastq IgBLAST AIRR tsv 議論サマリー

作成日: 2026-05-03  
対象: RG社PCRシーケンス由来BCRレパトア解析データ用GUI / CLI

## 最終注釈

2026-05-03時点の最終判断として、RG社Excelとの比較では、総リード数が一致し、Vコール、Jコール、`junctionAA` / `final_junction_aa` で定義したユニークリード数も概ね近い結果でした。そのため、この版では「アリルなしV候補セット + アリルなしJ候補セット + `final_junction_aa`」をユニークリードの定義とし、これ以上複雑化せずいったん完成版とします。

RG社Excelとの差分および完成判断は [`RG_company_excel_comparison_note_20260503.md`](RG_company_excel_comparison_note_20260503.md) に追記しました。

## 目的

この文書は、RG社データ用の `RG Paired Fastq IgBLAST AIRR tsv` について、これまでの議論で決めた解析方針、採用した理由、将来の検討点をまとめるものです。

共同研究者がGUIを使って解析したり、今後このリポジトリを修正したりする際に、どのような考え方で現在の仕様になったかを追跡できるようにすることを目的とします。

## 解析対象と前提

- 対象はRG社で実施されたPCRシーケンス由来のBCRレパトア解析データです。
- RG社データはUMIなしデータとして扱います。
- R1はC領域側、R2はV領域側から読むpaired-end readとして扱います。
- FASTQは標準設定ではトリムしません。
- R1とR2はマージしません。
- 大きなFASTAを一度にIgBLASTへ渡すと停止したように見えることがあるため、標準では `IgBLAST batch size = 10000` として分割実行します。
- IgBLASTの参照データはIMGT由来のIGHV/IGHD/IGHJをIgBLAST用に整形したBLAST DBを使う想定です。
- IgBLASTの確認済みバージョンは `igblastn 1.21.0` です。

## 基本コンセプト

現在のRG版は「正確性を厳しく追い込む版」ではなく、「まず情報を落とさない暫定版」です。

R1/R2のどちらかに有用な情報がある場合、それを簡単に捨てないことを重視しています。一方で、共同研究者がすぐに扱える集計表も必要なため、追跡用データと集計用データを分けています。

## 出力ファイル

GUIでRunすると、Resultsフォルダに次のファイルを作成します。

| ファイル | 位置づけ |
|---|---|
| `<sample>.airr.tsv` | IgBLAST `-outfmt 19` の生AIRR TSV。R1/R2両方を含む |
| `<sample>.R1.airr.tsv` | R1だけを抜き出したAIRR TSV |
| `<sample>.R2.airr.tsv` | R2だけを抜き出したAIRR TSV |
| `<sample>.integrated.tsv` | read pair単位の追跡用サマリー。全pairを残す |
| `<sample>.integrated_counts.tsv` | 集計対象pairだけをユニークリード単位にまとめたTSV |
| `<sample>.integrated_counts.xlsx` | `integrated_counts.tsv` と同じ内容のExcelファイル |

## R1/R2をマージしない理由

一般的にはR1/R2をマージしてからIgBLAST解析する方法があります。しかし、マージ条件によってreadが落ちること、マージできないreadにも有用なV/J/CDR3情報が残っている可能性があることを考慮し、RG版ではR1/R2をマージしない方針にしました。

この方針では、R1とR2は独立したIgBLAST queryとして解析されます。query名には `read_id|R1`、`read_id|R2` のようにread側を付けます。

## 統合TSVの役割

`integrated.tsv` はAIRR標準そのものではなく、共同解析で確認しやすくするための追跡用サマリーです。

このファイルには、R1/R2それぞれの値、暫定的なfinal値、どちらのreadを採用したか、採用理由、集計に入れるかどうかを残します。

重要なのは、`integrated.tsv` では全pairを残すことです。生AIRR TSVであるR1/R2 TSVと合わせて確認すれば、後から別の統合ルールで再解析できます。

## final_junction_aa の暫定ルール

`junction_aa` は最終的なユニークリード判定で重要なため、次のシンプルなルールで `final_junction_aa` を決めます。

| 状況 | 採用ルール |
|---|---|
| R1/R2で一致 | その値を採用 |
| 片方だけに値がある | 値がある方を採用 |
| R1/R2で不一致、長さが異なる | 長い方を採用 |
| R1/R2で不一致、長さが同じ | R2を採用 |

不一致を理由にpairを落とすことはしません。これは「落とさない」暫定版というコンセプトに基づきます。

## V/J/D/productive/junction の暫定ルール

- `v_call` はR2を優先します。R2はV領域側から読むためです。
- `j_call`、`d_call`、`productive`、`junction` は、原則として `final_junction_aa` に採用したread側を優先します。
- 採用したread側が空欄の場合は、もう一方のread側の値を使います。
- R1/R2の元値は `r1_*`、`r2_*` として残します。

## 集計対象にする条件

`integrated.tsv` では全pairを残しますが、`integrated_counts.tsv` / `.xlsx` に入れるのは `include_in_counts=true` のpairだけです。

`include_in_counts=true` になる条件は次の通りです。

- IgBLASTの `productive` がtrue相当であること
- V候補セットがあること
- J候補セットがあること
- `final_junction_aa` があること
- `final_junction_aa` にstop `*` がないこと
- `final_junction_aa` が `C` で始まること
- `final_junction_aa` が `W` または `F` で終わること
- `final_junction_aa` の長さが5-40 amino acidsであること

除外されたpairは `integrated.tsv` の `exclude_reason` に理由を残します。

## ユニークリードの定義

現在のRG版では、ユニークリードを次の3要素で定義します。

```text
unique_v_gene_set + unique_j_gene_set + final_junction_aa
```

### V/J候補セット

IgBLASTではVやJに複数候補が出ることがあります。現在のルールでは、候補を1つに無理に選ばず、候補セットとして扱います。

例:

```text
IGHV4-61*01,IGHV4-59*01
```

は、アリルを外してソートし、次のように扱います。

```text
IGHV4-59,IGHV4-61
```

### exact candidate set方式

候補セットはexact candidate set方式で扱います。つまり、候補セットが完全に同じ場合だけ同一とします。

例:

| V候補 | J候補 | junctionAA | 判定 |
|---|---|---|---|
| `IGHV4-59,IGHV4-61` | `IGHJ4` | `CVQGFDYW` | 同じ候補セットなら同一 |
| `IGHV4-61` | `IGHJ4` | `CVQGFDYW` | 上とは別 |

この方式を採用した理由は、RG社Excel内に、同じJ、同じCDR3でもV候補セットの違いで別行になっている例が確認されたためです。

### アリルを使わない理由

RG社や他社解析の納品形式では、アリルまで使わずgene levelで集計されていることが多いと判断しました。そのため、集計ではアリルを外します。

ただし、`integrated.tsv` にはIgBLASTのアリル付きcallを残します。必要であれば後からアリル込みの解析を再構築できます。

### D callをユニークリード判定に使わない理由

D遺伝子はコールされないことが多く、短い領域でもあるため、ユニークリード判定に入れると過剰に分割される可能性があります。

そのため、D callは `integrated.tsv` には残しますが、`integrated_counts.tsv` / `.xlsx` のユニークリード判定には使いません。

例:

| V | D | J | junctionAA | counts上の扱い |
|---|---|---|---|---|
| 同じ | Dあり | 同じ | 同じ | 同じユニークリード |
| 同じ | Dなし | 同じ | 同じ | 同じユニークリード |

この場合、`integrated_counts.tsv` では1行にまとまり、`read_pair_count=2` になります。Dの有無を確認したい場合は `integrated.tsv` を見ます。

## productive と canonical junctionAA

RG社Excelの `In frame` とIgBLASTの `productive` は完全に同じ判定ではありません。ただし、RG社納品Excelに近い見方をするため、集計表では `productive=true` かつcanonicalなjunctionAAを満たすpairを主に扱います。

ここでのcanonical junctionAAは、実用上のCDR3らしさを示す簡易条件です。

- Cで始まる
- WまたはFで終わる
- stop `*` を含まない
- 長さが5-40 amino acids

## RG社Excelとの比較から得た示唆

KKF103データを用いて、RG社作成Excelと本ツールの出力を比較しました。完全一致を目的にしたわけではありませんが、ユニークリード数などの大枠が近いことを確認しました。

一方で、RG社はR1/R2から得られる情報を最終的に1つの納品Excelに反映していると考えられます。本ツールでも、R1/R2の生AIRR TSVを残しつつ、暫定的な統合countsを作る設計にしました。

## 現在の設計の利点

- R1/R2の生AIRR TSVが残るため、世界中の研究者が標準フォーマットで再確認できます。
- `integrated.tsv` に判断過程が残るため、後からルールを変更して再解析できます。
- `integrated_counts.tsv` / `.xlsx` は1ユニークリード1行で、共同研究者が解析しやすい形式です。
- 正確性重視版を将来作る場合にも、現行版の出力を基礎データとして使えます。

## 現時点の注意点

- `integrated.tsv` はAIRR標準そのものではありません。追跡用の独自サマリーです。
- `integrated_counts.tsv` / `.xlsx` は暫定ルールによる集計結果です。
- conflictを厳密に扱う正確性重視版ではありません。
- IgBLASTは標準でbatch実行します。これは解析結果のルールを変えるものではなく、大きなFASTAを安定して処理するための実行方法です。
- `productive_true_count` と `canonical_junction_aa_count` は、countsに入る時点でフィルタ済みのため、多くの場合 `read_pair_count` と同じ値になります。
- D callはcountsのキーには含めません。Dを使った解析が必要な場合は `integrated.tsv` から再解析します。

## 将来の検討点

- conflictをより厳しく扱う正確性重視版を別に作るか。
- R1/R2でV/J/junctionAAが不一致の場合の品質指標を導入するか。
- read qualityやalignment scoreを使ったfinal値選択を行うか。
- RG社Excelの列構成に近いsummary出力を追加するか。
- CPM社用のUMI collapse版と、どこまで共通化するか。

## 現在採用している立場

現在のRG版は、共同研究者が使いやすいこと、後から追跡できること、情報を過度に落とさないことを優先した暫定解析システムです。

生データに近いR1/R2 AIRR TSVを残し、その上でシンプルな統合countsを作る、という二層構造を採用しています。

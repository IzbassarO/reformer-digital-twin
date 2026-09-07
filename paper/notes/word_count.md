# Word count of the manuscript

Measured on the current `paper/main.tex` and `paper/supplementary.tex` with TeXcount 3.1.1
(2018 Oct 28), the copy shipped with TinyTeX. Nothing in the manuscript was changed to obtain these
numbers.

## Correction to the premise

**texcount was not counting `thebibliography` as body text.** Its default already excludes the
environment, so the body figures reported in earlier commits (7775 for `main.tex`) were already
free of the 42 embedded references. Three independent checks agree:

```
$ texcount -total main.tex                     # default
Words in text: 7775

$ texcount -total -incbib main.tex             # force the bibliography to be counted
Words in text: 8843

$ texcount -total <main.tex with the thebibliography block physically deleted>
Words in text: 7775
```

The default run and the physically-stripped file give the identical 7775, and
7775 + 1068 = 8843 closes against the `-incbib` run. The 1068-word difference is the bibliography.

`-incbib` changes only the body count; captions and table text are unaffected (673 either way).

## Invocations

Both documents are single files with no `\input`, so `-inc` is a no-op here; it is kept below only
because earlier measurements used it, and it gives identical output.

| quantity | invocation (run from `paper/`) | field read |
|---|---|---|
| body, bibliography excluded | `texcount -inc -total <doc>.tex` | `Words in text` |
| body, bibliography included | `texcount -inc -total -incbib <doc>.tex` | `Words in text` |
| headers | `texcount -inc -total <doc>.tex` | `Words in headers` |
| captions + table text together | `texcount -inc -total <doc>.tex` | `Words outside text (captions, etc.)` |

texcount has no option that separates figure captions from table text — it reports both in the one
"outside text" bucket. The split below is obtained by extracting the environments into two
throwaway files and counting each:

```
# every \begin{figure}...\end{figure} of <doc>.tex, wrapped in that document's own preamble
$ texcount -total /tmp/<doc>_figure.tex        # -> figure captions
# every \begin{table}...\end{table}, same wrapper
$ texcount -total /tmp/<doc>_table.tex         # -> table text (caption + cell contents)
```

The bibliography alone was measured the same way, by extracting the block into a wrapper file and
counting it with `-incbib`:

```
$ texcount -total -incbib /tmp/<doc>_bibonly.tex
```

Both splits are checked against the whole-document figure: captions + tables must equal the
"outside text" total, and body + bibliography must equal the `-incbib` body.

## `paper/main.tex`

| quantity | words |
|---|---|
| body, excluding `thebibliography` | **7775** |
| `thebibliography` alone (42 entries) | **1068** |
| headers | **128** |
| figure captions (10 figures) | **459** |
| table text (5 tables, caption + cells) | **214** |
| *captions + tables, cross-check against "outside text"* | *673 = 673* |
| *body + bibliography, cross-check against `-incbib`* | *8843 = 8843* |
| **journal count, "words including tables"** | **8448** |

where the journal count is body excluding references + captions + table text
= 7775 + 459 + 214.

Headers are excluded from that total on the reading that section titles are not "words including
tables"; adding them would give 8576.

## `paper/supplementary.tex`

| quantity | words |
|---|---|
| body, excluding `thebibliography` | **1584** |
| `thebibliography` alone (7 entries) | **192** |
| headers | **93** |
| figure captions (7 figures) | **281** |
| table text (9 tables, caption + cells) | **277** |
| *captions + tables, cross-check* | *558 = 558* |
| *body + bibliography, cross-check* | *1776 = 1776* |
| **journal count, "words including tables"** | **2142** |

## Two things worth knowing about what these numbers do and do not contain

**Table cells count as very few words.** texcount scores numeric cell contents as symbols, not
words, so a table's contribution is dominated by its caption. `main.tex` has 5 tables carrying 214
words in total, of which the great majority is caption prose. This is why moving a seven-column
table of Sobol indices out of the main text freed only about 20 words of table text.

**Equations are excluded.** Displayed and inline maths are reported separately (172 inline, 5
displayed in `main.tex`) and appear in none of the counts above. A journal that counts a displayed
equation as a fixed number of words would add to these figures.

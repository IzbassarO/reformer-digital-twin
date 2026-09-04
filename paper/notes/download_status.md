# Download status (2026-09-04)

All files are stored in `literature/`, which is git-ignored. PDFs must never be committed.

## Succeeded

| File | Source | Method | Result |
|---|---|---|---|
| `literature/latham2008_msce_thesis.pdf` | Latham (2008), M.Sc.E. thesis, Queen's University. Handle `1974/1650`, item page <https://queensu.scholaris.ca/items/5ca025eb-0936-46bf-a257-3af4c3ca6d95> | `curl` on the Scholaris bitstream URL <https://queensu.scholaris.ca/server/api/core/bitstreams/d703f47a-c138-4699-8e88-1a11e0773040/content> | HTTP 200, `application/pdf`, 2.2 MB, 279 pages |
| `literature/yeh2021_applsci_11_231.pdf` | Yeh (2021), *Appl. Sci.* 11(1):231, doi:10.3390/app11010231 | `curl` on <https://mdpi-res.com/d_attachment/applsci/applsci-11-00231/article_deploy/applsci-11-00231.pdf> (MDPI static asset host) | HTTP 200, `application/pdf`, 5.6 MB, 20 pages |
| `literature/nabavi2025_arxiv_2507.07641.pdf` | Nabavi, Guo & Wang (2025), arXiv:2507.07641 | `curl` on <https://arxiv.org/pdf/2507.07641> | HTTP 200, `application/pdf`, 1.5 MB, 36 pages |

## Failed

| Target | URL(s) tried | Result | Action |
|---|---|---|---|
| Yeh (2021) via the MDPI article page | <https://www.mdpi.com/2076-3417/11/1/231/pdf> | HTTP 403 (`text/html`, bot block) | Not needed: the same PDF was obtained from `mdpi-res.com` (see above). |
| Haidemenopoulos et al. (2021), *Eur. J. Mater.* 1(1):1-22, doi:10.1080/26889277.2021.1994841 (open access) | `https://www.tandfonline.com/doi/pdf/10.1080/26889277.2021.1994841`, `.../doi/pdf/...?needAccess=true`, `.../doi/epdf/...`, `.../doi/full/...` | HTTP 403 on all four; response body is a Cloudflare "Just a moment..." JavaScript challenge (captcha-style bot check). Unpaywall reports the only OA location as the publisher DOI itself, no repository copy. | **Manual download required**: open <https://doi.org/10.1080/26889277.2021.1994841> in a browser and save the PDF as `literature/haidemenopoulos2021_eurjmater.pdf`. No library access needed (article is open access). |

## Bibliography fetch notes (CrossRef content negotiation)

- 16 of 20 supplied identifiers resolved directly.
- `10.1016/j.msea.2013.05.052` resolves to an unrelated paper; the intended Whittaker, Wilshire & Brear (2013) paper was recovered from PII `S0921509313005716` as **10.1016/j.msea.2013.05.041**.
- `10.1016/j.ijpvp.2023.104909` resolves to an unrelated paper; the intended Guguloth & Roy (2023) paper was recovered from PII `S0308016123000558` as **10.1016/j.ijpvp.2023.104938**.
- `10.1016/j.compchemeng.2017.04.026` does not exist ("DOI Not Found"). Tran et al. (2017) was located by a CrossRef bibliographic query as **10.1016/j.compchemeng.2017.04.013**.
- `10.11503/nims.1020` and `10.11503/nims.1042` are Japan Link Center DOIs; doi.org returned HTTP 406 for BibTeX. Metadata was taken from the JaLC API and entered by hand as `@techreport` entries.
- arXiv BibTeX was fetched from <https://arxiv.org/bibtex/2507.07641>.

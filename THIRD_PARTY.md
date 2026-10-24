# Third-party components

CRISPRme+ distributes (in its Docker image and/or scorer environments) the
following third-party components. Their copyright notices and permission text
are reproduced below as required by their licenses. This is a living list;
components added later are appended here.

---

## CRISPR-Bulge

- **Used for:** the deep-learning off-target activity score (the CRISPR-Bulge
  scorer environment).
- **Source:** https://github.com/OrensteinLab/CRISPR-Bulge
  (reference: Yaish & Orenstein, *Nucleic Acids Res.* 2024, doi:10.1093/nar/gkae428)
- **License:** MIT

```
MIT License

Copyright (c) 2024 Orenstein Lab

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

---

## impg

- **Used for:** assembly-search haplotype-vs-haplotype coordinate lookup.
- **Source:** https://github.com/pangenome/impg
- **License:** MIT

```
MIT License

Copyright (c) 2024 Erik Garrison

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

---

## CFD score tables (reference data)

- **Used for:** the CFD (Cutting Frequency Determination) off-target score, computed
  by CRISPRme+'s own implementation (`PostProcess/cfd_score.py`).
- **Source:** the per-position mismatch/gap and PAM penalty values published by
  Doench et al., *Nat. Biotechnol.* 2016 (doi:10.1038/nbt.3437; supplementary data and
  Fig. 5c,d). These are published scientific measurements (data), re-serialized here as
  reviewable JSON (`PostProcess/cfd_mismatch_scores.json`, `cfd_pam_scores.json`; see the
  committed generator `scripts/generate_cfd_tables.py`).

---

## pyliftover

- **Used for:** assembly-search coordinate lifting (lifting off-target
  coordinates through chain files).
- **Source:** https://github.com/konstantint/pyliftover
- **License:** MIT

```
Copyright (c) 2013, Konstantin Tretyakov

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

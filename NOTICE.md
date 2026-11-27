# Third-party notices

This project (`hrdness`) is licensed under the PolyForm Noncommercial License 1.0.0 (see `LICENSE`).
It incorporates or depends on the following third-party material.

## 1. ssGSEA (Broad Institute) — `hrdness/ssgsea.py`

`hrdness/ssgsea.py` is a Python port of the single-sample GSEA projection code (formerly "ssGSEAProjection", R; authors of the
original method and code: Pablo Tamayo, Chet Birger et al.) distributed as the ssGSEA module for GenePattern
(https://github.com/GSEA-MSigDB/ssGSEA-gpmodule), which is made available under the following BSD-style license.
Method reference: Barbie DA et al., Nature 2009;462:108-112.

The ssGSEA module for GenePattern is distributed under the following BSD-style license:

Copyright © 2003-2020 The Broad Institute, Inc., Massachusetts Institute of Technology and Regents of 
the University of California.  All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided 
that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the 
following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the 
following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the names of the copyright holders nor the names of their contributors may be used to endorse or 
promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED “AS IS.”  THE COPYRIGHT HOLDERS MAKE NO EXPRESS OR 
IMPLIED REPRESENTATIONS OR WARRANTIES OF ANY KIND REGARDING THE SOFTWARE AND 
COPYRIGHT, INCLUDING, BUT NOT LIMITED TO, WARRANTIES OF MERCHANTABILITY, FITNESS 
FOR A PARTICULAR PURPOSE, CONFORMITY WITH ANY DOCUMENTATION, 
NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT 
DISCOVERABLE. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE 
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
DAMAGES (INCLUDING, BUT NOT LIMITED TO PROCUREMENT OF SUBSTITUTE GOODS OR 
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR 
TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF 
THIS SOFTWARE, EVEN IF ADVISED OF, HAVE REASON TO KNOW, OR IN FACT SHALL KNOW OF 
THE POSSIBILITY OF SUCH DAMAGE.

If, by operation of law or otherwise, any of the aforementioned warranty disclaimers are determined 
inapplicable, your sole remedy, regardless of the form of action, including, but not limited to, negligence 
and strict liability, shall be replacement of the software with an updated version if one exists.

Development of the ssGSEA module for GenePattern was funded, in whole or in part, by the National Institutes of Health.

The names of the copyright holders and their contributors are not used to endorse or promote `hrdness`.

## 2. TCGA data — `data/reference_TCGA/`

The reference expression tables are derived from open-access data of The Cancer Genome Atlas (TCGA) Research Network
(https://www.cancer.gov/tcga), obtained via the NCI Genomic Data Commons (GDC; RNA-seq) and UCSC Xena (Affymetrix HT_HG-U133A).
Use is subject to the NIH Genomic Data Sharing policy for unrestricted-access data: do not attempt to identify individual
participants, and acknowledge the data source and repositories in publications. HRD labels (`TCGA_OV_HRD_label_rnaseq_case.tsv`)
were derived as described in Takamatsu et al., Br J Cancer 2023;128:1095-1104 (doi:10.1038/s41416-022-02122-9).

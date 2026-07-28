.. _thresholding:

=========================================================
Mass-Tag Detection: Background Model and Thresholding
=========================================================

This document specifies the detection-limit model used to call photocleavable
mass-tag (PCMT) positivity in single-cell MALDI spectra, and the cell-type
assignment that follows from it.

It documents **the code as executed**: ``script/sc_dictionary.py`` and
``script/sc_dictionary.ipynb``. Every formula below is transcribed from those
files. Functions present in the module but not called by the notebook are listed
in :ref:`not-used` and are **not** part of this pipeline.

.. contents:: Contents
   :local:
   :depth: 2


Notation
========

.. list-table::
   :header-rows: 1
   :widths: 12 88

   * - Symbol
     - Meaning
   * - :math:`i`
     - cell index (one MS1 spectrum), :math:`i = 1 \dots N`
   * - :math:`j`
     - mass-tag index, :math:`j = 1 \dots M` (:math:`M = 5`)
   * - :math:`b`
     - acquisition file (batch) index; :math:`b(i)` is the file cell *i* came from
   * - :math:`x_{ij}`
     - extracted intensity of tag *j* in cell *i*
   * - :math:`t_j`
     - cell type (lineage) that tag *j* marks
   * - :math:`k`
     - threshold multiplier, ``K = 2.0``
   * - :math:`\tau_{bj}`
     - detection threshold for tag *j* in file *b*

Tags and lineages, as configured in the notebook (Cell 2):

.. list-table::
   :header-rows: 1
   :widths: 14 20 30

   * - Tag
     - m/z (Table S2)
     - Lineage :math:`t_j`
   * - NeuN
     - 1308.69
     - neuron
   * - GFAP
     - 1011.53
     - astrocyte
   * - MBP
     - 1365.71
     - oligodendrocyte
   * - NFL
     - 1345.72
     - neuron
   * - SYN-I
     - 1482.75
     - neuron


Pipeline order
==============

The order is load-bearing: quality filtering happens **before** any statistic is
computed, because every estimator below operates on whatever matrix it is handed
and cannot distinguish a blank acquisition from a real cell.

.. code-block:: text

   1. load spectra, attach spot positions
   2. QC: drop re-ablated block + near-empty spectra      <-- BEFORE thresholding
   3. extract tag intensities  x_ij
   4. seed pass (provisional positivity)
   5. negative set per tag
   6. background centre + negative-set detection rate
   7. mode decision (background / floor), GLOBAL per tag
   8. sigma per tag per file
   9. threshold, detection, classification


Step 1 — Intensity extraction
=============================

``marker_intensity_from_spectra(spectra, target_mz, tol=0.01)``

For each cell and tag, the intensity is the **maximum** peak intensity inside a
fixed window around the tag's m/z:

.. math::

   x_{ij} \;=\;
   \begin{cases}
     \max \{\, I_{ip} \;:\; |m_{ip} - \mu_j| \le \delta \,\} & \text{if any peak lies in the window} \\
     0 & \text{otherwise}
   \end{cases}

where :math:`m_{ip}, I_{ip}` are the m/z and intensity of peak *p* in cell *i*,
:math:`\mu_j` is the tag's monoisotopic m/z and :math:`\delta =` ``TOL`` :math:`= 0.01`
Da.

.. note::
   :math:`x_{ij} = 0` means *no peak was detected in the window*. Because the
   spectra are centroided, a zero is an absence of a peak, not a low intensity.
   This distinction drives the ``nonzero`` handling throughout.


Step 2 — Quality filtering (before any statistic)
=================================================

Two masks are combined (notebook Cell 4). A cell is kept only if it passes both.

**(a) Re-ablated block.** The file matching ``PARTII_TAG`` re-visited positions
already measured in an earlier run, so its leading acquisitions carry depleted
material. They are removed by acquisition order:

.. math::

   \text{drop } i \quad\text{if}\quad b(i) = b_{\text{PartII}} \;\wedge\; o_i < n_{\text{drop}}

with ``N_DROP = 299`` and :math:`o_i` the acquisition index within the file.

**(b) Near-empty spectra.** ``flag_blank_spectra(npeaks, batch, rel=0.10)`` keeps
cell *i* only if its peak count exceeds a fraction of its **own file's** median:

.. math::

   \text{keep } i \quad\Longleftrightarrow\quad
   n^{\text{peaks}}_i \;>\; \max\!\left( r \cdot \operatorname{median}_{\,i' : b(i')=b(i)} n^{\text{peaks}}_{i'},\; a_{\min} \right)

with :math:`r =` ``rel`` :math:`= 0.10` and :math:`a_{\min} =` ``abs_min`` :math:`= 0`.

The cutoff is relative because spectral richness differs by more than an order of
magnitude between acquisitions; a single absolute peak count would mean different
things in different files.


Step 3 — Seed pass
==================

The background must be measured in cells *lacking* a tag, which requires a
provisional notion of positivity. ``_negative_centres`` builds it when no seed is
supplied, using ``_one_threshold(..., method="median_mad")``:

.. math::

   \tilde{\tau}_j \;=\; \operatorname{median}\{\, x_{ij} : x_{ij} > 0 \,\} \;+\; k\,\sigma_j

.. math::

   s_{ij} \;=\; \mathbb{1}\!\left[\, x_{ij} > \tilde{\tau}_j \,\right]

where :math:`\sigma_j` is the non-zero MAD sigma of :ref:`Step 6 <step-sigma>`,
computed on the same column.

The seed is deliberately crude. It decides only **which cells** are treated as
background; the threshold *value* comes from their intensities.


Step 4 — Negative set
=====================

``_negative_mask(..., negative_set="marker")``, the setting used
(``NEG_SET = "marker"``, Cell 5):

.. math::

   \mathcal{N}_j \;=\; \{\, i \;:\; s_{ij} = 0 \,\}

i.e. every cell the seed did **not** call positive for tag *j* itself.

The module also implements ``negative_set="lineage"``, which excludes all cells
positive for any tag of the same lineage (and, if
``require_other_lineage=True``, additionally requires positivity for a different
lineage):

.. math::

   \mathcal{N}^{\text{lin}}_j \;=\;
   \Big\{\, i \;:\; \textstyle\bigvee_{j' : t_{j'} = t_j} s_{ij'} = 0 \,\Big\}
   \;\cap\;
   \Big\{\, i \;:\; \textstyle\bigvee_{j' : t_{j'} \ne t_j} s_{ij'} = 1 \,\Big\}

This was **not** the setting used in the documented run.


Step 5 — Background centre and detection rate
=============================================

``_negative_centres`` computes two quantities per tag, both over the negative set,
using **non-zero values only** (``nonzero=True``):

**Negative-set detection rate** — the fraction of background cells that still show
a peak at the tag's m/z:

.. math::

   d_j \;=\; \frac{1}{|\mathcal{N}_j|}\sum_{i \in \mathcal{N}_j} \mathbb{1}\!\left[\, x_{ij} > 0 \,\right]

**Background centre** — the typical intensity of a background peak, conditional on
one being detected:

.. math::

   c_j \;=\;
   \begin{cases}
     \operatorname{median}\{\, x_{ij} : i \in \mathcal{N}_j,\; x_{ij} > 0 \,\}
       & \text{if } |\{i \in \mathcal{N}_j : x_{ij} > 0\}| \ge n_{\min} \\[4pt]
     0 & \text{otherwise}
   \end{cases}

with :math:`n_{\min} =` ``min_cells`` :math:`= 20`.

.. warning::
   :math:`c_j` is **conditional on detection**: cells in the negative set with no
   peak contribute nothing rather than contributing a zero. Including the zeros
   would place the median on the zero pile whenever more than half the background
   cells lack a peak, collapsing both :math:`c_j` and :math:`\sigma_j` to 0 and
   admitting every cell as positive.


Step 6 — Mode: background or floor
==================================

The threshold takes one of two forms depending on whether a background exists to
be measured at all. The decision uses :math:`d_j` computed on the **full matrix**
(all files pooled), so it is made **once per tag**:

.. math::

   \text{mode}_j \;=\;
   \begin{cases}
     \texttt{background} & d_j \ge f_{\text{ubiq}} \\
     \texttt{floor}      & d_j < f_{\text{ubiq}}
   \end{cases}

with :math:`f_{\text{ubiq}} =` ``ubiquitous_frac`` :math:`= 0.80` (default).

Interpretation:

**background mode** (:math:`d_j \ge 0.80`)
   A peak appears at this m/z even in cells lacking the tag, so a real ubiquitous
   background is present and must be exceeded. The centre is its measured median.

**floor mode** (:math:`d_j < 0.80`)
   The tag is absent from most background cells, so its background never cleared
   the instrument's peak-detection floor — centroiding already discarded it and
   the observable centre is zero.

.. note::
   The mode is decided globally rather than per file because it is a step
   function of :math:`d_j`: a tag near the cutoff would otherwise change regime
   between acquisitions, and its threshold with it. Per-file adaptation is
   retained in :math:`c_{bj}` and :math:`\sigma_{bj}` below.


.. _step-sigma:

Step 7 — Sigma
==============

``estimate_background_sigma(x, nonzero=True)``. Let
:math:`\mathcal{P}_{bj} = \{\, i : b(i) = b,\; x_{ij} > 0 \,\}` be the cells in
file *b* with a peak at tag *j*. Then

.. math::

   \tilde{m}_{bj} \;=\; \operatorname{median}_{\,i \in \mathcal{P}_{bj}} x_{ij}

.. math::

   \sigma_{bj} \;=\; 1.4826 \cdot \operatorname{median}_{\,i \in \mathcal{P}_{bj}} \left| x_{ij} - \tilde{m}_{bj} \right| \;+\; \varepsilon

with :math:`\varepsilon = 10^{-12}`. The constant 1.4826 rescales the median
absolute deviation to a standard-deviation equivalent for a Gaussian.

.. important::
   :math:`\sigma_{bj}` is computed over **all cells in the file** (restricted to
   non-zero values), not over the negative set. The centre comes from the negative
   cells; the spread comes from the whole file.


Step 8 — Threshold and detection
================================

Combining the mode (global) with the centre and sigma (per file):

.. math::

   \tau_{bj} \;=\;
   \begin{cases}
     c_{bj} + k\,\sigma_{bj} & \text{mode}_j = \texttt{background} \\[4pt]
     k\,\sigma_{bj}          & \text{mode}_j = \texttt{floor}
   \end{cases}

where :math:`c_{bj}` is the background centre of Step 5 recomputed **within file
b**, and :math:`k =` ``K`` :math:`= 2.0`.

Both branches are the same detection-limit form
:math:`\tau = \text{centre} + k\sigma`; only the measurable centre differs, being
zero in floor mode.

A cell is called positive for a tag when its intensity exceeds that file's
threshold:

.. math::

   P_{ij} \;=\; \mathbb{1}\!\left[\, x_{ij} > \tau_{b(i)\,j} \,\right]


Step 9 — Cell-type assignment
=============================

``classify_lineage(present, marker_names, marker_type)``.

**Lineage positivity.** A cell is positive for a lineage if **any** of that
lineage's tags fired:

.. math::

   L_{it} \;=\; \bigvee_{j \,:\, t_j = t} P_{ij}
   \qquad\qquad
   n_i \;=\; \sum_t L_{it}

**Lineage call.**

.. math::

   \text{lineage}_i \;=\;
   \begin{cases}
     \texttt{unassigned} & n_i = 0 \\
     t \text{ such that } L_{it} = 1 & n_i = 1 \\
     \texttt{dual} & n_i \ge 2
   \end{cases}

``dual`` denotes cross-lineage co-expression (cross-reactivity). Co-expression of
several tags **within** one lineage is not a dual — all three neuronal tags mark
the neuron lineage in their own right.

**Subtype.** For cells assigned to the lineage of the subtype marker
(SYN-I → neuron):

.. math::

   \text{subtype}_i \;=\;
   \begin{cases}
     \texttt{PS}  & \text{lineage}_i = t_{\text{SYN-I}} \;\wedge\; P_{i,\text{SYN-I}} = 1 \\
     \texttt{NPS} & \text{lineage}_i = t_{\text{SYN-I}} \;\wedge\; P_{i,\text{SYN-I}} = 0 \\
     \texttt{""}  & \text{otherwise}
   \end{cases}


Parameters
==========

.. list-table::
   :header-rows: 1
   :widths: 22 12 66

   * - Parameter
     - Value
     - Role
   * - ``TOL`` (:math:`\delta`)
     - 0.01 Da
     - half-width of the extraction window
   * - ``K`` (:math:`k`)
     - 2.0
     - threshold multiplier on :math:`\sigma`
   * - ``ubiquitous_frac`` (:math:`f_{\text{ubiq}}`)
     - 0.80
     - background/floor cutoff on :math:`d_j`
   * - ``NEG_SET``
     - ``"marker"``
     - negative set = cells not seed-positive for that tag
   * - ``min_cells`` (:math:`n_{\min}`)
     - 20
     - minimum non-zero negatives before :math:`c_j` is trusted
   * - ``N_DROP``
     - 299
     - leading acquisitions dropped from the PartII file
   * - ``rel`` (:math:`r`)
     - 0.10
     - blank-spectrum cutoff, fraction of the file's median peak count
   * - ``nonzero``
     - ``True``
     - all medians/MADs computed over non-zero values only


.. _not-used:

Implemented but not used in this run
====================================

The module contains further estimators that the notebook does not call. They are
listed here so the documented pipeline is not confused with them:

* ``otsu_threshold`` — bimodal split of the log1p non-zero histogram
* ``mixture_threshold`` — 1- vs 2-component Gaussian mixture with a BIC and
  Ashman's-D test for bimodality
* ``blank_window_stats`` / ``blank_window_threshold`` — background measured in
  signal-free m/z windows flanking each tag
* ``isotope_check`` — M+1/M isotope-ratio confirmation
* ``negative_cell_threshold`` — the single-batch threshold routine
  (``auto_marker_threshold_per_batch`` calls ``_negative_centres`` directly)
* ``classify_ssms`` — the earlier two-level scheme, superseded by
  ``classify_lineage``
* ``tic_normalize``, ``pqn_normalize``, ``batch_correct_nonneg`` — normalisation
  variants; no per-cell normalisation is applied in this pipeline


Discrepancy in the executed notebook
====================================

``GLOBAL_MODE = False`` is defined in Cell 2 but never read. The string occurs
exactly once in the notebook — in its own assignment — and Cell 5 calls
``auto_marker_threshold_per_batch`` unconditionally. That function always decides
the mode globally (Step 6).

**The documented run therefore used the global mode decision**, and the flag has
no effect. Either wire the flag to a branch or remove it.

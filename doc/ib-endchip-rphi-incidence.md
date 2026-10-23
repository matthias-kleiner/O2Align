---
name: ib-endchip-rphi-incidence
description: "IB end-chip dy-vs-Z excess (Lr0 |Z|>12) is an rphi-only term linear in tan(phi_inc) with no dz counterpart, i.e. not a rigid-body/geometric DOF; doc/handout_endchip_dy_vs_z.md's pooled RY \"endChipTilt\" has the wrong shape"
metadata:
  node_type: memory
  type: project
  originSessionId: 0fc83eca-09e7-412a-8690-fdc80e176fef
  modified: 2026-10-22T21:25:47.017Z
---

Found 2026-10-22 from the solA2v10 and v10n checkResid `res` trees (pp B0 562862, OO B+ LHC25ae):
- The Lr0 dy excess at |Z|>12 (+13 µm at +13) sits in the outer half of chip0/chip8, populated only by |zv|>~11 tracks. The same tracks give a phi-independent "vortex" vertex dy of −20…−40 µm (the page-1 vertex arch). Same in V10n and MP2, and in B0 and B+.
- Per chip, dy = a + b·tan(phi_inc) (the turbo tilt makes <tan>≈0.27, hence a phi-averaged mean). Lr0 b runs from +25 µm (c0) to −7 µm (c8); Lr1/Lr2 have similar, dataset-dependent terms.
- The dz counterpart that any normal displacement or rotation must give (dz = δX·tanλ_inc) is absent. The calibration check holds: fitted TY changes V10n→MP2 show up equally in the dy- and dz-based estimates (slope 0.72). Rigid-body fits per chip give a consistent ωX (in-plane rotation) and roughly a consistent pitch from dy and dz, but δX0 and ωZ disagree strongly.
- The v10 fit has well-determined chip TY/RY (errors 0.5 µm / ~45 µrad). Chip0 RY = +411 µrad made the checkResid in-plane slope worse (+108 → +441 µrad), because RY absorbs part of the unmodelled term. The handout's "underdetermined" claim is wrong.

**Why:** a non-rigid, rphi-only incidence effect (cluster centroid vs incidence / data–MC response difference, or an equivalent in-plane rphi stretch; primaries cannot separate the two). No rigid-body DOF, weighting or pooled RY can fix it. The MC reference (checkResidMCRef) is flat.

**How to apply:** before adding DOFs, split by cluster size/pattern (needs checkResid to store it) and compare with MC. If it is cluster-level, the remedy is a y-only tan(phi_inc) term per chip (or a dictionary/centroid correction), or de-weighting the end hits. Lr0 is an extrapolation residual in checkResid (the outward prediction holds only the vertex), so per-layer attribution among Lr0–Lr2 is not unique. Scratch tools: flatIB.C + numpy fits (sensDump.C for alpha/zc/refX), rebuild if needed. Related: [[its-cluster-x-fix]].

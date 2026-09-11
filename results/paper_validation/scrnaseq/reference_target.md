# Reference Target Specification: Pancreatic Endocrinogenesis E15.5

## 1. Study and Data Origin
- **Study:** Bastidas-Ponce et al. (Development 2019; DOI: 10.1242/dev.173849)
- **Repository:** https://github.com/theislab/pancreatic-endocrinogenesis
- **Dataset File:** `data/Pancreas/endocrinogenesis_day15.h5ad`
- **Cell Count:** 3,696 cells
- **Feature Count:** 27,998 genes
- **Developmental Stage:** Murine embryonic day 15.5 (E15.5)

## 2. Author Ground-Truth Annotations
The AnnData object includes author cell-type annotations in `obs['clusters']`:
1. **Ductal** (N = 916): Pancreatic ductal progenitor pool
2. **Ngn3 low EP** (N = 262): Early endocrine progenitor with low *Neurog3* expression
3. **Ngn3 high EP** (N = 642): Committed endocrine progenitor with peak *Neurog3* and *Hes6*
4. **Pre-endocrine** (N = 592): Transitional cells expressing *Pax4*, *Neurod1*, *Fev*
5. **Beta** (N = 591): Insulin-producing mature endocrine cells (*Ins1*, *Ins2*, *Pdx1*)
6. **Alpha** (N = 481): Glucagon-producing mature endocrine cells (*Gcg*, *Arx*)
7. **Delta** (N = 70): Somatostatin-producing cells (*Sst*, *Hhex*)
8. **Epsilon** (N = 142): Ghrelin-producing cells (*Ghrl*)

## 3. Canonical Marker Hierarchy
- Progenitors: *Sox9*, *Hes1*
- Early commitment: *Neurog3* (low), *Fev*
- Peak commitment: *Neurog3* (high), *Hes6*
- Pre-endocrine transition: *Pax4*, *Neurod1*
- Beta lineage: *Ins1*, *Ins2*, *Pdx1*, *Nkx6-1*
- Alpha lineage: *Gcg*, *Arx*, *Mafb*
- Delta lineage: *Sst*, *Hhex*
- Epsilon lineage: *Ghrl*

## 4. Endocrine Differentiation Trajectory
The canonical developmental trajectory proceeds along:
Ductal / Ngn3 low EP → Ngn3 high EP → Pre-endocrine → Beta / Alpha branching.
Evaluated via PAGA (partition-based graph abstraction) and DPT (diffusion pseudotime) from the early progenitor root.

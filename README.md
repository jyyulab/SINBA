
<!--README.md is generated from README.Rmd. Please edit that file -->

# SINBA

<!-- badges: start -->
<!-- badges: end -->

## Background

*S*ynergy *i*nference by data-driven *n*etwork-based *B*ayesian
*a*nalysis(**SINBA**) is used to bridge the novel discoveries from
network-based systems biology to reverse pharmacology. SINBA will
tremendously reduce the time and resource cost of conventional
combination screening. SINBA also integrated a blood-brain barrier
penetrated drug-target database, which will expedite the drug discovery
for CNS tumors.

## Package structure

<center>

<figure>
<img src="[https://github.com/jyyulab/SINBA/blob/main/SINBA_1.0.png]"
alt="SINBA(1.0)" />
<figcaption aria-hidden="true">SINBA(1.0)</figcaption>
</figure>

</center>

## Installation

------------------------------------------------------------------------

Require R \>= 3.5.0.

1.  install from github

``` r
devtools::install_github("jyyulab/SINBA",auth_token = "your auth_token", lib="your lib path")
```

2.  install from local file

``` r
pkg.dir <- "/Volumes/project_space/yu3grp/software_JY/yu3grp/git_repo/SINBA" #"/research/projects/yu3grp/software_JY/yu3grp/git_repo/SINBA_1.O"
devtools::install_local(sprintf("%s/SINBA_0.0.1.0.tar.gz",pkg.dir),lib="your lib path")
```

## Documentation

Instruction, documentation, and tutorials can be found at:

- <[https://jingl87.github.io/SINBA/](https://jingl87.github.io/SINBA/index.html)>

## Features

- **Comprehensive Drug-Gene Interaction Resources:** Includes six
  curated drug-gene interaction databases, featuring the blood-brain
  barrier penetrant drug database (`BPdb`), which supports the selection
  of CNS-specific drugs.  
- **Integrated Gene-Gene Interaction Network:** Facilitates the
  inference of disease or subgroup-specific driver genes. These drivers
  may not have identifiable genetic alterations and can be “hidden”
  regulators affected by epigenetic, post-transcriptional, or
  post-translational mechanisms.
- **Synergy Discovery Platform:** Harnesses the strengths of *in silico*
  prediction and *high-throughput screening* to identify novel
  synergistic drug pairs.

## Credits and historical remarks

The companion packages and data portals

- [NetBID](https://github.com/jyyulab/NetBID)
- [SJARACNe](https://github.com/jyyulab/SJARACNe)
- [Medulloblastoma SINBA
  portal](https://yulab-stjude.shinyapps.io/SINBA_MB/)
- [Medulloblastoma patient
  scRNA-seq](https://scminer.stjude.org/study/GSE155446)
- [Cerebellum Glutemetergic lineage developmental
  reference](https://scminer.stjude.org/study/Glutamatergic-lineage)

## References

<table>
<tr>
<td>
A. Khatamian, E. O. Paull, A. Califano and J. Yu Bioinformatics 2019:
SJARACNe: a scalable software tool for gene network reverse engineering
from big data.
<a href='https://academic.oup.com/bioinformatics/article/35/12/2165/5156064'>PMCID:
PMC6581437</a>
</td>
</tr>
<tr>
<td>
X. Dong, L. Ding, A. Thrasher, X. Wang, J. Liu, Q. Pan, et al.Nat Commun
2023: NetBID2 provides comprehensive hidden driver analysis.
<a href='https://www.ncbi.nlm.nih.gov/pubmed/37142594'>PMCID:
PMC10160099</a>
</td>
</tr>
</table>

# OncoSplicing

The code related to the database OncoSplicing (www.oncosplicing.com).

The OncoSplicing database was initially released in 2019 for visualization of survival associated alternative splicing for cancers in the TCGA SpliceSeq project, which was presented in our previous paper (https://doi.org/10.1038/s41388-019-0910-7). The code under **OncoSplicing/v1/** were two pieces of code related to survival analysis and consuse clustering analysis used in the Pan-cancer analysis. 

Recently, we have updated the OncoSplicing database. Code under **OncoSplicing/v2/** were used in visualization of clinically relevant alternative splicing in the SpliceSeq and SplAdder project.

Given the similarity of code used in these two projects, the code related to the SplAdder project were more representative, including:
1. spladder_kmplot.R  
  code to perform Kaplen Meire plot for overall survial (OS) or progression free interval (PFI) 

2. spladder_tnplot.R  
  code to plot PSI difference between tumour and adjacent normal or GTEx normal tissues       

3. spladder_ciplot_survival.R   
  code to perform Kaplen Meire plot for OS, PFI, DFI (disease free interval) or DSS (disease specific survial)

4. spladder_ciplot_boxplot.R  
  code to plot PSI difference between two groups of a clinical indicator 

5. spladder_panplot.R   
  code to plot PSI distribution in pan-cancer view       

6. spladder_pancox.R  
  code to plot result of CoxPH analysis in pan-cancer view 

7. spladder_pandiff.R  
  code to plot PSI difference in pan-cancer view 


Code related to the SpliceSeq project can also be obtained under the directory OncoSplicing/v2/spliceseq/, including:
1. spliceseq_kmplot.R
2. spliceseq_tnplot.R
3. spliceseq_ciplot_survival.R
4. spliceseq_ciplot_boxplot.R
5. spliceseq_panplot.R
6. spliceseq_pancox.R
7. spliceseq_pandiff.R

cite:

![1. Zhang Yangjun, Yan Libin, Zeng Jin, Zhou Hui, Liu Haoran, Yu Gan, Yao Weimin, Chen Ke, Ye Zhangqun, Xu Hua*. Pan-cancer analysis of clinical relevance of alternative splicing events in 31 human cancers. Oncogene.2019 Oct;38(40):6678-6695.](doi:10.1038/s41388-019-0910-7)

![2.OncoSplicing: an updated database for clinically relevant alternative splicing in 33 human cancers](www.oncosplicing.com)

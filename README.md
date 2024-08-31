# OncoSplicing

The code related to the database OncoSplicing (www.oncosplicing.com).

The OncoSplicing database was initially released in 2019 for visualization of survival associated alternative splicing for cancers in the TCGA SpliceSeq project, which was presented in our previous paper (https://doi.org/10.1038/s41388-019-0910-7). The code under **OncoSplicing/v1/** were two pieces of code related to survival analysis and consuse clustering analysis used in the Pan-cancer analysis. 

The code under **OncoSplicing/v2/** were used in visualization of clinically relevant alternative splicing in the SpliceSeq and SplAdder project.

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


Code related to the SpliceSeq project can also be obtained under **OncoSplicing/v2/spliceseq/**, including:
1. spliceseq_kmplot.R
2. spliceseq_tnplot.R
3. spliceseq_ciplot_survival.R
4. spliceseq_ciplot_boxplot.R
5. spliceseq_panplot.R
6. spliceseq_pancox.R
7. spliceseq_pandiff.R



Recently, we have updated the OncoSplicing database. Code under **OncoSplicing/v3/** were used in visualization of the relationship between RBPs and alternative splicing in cancers.
1. mapas_plot.R
    code to perform tracks plot to visulizaion the relative location of RNA binding motifs or eCLIP-seq peaks 
    of RBPs to the structured AS event
   
2. encode_plot.R
    code to perform tracks and sashimi plot to visulizaion signal difference between normal contorl and RBP 
    perturbation samples
   
3. encode_rbp.R
    code to perform volcano plot for visulizaion of the PSI difference and significance of the queried AS event 
    in different RBPs perturbation experiments.
   
4. encode_event.R
    code to perform volcano plot for visulizaion of  the PSI difference and significance of all regulated AS 
    events after perturbation of the queried RBP
   
5. coexp_table.R
     code to perform correlation analyses between the queried AS event and RBP in TCGA cancers or GTEx tissues
   
6. coexp_plot.R
     code to perform scatter plot for visulizaion of correlation result between the queried AS event and RBP in 
     a cancer/tissue type
    
7. script/encode_vis.R
     code to perform tracks plot that was refferenced by the mapas_plot.R and encode_plot.R
     this code was a modification of BioSeqUtils, full code can be found at https://github.com/junjunlab/BioSeqUtils.
    
9. script/geom_arch.R
     code to perform tracks plot that was refferenced by the mapas_plot.R and encode_plot.R



citation:

1. Zhang Yangjun, Yan Libin, Zeng Jin, Zhou Hui, Liu Haoran, Yu Gan, Yao Weimin, Chen Ke, Ye Zhangqun, Xu Hua*. Pan-cancer analysis of clinical relevance of alternative splicing events in 31 human cancers. Oncogene.2019 Oct;38(40):6678-6695. 

2. Zhang Yangjun, Yao Xiangyang, Zhou Hui, Wu Xiaoliang, Tian Jianbo, Zeng Jin, Yan Libin, Duan Chen, Liu Haoran, Li Heng, Chen Ke, Hu Zhiquan, Ye Zhangqun, Xu Hua*. OncoSplicing: an updated database for clinically relevant alternative splicing in 33 human cancers. Nucleic Acids Res. 2021 Sep 23:gkab851. doi: 10.1093/nar/gkab851.

# PROJECT: Simulation of GC Content Evolution

The association between GC content and protein disorder percentage in bacteria was uncovered 
and whether GC content defines the percentage of IDPs remains known. To tackle the problem, 
this project uses a simulation model to introduce neutral mutations and compare the resulting 
relationship.

# Materials and Methods

**Data Preparation**

20 bacterial pathogens were selected because of their association with global health loss according to a recent study (Ikuta et al., 2022). All bacterial species selected are classified as gram-negative, in case the structural differences between gram-positive and negative bacteria lead to confusion in analysis. 

  All data consisting of gene and protein sequences were downloaded from the database proGenomes3 (https://progenomes.embl.de/genome.cgi), which incorporated well-annotated proteins and corresponding coding sequences. For each species, an example sample was downloaded for analysis.
  
**Prediction of Intrinsically Disordered Proteins and GC Content Calculation**

With protein sequence, IUPred2A (https://iupred2a.elte.hu/plot) was used to predict intrinsically disordered proteins based on the estimation of energy generation to identify ordered and disordered residues (Mészáros et al., 2018).
IUPred2A has smoothed the energy predicted for each residue in the amino acid sequence through a specific window size (Mészáros et al., 2018), so here the protein disorder percentage was calculated by a single residue. Residues that have an IUPRED score over 0.5 are regarded as disordered (Mészáros et al., 2018), and they are summed for a subsequent division for percentage calculation.
  Corresponding protein-coding sequences were counted GC as stored values.
  
**Linear Regression and Identification of Influential Outliers**

The association between GC content and protein disorder was modeled using linear regression. For applying linear regression, many assumptions are considered and tested like exhibiting linear relationship and normality of residuals. Cook’s Distance (Cook, 2000) was used to identify influential outliers. 

Additionally, R-squared, p-value in regression analysis, and Pearson or Spearman correlation coefficient served as fitting evaluations. 

**Simulation: Introduction of neutral mutations**

According to the standard codon table, one amino acid can be encoded by different codons. Using codon-to-aa and aa-to-codon dictionaries, gene sequences were sliced by 3 bases (if not divisible, the gene would be discarded) for matching corresponding amino acids and finding their neutral mutations.

  The simulation model has three modes. First, set the different probability of neutral mutation (exclude original codon), influencing the coverage of mutation. Three cases of randomly introducing mutations in low probability (P=0.05), medium probability (P=0.5), and high probability (P=0.95) were simulated around 500 times to identify appropriate linear regression. The second mode is to manually introduce mutations that lead to the best and worst fitting. Different codons encoding the same amino acid were sorted by GC content and then the highest and lowest GC content for each sample can be calculated. Referring to the visualization of the original data, the worst case can be manually set extreme GC content to deviate from the regression line, and the best case is the opposite. The third mode detaches the original DNA sequence, considering all possible codons equally. In that case, according to every amino acid of the given sequences, all possible codons including the original codon were selected randomly to observe the theoretical results.
After changing GC content, the relationship with protein disorder was analyzed with various evaluations.


# Sharing core codes
My complete simulation model can achieve many functions (3 modes) and I shared the core codes (other R codes are not shown) publicly. 
They are not so formally or well packaged. I just run it for my personal convenience.

But I still do that because of its idea and the potential meaning for the evolutionary field.

Hope you can understand it, and feel free to contact me (**yutongh.20@intl.zju.edu.cn**) with any ideas or questions. 

# BIRdCAGE
This repository collects Matlab code I used during my PhD to simulate complex chromosomal rearrangements in cancer genomes. An overview of the methodology and presentation of some results can be found here: https://www.researchgate.net/publication/323656790_Inferring_cancer_genome_evolutions_from_paired-end_and_copy_number_data 

# Motivation:

Cancer genomes are characterized by complex structural and copy number alterations , which have been accumulated in the tumour cell population over many years. There is a great interest in trying to reconstruct the resulting aberrant karyotypes, as well as inferring the underlying molecular mechanisms. In this study, we implemented an algorithmic approach for the reconstruction of cancer digital karyotypes, generalizing different classes of DNA repair mechanisms through the use of two main operations: Cut and Repair (CR) and Break-Induced Replication (BIR). The algorithm models observed paired end reads as a series of rearrangement operations, thus constructing multiple steps of structural and copy number modifications of the genome configuration. The agreement of each resulting structures with the read depth and allelic data from the tumour sample is then assessed, assigning a score to each solution; these two types of information provide an efficient way of selecting the evolutionary hypotheses which best explain the tumour data, as indicated by our simulations. We tested our evolutionary model on both simulated and real data from a breast cancer sample. Results demonstrates the power of this approach in identifying reliable explanations. By modelling a large collection of known rearrangement mechanisms, our model has a potential in shedding light on the repair mechanisms leading to complex rearrangements.

# Basic usage:


The main function to construct a set of evolutions is

Struct=Gen_simul_BIR10281_MF_fast_TD_track3_modified(genA,rearredge,telsmat,MFvet,varargin);


The function takes 4 mandatory input variables:

1) genA. A single row matrix representing the starting genome, where each integer represent a DNA segment and zeros separate adjacent chromosomes. It has structure [0 chromosome_1 0 chromosome_2 0 ... 0 chromosome_n 0] .   For example, two homologous copies of a chromosome made up of 4 segments are represented as genA=[ 0 1 2 3 4 0 1 2 3 4 0]. A genome containing two non-homologous chromosomes divided into 4 and 3 segments, respectively, can be represented as genA=[ 0 1 2 3 4 0 5 6 7 0] . Non homologous chromosomes are not allowed to share any integer symbol.

2) rearredge. 2 by n matrix of integers, representing the n somatic connections to be added during evolution. Every row takes the form [+/- x, +/-y]  , with x<=y, where x (resp. -x) represent the left (resp. right) side of node (i.e. DNA segment) x, and similarly y and -y represent the left and right side of node y, respectively.   Connections will be added according to the order of rows (row 1, then row2, etc).  For example, [3,4;-2,2] represent two connections: one connecting the right side of node 3 to the right side of node 4, the other connecting the left side of node 2 to the right side of the same node.  It is important to notice that the set of integers in genA must include the set of integers in rearredge.  

3) telsmat.   2 by v matrix of integers, where v is the number of distinct chromosomes in the original genome.  Every row takes the form [x y] , with x<=y, where x and y represent the nodes carrying the telomeres of a particular chromosome.  For example, if we have genA=[ 0 1 2 3 4 0 5 6 7 0], we define telsmat as [1,4;5,7]: one set of telomere nodes [1,4] defining the first chromosome, and another set [5,7] defining the other one.

4) MFvet. A single row matrix of the same size as genA, useful to keep track of the copy number changes of homologous copies. It has the same structure as genA, but all DNA segments of the same chromosome are defined as the same single integer number. For example, given genA=[ 0 1 2 3 4 0 1 2 3 4 0] we can define MFvet as [0 1 1 1 1 0 2 2 2 2 0] . This would allow to keep track of allele specific copy number changes. The algorithm  adds rearrangement operations to MFvet, generating a new matrix where the multiplicity of each integer reflect the copy number changes of that particular chromosome.  Alternatively, it can be simply defined as MFvet=genA. 


optional variables:

"random"  : used to decide whether to randomly sample the space of evolutions

'fr_el' : at the end of the evoltion, decides whether to keep (value of 0) or eliminate (value of 1) rearranged chromosomes carrying somatic telomeres [default=1] 
'scorevet'  : vector of length 3 assigning a weight for each class of operations


#######################


output:

A Structure array with the following fields:

Evolutions: cell array of consistent evolutions. Each cell contains an ordered collection of matrices, each representing an evolutionary stage of the same evolution. 

Rearr_matrix:  cell array of matrices, each corresponding to a different evolution. Rows represent evolutionary stages, columns are position in the final rearranged genome. For each row (evolutionary step), zero and ones indicate whether the corresponding DNA segment has already been implicated in a rearrangement (one) or not (zero). 

Evolutions: cell array of MFvet matrix evolution. Each cell correspond to an evolution, and contains the final MFvet matrix obtained by applying all rearrangement operations. 


Example to run Gen_simul_BIR10281_MF_fast_TD_track3_modified


genA= [ 0 1 2 3 4 0];
rearredge=[-2,2;3,3];
telsmat=[1,4];
MFvet=[ 0 1 1 1 1 0];


Struct=Gen_simul_BIR10281_MF_fast_TD_track3_modified(genA,rearredge,telsmat,MFvet)

Struct = 

Evolutions: {[0 1 2 3 4 0]  [0 1 2 2 3 4 0]  [0 1 2 2 3 -3 -2 -2 -1 0 4 0]  [0 1 2 2 3 -3 -2 -2 -1 0]}
     Final_genomes: [0 1 2 2 3 -3 -2 -2 -1 0]
      Rearr_matrix: [3x10 double]
    Allelic_matrix: [0 1 1 1 1 -1 -1 -1 -1 0]


%Visualising evoutions:

Evols={Struct.Evolutions};

%Step 2 of evolution 1:

Evol1{1}{2}

ans =

     0     1     2     2     3     4     0



#We can also generate RANDOM evolutions:
[evolgood,same_evol,genevolfilt,mapfiltnew2,CNmat,CNtarget,best_score,ind_best_score,scorefilt,ind_scorefilt,ind_scorefilt_p,statsmat,rearrmat,opmatnew,opindmat,cont,rmat,MFgood]=cn_finder_b91b_all_perms_MF_new_score_10281_2_nfr_comp_33_track(nevents,nsamples,fast,varargin)



[genevolfinal,usedbpmat,genfinalfilt,genfinalfilt2,evol_mapfinal,evol_mapfilt,evol_mapfilt2,genfinalfiltconn,opmatfinal,ind1,ind2,templvet,opindmat,opvet2,vetfilt3,n_empty,ind22bpmatfinal,CD2,replvet,breakvet2,genevolfinalfr,which_reused,scorevolfinal,scorefinalfilt2]=Gen_evol_BIR10281_noMF_valsmat(genA,rearredge,telsmat,valsmat,varargin)

input:
valsmat: array of copy number information, where position "i" contains copy number data for the ith DNA segment.

output:
genevolfinal: array containing all genomes. Columns corresponds to evolutionary stages (column 1 is the reference genome), rows correspond to different evolutionary choices 
usedbpmat : matrix indicating the use of breakpoints. 1 indicates "used", zero indicates "not used"
genfinalfilt:   array containing filtered genomes 
genfinalfilt2:  array containing filtered genomes 
evol_mapfinal
evol_mapfilt
evol_mapfilt2
genfinalfiltconn
opmatfinal
ind1
ind2
templvet
opindmat
opvet2
vetfilt3
n_empty
ind22bpmatfinal
CD2
replvet
breakvet2
genevolfinalfr
which_reused
scorevolfinal
scorefinalfilt2


###########
[score,scorevet,stats,fityCN,Xfig,Yfig,zeromat] =cn_score_2_weight_d(cn,valsmat,r_tels,do_mean,fig)


[evol_filt2,evol_filt2b,MFevol_filt2,MFevol_filt2b,trackevol_filt2] =map2evol_b_MF_track(genevol,MFevol,mapfilt2,trackevol)

[genevolfinal,genfinalfilt2,evol_mapfinal,evol_mapfilt2,opmatfinal,vetfilt3,ind22bpmatfinal,CD2,replvet,breakvet2,genevolfinalfr,which_reused,scorevolfinal,scorefinalfilt2,MFmatfinal,MFmatfilt2,trackevolfinal,trackmatfilt2,evol_filt2]=Gen_simul_BIR10281_MF_fast_TD_track31(genA,rearredge,telsmat,MFvet,varargin)






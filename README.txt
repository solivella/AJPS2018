########################################
##      Replication archive for       ##
##      "Tree-based models for        ##
##      political science data"       ##
##      By J. M. Montgomery           ##
##        (jacob.montgomery@wustl.edu)##
##      and S. Olivella               ##
##          (olivella@unc.edu)        ##
##                                    ##
##      This version: February 2017   ##
##                                    ##
########################################

README.txt


This folder contains all replication materials
for "Tree-based models for political science
data". The following is a description of all
included files.

All analyses were conducted using the following
software versions: R 3.3.1; Plyr 1.8.4;
Lattice 0.20-34; doMC 1.3.4, gbm 2.1.1,
BayesTree 0.3, rpart 4.1.  


- (1)Figures_1_2_3.R:
	R code used to produce Figures 1, 2 and 3
    	in the main text.

- (2)Synth_data_for_5_1.R:
	R code used to generate the synthetic
   	 data used in section 5.1 of the main
    	 text.

- (3)CVFunction.R
	R code defining the cross-validation
	function used to select best models
	for contest described in section 5.1.

- (4)Analyses_section_5_1.R:
	R code used to define, estimate and cross-
    	validate all models discussed in section
    	5.1. SPECIAL NOTE: This code was developed
    	to be executed on a distributed-memory
    	cluster, with access to over 400 processing
    	cores. It takes about 8 hours to complete,
    	and relies on a batch file for submission
    	to the LSF scheduler (viz. contest.job; see
    	below). The code ultimately produces a
    	set of "best" models, selected on the
    	basis of a cross-validation exercise (see
   	main text for more details). The resulting
    	models are available as part of this replication
    	archive (see BestModels(finegrid).RData below).

- (5)Figures_4_and_SI1.R:
	R code used to produce figure 4
	in the main text and SI-1 in appendix. 

- (6)Analysis_sect_5_2_Tables_2_3.R:
	R code to conduct all estimations discussed
	in section 5.2 of the main text. This includes
	all weight estimations using GBM and BART,
	as well as all calculations necessary for
	producing tables 2 and 3 of the main
	text. Currently set up to run on multiple cores
	using the doMC backend.

- (7)Analysis_sect_5_3.R
	R code to perform the replications and
	extensions discussed in section 5.3 (i.e.
	poststratified turnout and vote intention
	models). SPECIAL NOTE: This code will produce
	very large intermediate objects (particularly, the
	object containing BART estimates of vote intention
	is over 9Gb large), so care should be exercised
	when running it. While not part of this archive,
	these intermediate objects are available upon
	request. Other objects needed by other script
	files are saved as .CSV objects, and are
	provided as part of this archive (see
	MRP_BART_Expansion_Preds.cvs below).

- (8)Figure_5.R
	R code used to produce figure 5 in the main
	text.

- BestModels(finegrid).RData
	Compressed R object files containing
	best (cross-validated) models in contest
	described in section 5.1.

- census-pums-pop-2000-04-08.dat
	Census counts by demographic group used in the
	replication and extension of GG. 
	Source: Ghitza, Yair; Gelman, Andrew, 2012, "Replication data for: Deep
        Interactions with MRP: Election Turnout and Voting Patterns Among Small
        Electoral Subgroups", hdl:1902.1/18580, Harvard Dataverse, V3

- Codebook.pdf
	PDF describing all variables used and documenting
	their source.

- contest.job
	LSF batch file that executes (4)Analysis_section_5_1

- cps2000-04-08-DKs.dat
	Survey responses to CPS questions about turnout and 
	demographics. 
	Source: Ghitza, Yair; Gelman, Andrew, 2012, "Replication data for: Deep
        Interactions with MRP: Election Turnout and Voting Patterns Among Small
        Electoral Subgroups", hdl:1902.1/18580, Harvard Dataverse, V3 

- dem.panel.RData
	Compressed R object file containing
	negative campaign events and election
	outcomes.
	Source: Blackwell, Matthew, 2013, "Replication data for: A
        Framework for Dynamic Causal Inference in Political Science", hdl:
        1902.1/19801, Harvard Dataverse, V2

- mapreduce_1.2.6.tar.gz
	Tarball of package used by Ghitza and Gelman,
	used for replication purposes. See file 
	(7)Analysis_sect_5_3.R for installation notes.

- MRP_BART_Expansion_Preds.csv
	Comma-separated values file containing
	post-stratified predictions of turnout
	and vote intentions, based on BART model.

- MRP-20120715.RData
	R object file containing results of Ghitza
	and Gelman (2012)’s analysis. Used for replication
	purposes. 

- README.txt:
	This file. 

- ResultsTable3.RData:
	R Data file containing object ``out’’,
	which contains results of bootstrapped
	standard errors reported in Table 3 off manuscript. 

- state-stats.dat
	Dataset of state-level covariates by year, including population, income,
	total population and nr. of votes cast.  
	Source: Ghitza, Yair; Gelman, Andrew, 2012, "Replication data for: Deep
        Interactions with MRP: Election Turnout and Voting Patterns Among Small
        Electoral Subgroups", hdl:1902.1/18580, Harvard Dataverse, V3 

- votechoice2000-04-08.dat
	Survey responses to NAES survey questions about
	vote intentions and demographics. 
	Source: Ghitza, Yair; Gelman, Andrew, 2012, "Replication data for: Deep
        Interactions with MRP: Election Turnout and Voting Patterns Among Small
        Electoral Subgroups", hdl:1902.1/18580, Harvard Dataverse, V3	
	
	
	




            

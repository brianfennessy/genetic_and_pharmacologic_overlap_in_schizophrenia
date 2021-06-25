# Mechanistic overlap between genetic risk and pharmacologic modulation of psychosis

Hello! Thank you for reading. We hope you enjoyed our paper, "Mechanistic overlap between genetic risk and pharmacologic modulation of psychosis." If you are looking to reproduce our results, there are several datasets that you need to download first. The data sets that we start with are listed below. A link is provided where you can download each data set. We have also indicated whether the data is open access, or if you must register with the respective organization to access the data.

While we provide the links to each data set, all the raw files you will need are included in the first section of our code, "1_inhouse_functions_and_setup." To recreate our results, download all nine sections and then work through each section in numerical order. Directories and subdirectories will be created for you. All you need to do is run the code in order!

===========================================

Data: VigiBase

Version: December 1st, 2016

Source: https://www.who-umc.org/whodrug/access-tools/download-area/

Access: Valid WHODrug Global subscription and a personal UMC username and password required

===========================================

Data: Side Effect Resource (SIDER)

Version: 4.1

Source: http://sideeffects.embl.de/download/

Access: Public

===========================================

Data: Anatomical Therapeutic Chemical Classification System (ATC)

Version: 2014

Source: https://www.whocc.no/atc/application_form/

Access: By request

===========================================

Data: DrugBank

Version: 5.1.1

Source: https://go.drugbank.com/releases

Access: Account required

===========================================

Data: Similarity Ensemble Approach (SeaChange)

Version: January 13th, 2014 release

Source: https://seachangepharma.com/

Access: Private - sent from company individually

===========================================

Data: RxNorm

Version: September 4th, 2018 release

Source: https://www.nlm.nih.gov/research/umls/rxnorm/docs/rxnormarchive.html

Access: Account Required

===========================================

Data: Medical Dictionary for Regulatory Activities (MedDRA)

Version: v20

Source: https://apps.meddra.org/selfservice/

Access: Organization access required

===========================================

Data: Psychiatric Genomics Consortium 3 (PGC3)

Version: PrediXcan - Dorsolateral Prefrontal Cortex

Source: https://www.medrxiv.org/content/10.1101/2020.09.12.20192922v1

Access: SNPs are public (supplementary table 3)

===========================================

Data: Schizophrenia Exome Meta-Analysis Consortium (SCHEMA)

Version: July 29th, 2019 release

Source: https://schema.broadinstitute.org/results

Access: Public

====================================================================================================================================================================================

In terms of the code for this project, our approach can be broken down into nine major sections including an initial section for in-house functions. These sections are outlined and briefly described below. The section numbers correspond to folders in the GitHub directory. In each section folder we have provided the code we used to produce the results in our paper.

Section One: In-house functions created by Alexander Charney and all raw data necessary to recreate our results.

Section Two: Standardizing drug, target, and mechanism data.

Section Three: Identifying ATC antipsychotic drugs.

Section Four: Identifying VigiBase propsychotic drugs.

Section Five: Identifying propsychotic and antipsychotic drug targets.

Section Six: Identifying significant propsychotic and antipsychotic drug targets (those hit more frequently than by chance) and evaluating overlap between propsychotic and antipsychotic significant targets (both with and without central nervous system drugs).

Section Seven: Identifying significant propsychotic and antipsychotic drug-target mechanisms (those occurring more frequently than by chance).

Section Eight: Evaluating overlap between propsychotic or antipsychotic drug targets with rare or common variants implicated in schizophrenia risk.

Section Nine: Creation of final results table.

====================================================================================================================================================================================


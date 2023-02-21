# PropaNet2
Propanet with KCCA GRN

Advantages of KCCA base GRN construction
- Complex relationship of multiple transcription factors (TFs) regulating multiple target genes (TGs) are modeled
- Different combinations of co-working TFs specific to certain conditions (e.g. time-series responses against heat stress) are modeled
- Our method outputs GRN in a feasible size of sub-networks, where regulatory relations from co-working TFs and co-regulated TGs are suggested as a subnetwork.

PropaNet discovers the dynamics of TF networks against external stress such as heat stress or cold stress.  
Given time-series **gene expression profile data** and **template network**, PropaNet retrieves **master regulatory TFs** in each time-point.


## Installation

## Dependency

## Usage
```
snakemake --configfile config.yaml --cores 1
```

## Input File Format
There are three types of input files, **gene expression profile data**, **DEG binary data**, **protein-protein interaction data**, **List of transcription factor**, and **template network**. 
Each data should take the following format:

#### Gene expression profile data
Time-series gene expression data has to be stored in a single tab-delimited matrix form.
There should be a header as the first row and it should take the following format:  

| Gene_id | Condition-time1  | Condition-time2  | Condition-time3  | Condition-time4  | ... |
| ------- | :--------------: | :--------------: | :--------------: | :--------------: | :-: |
| Gene1   | expression level | expression level | expression level | expression level | ... |
| Gene2   | expression level | expression level | expression level | expression level | ... |
| Gene3   | expression level | expression level | expression level | expression level | ... |
| ...     | ...              | ...              | ...              | ...              | ... |

###### Example)
```
gene	Heat-Shoots-0.25h	Heat-Shoots-0.5h	Heat-Shoots-1h	Heat-Shoots-3h
AT1G01010	4.91798979107909	4.81768901700114	5.67991373730826	6.02273932823551
AT1G01030	4.45806149779191	4.53272006731038	5.36916581768974	5.51771002546865
AT1G01040	3.8552645118378	3.87640592858555	4.42911094757011	4.83351626001844
AT1G01050	-2.84787294715136	-2.64691204175781	-3.05377896002984	-3.87740186446588
AT1G01060	-0.899222494190051	0.609781349298916	-0.985720858842949	1.57868637136184
AT1G01070	-0.899064779308411	0.609662493243061	0.985561143281074	1.57841350780699
AT1G01080	-0.89887322954867	-0.609596150621448	0.98554170116417	1.57833846949708
AT1G01090	0.898799238909556	0.609558025275848	0.985367888674293	-1.57830555053209
AT1G01100	-0.898776072706743	-0.609434153689444	-0.985195207884564	1.57828766511001
```

#### DEG binary indication data
DEG binary indication data is a matrix that indicates which genes in our expression data are differentially expressed genes(DEGs).  
If a value in the DEG binary indication data matrix is "1", it would mean that the corresponding gene at the same location in our expression profile data is differentially expressed. 
The following is our example DEG binary data.
```
gene	Heat-Shoots-0.25h	Heat-Shoots-0.5h	Heat-Shoots-1h	Heat-Shoots-3h
AT1G01010 1	1	-1	-1
AT1G01030	0	0	0	-1
AT1G01040	0	0	0	0
AT1G01050	0	0	0	0
AT1G01060	0	0	-1	-1
AT1G01070	0	0	0	-1
AT1G01080	0	0	0	-1
AT1G01090	0	0	0	-1
AT1G01100	0	0	0	0
...
```

#### Template network
The template network file should be comprised of 3 columns 
There should be _no header_ in the first row.

| Source gene  | Target gene  | Weight       |
| :----------- | :----------- | :----------- |
| Source gene1 | Target gene1 | Edge weight1 |
| Source gene2 | Target gene2 | Edge weight2 |
| Source gene2 | Target gene2 | Edge weight3 |
| ...          | ...          | ...          |


#### Protein-Protein interaction
The protein-protein interaction data file should be comprised of 2 columns 
There should be _no header_ in the first row.

| Source gene  | Target gene  |
| :----------- | :----------- |
| Source gene1 | Target gene1 |
| Source gene2 | Target gene2 |
| Source gene2 | Target gene2 |
| ...          | ...          |


---

## Workflow
1. 

1. Construct weighted network (Optional)
2. Extract optimal transcription factors (Main process of PropaNet)
3. Extract target genes
4. Make final result


## Results
The output is a network(in the form of an edge list) comprised of the resulting TFs/TGs that are found by PropaNet.

---


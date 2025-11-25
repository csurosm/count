This folder contains Count analyses over kingdom-level data sets in 269 Archaea

* tree269-gca-names-txids.txt gives the species names, NCBI assembly and taxonomy identifiers from Dombrowski et al 2020 [https://dx.doi.org/10.1038/s41467-020-17408-w]
* tree269-baker-alti-gtdb.tre reflects the tree from Baker et al 2025 [https://dx.doi.org/10.1038/s41564-025-02024-5]; tree269-baker-alti-names.tre gives species names 
* arc269v24-ar14asd24multi.csv is the arCOG-format annotation for domain occurrences across all genes 
* arc269v24-ar14asd24multi.txt is the corresponding input profile table for Count

datasets-sims-ED194-E114-D80-P75.countxml.gz is a saved session set from Count with analyses over the four kingdom-level genome subsets 
* ED194: Methanobacteriati/Euryarchaeota + Nanobdellati /DPANN
* E114: Methanobacteriati/Euryarchaeota
* D80: Nanobdellati/DPANN 
* P75: Proteoarchaeota
The sessions include simulations and inferences with the optimized models, and inferences with the optimizations over the simulated datasets. Also are include the filtered datasets over families with inferred ancestral presence at the root. 

cryptic-t94-min4.countxml.gz and cryptic-e114-min4.countxml.gz are analyses of cryptic inheritance: whether a family is present at node X when inferred from a small dataset A and a larger dataset B: here X is T28, M14 and H51 and B is their combined dataset T94. 
 
cryptic-e114-min4.countxml.gz shows cryptic inheritance at X=THM node when inferred from A=T94 or B=E114. 

You can open a Count session either from the GUI or from command-line with -load
java -jar CountXXV.jar -load datasets-sims-ED194-E114-D80-P75.countxml.gz & 

(Gzip compression is recognized in Count input files by the extension gz.)

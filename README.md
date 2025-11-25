Count 25 - beta
================

The Count is a software package for the evolutionary analysis of phylogenetic profiles
(homolog family sizes), or other numerical census-type characters (copy numbers) along a phylogeny.

Availability
------------
Count is written entirely in Java, and, thus, can be used in
different operating systems, including Mac OS X, Microsoft
Windows, and various Unix/Linux versions. The software is
packaged in a JAR file, and can be executed in Java versions
1.13 and above.

Download the latest JAR file from under Releases (-->)

Licensed under the Apache License, Version 2.0: 
see http://www.apache.org/licenses/LICENSE-2.0.

You can launch the Count GUI as `java -jar CountXXV.jar`, but various functionalities 
are also available from the command line as `java -cp CountXXV.jar ...`. See the command-line builder below.

Using Count
===========
Count implements a number of parsimony and probabilistic methods for inferring copy number evolution 
along a phylogenetic tree from **phylogenetic profiles**. A profile is formed by 
the number of copies at each terminal node (leaf). 

Errors 
------
It may happen that something goes wrong. Count displays a
window with the error message in such cases. If you think
that the error was caused by a programming bug, or you
would like to ask me for help with it, then send me an
email, with the detailed technical error message included in the
message body. The technical details are shown only when
you select the corresponding checkbox.

Sessions and the work area 
--------------------------
Count operates with **sessions**, defined by a set of terminal taxon names and one more phylogenies 
(trees) 
over the same leaf set.
The workspace has two browsers on the left: a model browser for trees and probabilistic models,
and a data browser for data sets (tables) 
and ancestral inferences by parsimony or posterior reconstruction. 
On the right of the work area, you have the selected browser item on display: a tree, a table, or 
an ancestral reconstruction.
On the top, you have the pull-down menu and a toolbar: actions refer to the currently displayed 
item in the work area (so they get enabled or disabled as you select different items in the browsers). 

Starting 
--------
As a first step, you can either
 * start a **new session** by loading a phylogeny in Newick format, 
	which will be the main tree associated with the session; or
 * load a table in a **no-tree session** or **import a table** of family-sizes (copy-numbers), 
and assume a simple initial tree (star tree, random tree, or UPGMA calculated from copy number profiles); or
 * **load a session** previously saved by Count (XML format).

Tables 
------
Count works with a family size **table**: 
a set of genomes described as a multiset of families their genes belong to. 
Such a table is tab-delimited text file that starts with a header line 
that lists the taxon names as column headers. The first column has the family names. 
Additional columns can be loaded as family annotation columns. 
For a tree with four leaves named A,B,C,D: 
<pre>
Family	A	B	C	D	
og001	1	0	3	1
og002	2	0	0	0
...
</pre>
You can either **load a table** in the tab-delimited format, or **import a table**
from COG-style csv files or MCL clustering output. 
You can derive further tables by **filtering the rows** of an existing table 
(double-click on a cell for condition on values in that column), 
or by converting it to a **binary** presence-absence table.
Clicking on a column header sorts the table rows by that column. 

Trees
-----
Every session has a main phylogeny, and can include other rooted trees over the same leaf set. 
You can either **load a tree** in Newick format, or **build a tree**, or modify one by hand, using the 
**edit tree** function. (Tres can be built with a simple UPGMA or Neighbr-Joining procedure from profiles, 
or set to random or star topology.
 Edit operations include rerooting, edge-contraction, and subtree-prune-and-regraft.)
 
 
Rate models 
------
Rate models equip the selected phylogeny's edges with a 
linear birth-and-death process defined by rates and length: loss, duplication and gain rates. 
The probabilistic model of evolution assumes 
a continuous-time linear birth-death process  with its own parameters on every edge: 
 * <var>t</var> length (is measured by normalizing for loss rate)
 * <var>e</var> per-copy loss rate (assumed to be =1.0 for scaling the edge lengths)
 * <var>d</var> per-copy duplication (assumed to be at most <var>e</var> usually)
 * <var>r</var> (copy-independent) gain rate 
The parameters define the rates 
for a set of <var>n</var> copies, births ((<var>n</var>&rarr;<var>n</var>+1) 
arrive by a rate of 
<var>r</var>+<var>d</var><var>n</var> 
and deaths by <var>e</var><var>n</var>.  


Phylogenies with multifurcations are OK, but a completely resolved model 
uses a binary tree. 

You can **load a rate model** (previously saved from Count), or make one by **model optimization**. 
(Models are optimized by maximizing their likelihood). 
You should set the observation (aka ascertainment) bias for your input table, 
depending on how the homolog families were constructed: 
whether they have a minimum of 0 (if even the families without any members are known), 
1 (if the table covers all genes from all genomes),
 or 2 members (if families correspond to alignments).
Make sure you filter your table for observation bias on the membership count (sort by clicking on the `#mem` header to find such a cell): 
families with at least 1 member, or families with at least 2 members, etc. if you work with rate models.

You can also **simulate a table** using a rate model. The simulation keeps the true ancestral history, 
which can then be compared to an ancestral reconstruction. 

Ancestral reconstructions
-------------------------
Given a table, you can do **ancestral reconstructions**. 
by either **numerical parsimony** (a generalized Wagner penalization), 
or by **Dollo parsimony** if the data is binary. 
If you have a rate model, then you can also reconstruct the ancestral genomes by **posteriors**, 
giving probabilities and expectations for family memberships; 
the reconstruction includes correction for the observation bias of minimum copies.

Saving and exporting
--------------------
You can **save** phylogenies, rate models, and tables, or **export** ancestral reconstructions.
 With the exception of phylogenies, all other data is in tab-delimited files. 
 You can also **save all** your sessions if you want to come back to them later.
 (The sessions are saved in an XML-format file.)

Tooltips 
--------
Tooltips give further instructions and informations on the displayed elements in the GUI.


Command-line usage and multithreading 
-------------------------------------
Time-consuming actions, like rate model optimization and ancestral reconstruction
 are launched in the background of the GUI using multiple threads, and you can follow the 
 progress in the GUI. 
 For model optimization, you may prefer to run it from the command-line: you have a  
 command builder in the GUI that constructs the one-line command or script 
 to do it. 
(Command-line execution writes log messages on the standard output (preceded by `#`), 
but the main result can be written in a separate output file if you set it in the 
command builder.)  


Limitations 
===========
Limitations in this beta version:
 *	optimization with rate variation across families (rate categories) is not supported
 *	optimization with homogeneous rates (same duplication or gain across lineages) is not supported
 *	help menu is not available

#MolNetwork

####A tool to genereate a molecular networks

###Overview:
This is a command line python tool for generating molecular networks. At the moment this is accomplished using Morgan (circular) feature counts with a radius of 4 bonds. These fingerprints are generated for each input molecule and then pairwise Dice similarity scores are caclulated for all compounds. Those scores which are above the specified cutoff are included as edges in a resulting tab deliminated output edge list file. This can be easily loaded into network visualization software such as cytoscape or gephi.

####Dependencies:
* Written for python 3.5 (it's time to move on all you 2.7 users...)
* [RDKit](http://www.rdkit.org/)

  This can be a bit of a pain to install from source. Highly reccommend using the [Annaconda Python](https://www.continuum.io/downloads) distribution and installing using the [conda recipe](http://www.rdkit.org/docs/Install.html).

* [progressbar2](https://pypi.python.org/pypi/progressbar2)

  There are potentially some long running loops in this program so it's nice to have some idea of when they're going to finish. progressbar2 does a pretty job of it. This is a bit of luxury dependency, might get removed down the line...

###Features:

Right now things are fairly basic. This is really meant to be a way to vizualize hits from screening campaigns but can be applied to general network generation use.

* Input files can be .smi files or .sd(f) files.
  * .smi file format has lines of : `smiles<tab>name,`
  * .sd or .sdf file will be detected from the extension
    *if no name attribute is in the sd file then the index will be used as the name for the output
* Full Pairwise Dice Scoring

  The default action is to compute all pairwise Dice scores (non-redundant and no self-loops) and output those above the default threshold score of 0.8

* Threshold can be changed w/ `--thresh <x>` (x=float)
* Any resulting singletons can be included as self loops w/ `--singletons`
* Negative SAR analysis

  The chemical space around hits can be explored if the full library of compounds is known. This is accomplished by the --negsar <full_library_file> argument. Hits are compared to all members of the full library and any scores over the threshold are included in the edge list file. An additional text file is generated (<outfile_hitMap.txt) that is a binary map of whether a node name is a hit (1) or a non-hit (0) that is of similar chemical composition to the hits.

###Experimental Featueres:

* Pairwise sampling for psuedo clustering

  For large datasets (~ >10k molecules) generation of the fingerprints and calculating all possible scores becomes infeasible. For *m* compounds there are *m*<sup>2</sup>/2 - *m* possible comparisons. By using the `--samplesize <n>` (n=int) argument you can randomly sample n compounds to compare against each compound in the infile. I'm not really sure if this is useful but it might be worth playing with.

###Future Feature Ideas:

* Multiprocessing
* Networkx / Gephi output support
* ???

###Usage:

#####Basic:

`python MolNetwork.py <infile> <outfile>`

#####With singletons included and a lowered threshold:

`python MolNetwork.py <infile> <outfile> --thresh 0.7 --singletons`

#####Negative SAR option used:

`python MolNetwork.py <infile> <outfile> --negsar <full_library_file>`

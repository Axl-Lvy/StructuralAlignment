# StructuralAlignment

This repo is a start of a protein alignment using structural informations project.

## Needleman-Wunsch algorithm

The Needleman-Wunsch algorithm is written in [Needle.py](Needle.py). There are differents variations of this algorithm, which are :
* gap, substitution and identification score
* scoring matrix (example : blosum62)
* affine gap penality, which can be from fixed scores or a scoring matrix

## Guide tree

In [Needle.py](Needle.py) you can also find a method based on UPGMA (explained on wikipedia [here](https://en.wikipedia.org/wiki/UPGMA)) to make an alignment of many structures.

## Structural information

At the end of [the notebook](Sequence_Profile_Alignment.ipynb) you can find a code that take superpose 2 proteins given in 2 differents orientations. This is the first step to use the structural informations to align the sequences.

## To continue

Using 2 proteins structures, the score of each step in the Needleman-Wunsch can be updated. For example, add a function of the distance between two amino acids.
# FragQuality
Measures the RMSD for each fragment at each protein structure's position and plots the lowest RMSD fragment for each position

## Description
This is a PyRosetta script that measures the RMSD for each of the 200 9-mer fragments (from the [Robetta Server](http://robetta.org/)) at every position in a protein's structure, finds every position's fragment's lowest RMSD value, then plots these (the lowest RMSD value at each positon) against the position number. This gives an indication of the quality of the fragments that will be used in folding simulations, and thus assists in determining the quiality of an Abinitio folding simulation.

## Requirements
1. It goes without saying that you need to download and compile [PyRosetta](http://www.pyrosetta.org/) to use this script.
2. You will need to install GnuPlot for the script to generate the PDF plot at the end:

`sudo apt install gnuplot`

## How To Use
1. This script, the protein structure (FILENAME.pdb), and the 9-mer fragment file (FRAGFILE) should all be in the same working directory.
2. Run using the following command:

`python3 FragQuality.py FILENAME.pdb FRAGFILE`

3. Output will be a file (RMSDvsPosition.dat) that contains the position (column 1) and its lowest RMSD value (column 2).
4. Output will also include a plot (RMSD vs POSITION) in PDF.

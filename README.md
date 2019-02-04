# Codon evolution

These scripts and datasets were created for the codon evolution rotation project in the Wilke lab, Winter 2018-2019

**Broad scope:** *in-silico* directed evolution experiment that utilizes pinetree to simulate transcription and translation rates. 

**Rationale:** presence of codon ramps at the 5’ and 3’ end of mRNAs could be a biological mechanism for controlling speed of protein production in highly translated genes (Tuller and Zur 2014). In addition, ramps are observed in Ashley’s program that may be of biological significance, or an artifact of the program (TASEP ribosomal model used in her translational modeling tool).

**Requirements**
* Python 3.7.0
 * NumPy v1.15 
 * Pandas
* Pinetree https://github.com/benjaminjack/pinetree

**Usage**
To begin simulation of a transcript with a slow rate of 0.5 (-m option) and a fast rate of 1.0 (-r option) for one generation (-g option), navigate to script directory and enter:

```bash
python3 evolve_variable_rates.py -g 1 -m 0.5 -r 1.0
```

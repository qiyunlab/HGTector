# About this example

The input dataset (`gsul.txt`) contains 100 protein-coding genes from the genome of _Galdieria sulphuraria_, a unicellular red alga found in hot sulphur springs. The genome of _G. sulphuraria_ and its HGT pattern was reported in the following paper:

> Schönknecht, G. _et al_. [Gene transfer from bacteria and archaea facilitated evolution of an extremophilic eukaryote](https://science.sciencemag.org/content/339/6124/1207.long). _Science_ **339**, 1207-10 (2013)

To run this example, execute:

```bash
hgtector search -i gsul.txt -o .
hgtector analyze -i gsul.tsv -o .
```

The sample output files are provided in `output`. They were generated on 2019-10-16, based on the NCBI nr database by time. The file `hgts/gsul.txt` was manually modified to include potential donors (a feature of more recent versions of HGTector).

Detailed instruction for running this example and interpreting outputs is provided in [first run](../doc/1strun.md).

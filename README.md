# BioJSON

**BioJSON** is a collection of rough but functional hacks to pull biological flatfiles kicking and screaming into the world of simple data exchange standards. 

It uses pyparsing as the recursive descent parsing engine.  Currently only the genbank file parser is present, parsers for BLAST, HMMER, and other formats will 
be added as needed

### Genbank Parser

I've used this to parse entire metazoan genomes, so I believe that it should work on the vast majority of NCBI [genbank](http://www.ncbi.nlm.nih.gov/genbank/) files.  
I've also used it to import genbank style files produced by the dna editing programs [Vector NTI](http://www.invitrogen.com/site/us/en/home/LINNEA-Online-Guides/LINNEA-Communities/Vector-NTI-Community/vector-nti-software.html) and [ApE](http://biologylabs.utah.edu/jorgensen/wayned/ape/)

Currently only the sequence and feature table are mapped to JSON in a nice way as these are the features of primary interest to synthetic biologists hacking sequences.


# TDNA_primer_design
Command Line tool to generate genotyping primers for Arabidopsis T-DNA lines

This module essentially duplicates the functionality of the SALK Signal T-DNA primer design tool (http://signal.salk.edu/tdnaprimers.2.html)
Running setup.py will download the necessary files including the Arabidopsis genome fasta and T-DNA insertion locations from SALK website.
After the initial file downloads the tool can be run without an internet connection.

A correspsonding REST API is currently partially implemented, available at https://wh725mphi9.execute-api.us-east-1.amazonaws.com/alpha . See below for endpoints. Full documentation coming soon.

/fasta 
  /by-gene
    /{AGI}
  /by-coords
    ?start=[start coord]&chr=[chromosome number]&bp_upstream=[basepairs upstream to include]&name=[optional name for sequence]&end=[end coords]&orientation=[+ or - strand]&bp_downstream=[basepairs downstream to include]
  /by-tdna
    /{SALK_ID}
/genes
  /get-coords
    /{gene}
/tdna
  /get-coords
    /{SALK_ID}
  /get-fasta
    /{SALK_ID}

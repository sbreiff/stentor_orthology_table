# stentor_orthology_table
For assigning Stentor gene models into ortholog groups

This repository contains a script for categorizing Stentor OrthoMCL results and organizing into a table, as well as text files containing lists of ortholog groups with specific phyletic patterns.
Resulting data can be found in Figure 2B, Figure S2A-B, and Table S2 of the [Stentor genome paper](http://www.cell.com/current-biology/fulltext/S0960-9822(16)31539-1).

Basic description of workflow:
1. OrthoMCL was run on the Stentor predicted proteins to organize them into ortholog groups
2. Lists of ortholog groups that comprised specific phyletic patterns were created
3. When a sequence didn't match any OrthoMCL group, orthologs were looked for in other ciliate proteomes
4. When a sequence didn't match any ortholog group, including in other ciliate proteomes, Stentor paralogs were looked for
5. A table was generated containing each Stentor gene, what ortholog group it belongs to (if any), and the group's phyletic pattern

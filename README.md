# proteinAnalysis2

`proteinAnalysis2` is an updated version of `proteinAnalysis (version 1)` Python package developed at **Biochemistry & Bioinformatics Lab, Department of Applied Sciences, Indian Institute of Information Technology Allahabad (IIIT-A), Devghat, Jhalwa, Prayagraj-211015, U.P., India**.

This package is written in pure Python to compute the basic composition of a protein sequence.

## Introduction

Understanding the compositional details and physicochemical properties of a protein is crucial for any protein-related study. These properties help determine the type (hydrophobic or hydrophilic), activity, and behavior of a protein. `proteinAnalysis2` provides an easy way to compute various compositional and physicochemical properties of a protein sequence.

### **Development & Publication**  
The development of this package was part of the **Anti-Dengue Peptide Database (ADPDB)** project, which has been published in the **International Journal of Peptide Research and Therapeutics** by **Springer Netherlands**. Click on the DOI in the citation section below to see the article.

## Requirements  

- Python 3.10 or higher  
- Pre-installed `peptides` and `Bio` packages  

## Updates  

The `proteinAnalysis2` package includes the following improvements over version 1:  

1. **Amino Acid Composition**: Computes amino acid counts and frequencies (previously unavailable).  
2. **Molecular Formula & Hydrophilicity**: Calculates the molecular formula and hydrophilicity of a protein sequence.  
3. **Single Function for All Properties**: Retrieves all compositional details and physicochemical properties in a dictionary format without calling multiple functions.  
4. **No Sequence Length Limitation**: Works with protein sequences of any length.  
5. **Optimized Performance**: Internally utilizes functions from `peptides` and `Bio` packages.  

## Computable Properties  

`proteinAnalysis2` can analyze the following:  

1. Length  
2. Molecular formula  
3. Amino acid counts  
4. Amino acid frequencies  
5. Missing amino acids  
6. Most occurring amino acids  
7. Least occurring amino acids  
8. Hydrophobic amino acid count  
9. Hydrophilic amino acid count  
10. Basic amino acid count  
11. Acidic amino acid count  
12. Modified amino acid count  
13. Modified amino acid frequencies  
14. Molecular weight  
15. Aliphatic index  
16. Instability index  
17. GRAVY (Grand Average of Hydropathicity)  
18. Hydrophilicity  
19. Hydrophobic moment  
20. Net charge  
21. Isoelectric point  
22. Structural class  
23. Mass shift  
24. Aromaticity  
25. Secondary structure fraction  
26. Molar extinction coefficient (Cysteine | Cysteine)  

## Download the Package  

The package can be downloaded by right-clicking and selecting **'Save link as...'**  

## Usage Example  

```python
# Import the package
from proteinAnalysis2 import proteinAnalysis2 as pan2

# Provide a protein sequence
seq = "ZAQALLPIASYAGLTVSAPVFAALVTVYGAYALYRYNIRRRENSYQRIRSDHDSHSCANNRGWCRPTCFSHEYTDWFNNDVCGSYRCCRPGRR"

# Create an object
obj = pan2.proteinAnalysis2(seq)

# Compute all compositional and physicochemical properties in dictionary format
print(obj.all_comp_physP())

# Get only the headers in list format
print(obj.header_info())
```

### **Citation**  
**Mondal, R.K., Anand, A.A. & Samanta, S.K., 2024.** ADPDB: A Comprehensive Knowledgebase of Manually Curated Peptides Against Dengue Virus. *Int J Pept Res Ther 30, 36 (2024).*  
**DOI:** [https://doi.org/10.1007/s10989-024-10610-0](https://doi.org/10.1007/s10989-024-10610-0)

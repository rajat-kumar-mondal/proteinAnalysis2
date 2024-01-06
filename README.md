# proteinAnalysis2
`proteinAnalysis2` is an updated version of `proteinAnalysis (version 1)` python package developed at Biochemistry & Bioinformatics Lab, Department of Applied Sciences, Indian Institute of Information Technology Allahabad (IIIT-A), Devghat, Jhalwa, Prayagraj-211015, U. P. India.

This package is written in pure python language to compute very basic composition of a protein sequence.

## Introduction
It is very important to know the compositional details and physicochemical properties for any kind of protein related study. The compositional details and physicochemical properties of a protein plays a crucial role to decide the type (hydrophobic or hydrophilic), activity and behaviour of a protein. This package can be use to find various compositional details and physicochemical properties of a protein.

## Requirements
python 3.10 or higher versions required. The package also requires pre-installied `peptides` and `Bio` package.

## Updates
proteinAnalysis2 package has the following updates:

   1. This version of package can compute the amino acid counts and frequencies where previous can't.
   2. This version of the package also able to calculate molecular formula and hydrophilicity of a protein sequence.
   3. The current version of the package also allow a user to get all the compositional details and physicochemical properties of protein sequence in dictionary by just calling a single function. So, their is no need to call multiple funstions.
   4. There is no limitation on size of the protein sequence.
   5. This packge internally use some functions of `peptides` and `Bio` package.

Following things can be computed by using  proteinAnalysis2

   1. Length
   2. Molecular formula
   3. Amino Acid counts
   4. Amino Acid frequencies
   5. Missing Amino Acid(s)
   6. Most Occurring Amino Acid(s)
   7. Less Occurring Amino Acid(s)
   8. Hydrophobic Amino Acid(s) Count
   9. Hydrophilic Amino Acid(s) Count
   10. Basic Amino Acid(s) Count
   11. Acidic Amino Acid(s) Count
   12. Modified Amino Acid(s) Count
   13. Modified Amino Acid(s) Frequencies
   14. Molecular weight
   15. Aliphatic index
   16. Instability index
   17. GRAVY
   18. Hydrophilicity
   19. Hydrophobic moment
   20. Net charge
   21. Isoelectric point
   22. Structural class
   23. Mass shift
   24. Aromaticity
   25. Secondary structure fraction
   26. Molar extinction coefficient (cysteine|cysteine)'

## Download the package
The package can be download by right click and then 'Save link as...'

## Usage example
```python
# import the package
import proteinAnalysis2 as pan2

# take a protein sequence
seq = "ZAQALLPIASYAGLTVSAPVFAALVTVYGAYALYRYNIRRRENSYQRIRSDHDSHSCANNRGWCRPTCFSHEYTDWFNNDVCGSYRCCRPGRR"

# create an object
obj = pan2.proteinAnalysis2(seq)

# computing compositional information and physicochemical properties in dictionary format
print(obj.all_comp_physP())

# know the headers only in list format
print(obj.header_info())
```

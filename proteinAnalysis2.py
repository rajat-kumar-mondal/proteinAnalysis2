from Bio.SeqUtils.ProtParam import ProteinAnalysis
import peptides
from collections import Counter
class proteinAnalysis:
    def __init__(self, seq) -> None:
        self.seq = seq
        resList = list('ARNDCQEGHILKMFPSTWYV'+'OUBZJX')
        resList.sort()
        resDict = {}
        for i in resList:
            resDict.update({i: self.seq.upper().count(i)})
        maxOcc = max([resDict[x] for x in resDict])
        minOcc = min([resDict[x] for x in resDict if resDict[x] > 0])
        missingAA = []
        self.maxOccAA = []
        self.minOccAA = []
        for x in resDict:
            if resDict[x] == 0:
                missingAA.append(x)
            if resDict[x] == maxOcc:
                self.maxOccAA.append(x)
            if resDict[x] == minOcc:
                self.minOccAA.append(x)
        self.missAA = list("".join(missingAA).replace('B','').replace('J','').replace('O','').replace('U','').replace('X','').replace('Z',''))
        self.phoCount = 0
        for i in list('GAMLIVFWP'):
            for j in self.seq.upper():
                if i == j:
                    self.phoCount += 1
        self.phiCount = 0
        for i in list('RNDCQEHKSTY'):
            for j in self.seq.upper():
                if i == j:
                    self.phiCount += 1
        self.basicCount = 0
        for k in list('HRK'):
            for j in self.seq.upper():
                if k == j:
                    self.basicCount += 1
        self.acidicCount = 0
        for k in list('DE'):
            for j in self.seq.upper():
                if k == j:
                    self.acidicCount += 1
        self.modAAdict = {}
        self.modAAdictFreq = {}
        for i in list('OUBZJX'):
            couRec = 0
            for j in self.seq.upper():
                if i == j:
                    couRec += 1
            if couRec > 0:
                self.modAAdict.update({i: couRec})
                self.modAAdictFreq.update({i: couRec/len(self.seq)})
    def missingResidues(self) -> str: return f'{"; ".join(self.missAA)}'
    def mostOccuringResidues(self) -> str: return f'{"; ".join(self.maxOccAA)}'
    def lessOccuringResidues(self) -> str: return f'{"; ".join(self.minOccAA)}'
    def hydrophobicAACount(self) -> int: return self.phoCount
    def hydrophilicAACount(self) -> int: return self.phiCount
    def basicAACount(self) -> int: return self.basicCount
    def acidicAACount(self) -> int: return self.acidicCount
    def modifiedAACount(self) :
        if len(self.modAAdict) > 0: return self.modAAdict
        else : return 0
    def modifiedAAFrequency(self):
        if len(self.modAAdict) > 0: return self.modAAdictFreq
        else : return 0
def mol_formula(protein_sequence):
    amino_acid_composition = {
        'A': {'C': 3, 'H': 7, 'N': 1, 'O': 2, 'S': 0},
        'C': {'C': 3, 'H': 7, 'N': 1, 'O': 2, 'S': 1},
        'D': {'C': 4, 'H': 7, 'N': 1, 'O': 4, 'S': 0},
        'E': {'C': 5, 'H': 9, 'N': 1, 'O': 4, 'S': 0},
        'F': {'C': 9, 'H': 11, 'N': 1, 'O': 2, 'S': 0},
        'G': {'C': 2, 'H': 5, 'N': 1, 'O': 2, 'S': 0},
        'H': {'C': 6, 'H': 9, 'N': 3, 'O': 2, 'S': 0},
        'I': {'C': 6, 'H': 13, 'N': 1, 'O': 2, 'S': 0},
        'K': {'C': 6, 'H': 14, 'N': 2, 'O': 2, 'S': 0},
        'L': {'C': 6, 'H': 13, 'N': 1, 'O': 2, 'S': 0},
        'M': {'C': 5, 'H': 11, 'N': 1, 'O': 2, 'S': 1},
        'N': {'C': 4, 'H': 8, 'N': 2, 'O': 3, 'S': 0},
        'P': {'C': 5, 'H': 9, 'N': 1, 'O': 2, 'S': 0},
        'Q': {'C': 5, 'H': 10, 'N': 2, 'O': 3, 'S': 0},
        'R': {'C': 6, 'H': 14, 'N': 4, 'O': 2, 'S': 0},
        'S': {'C': 3, 'H': 7, 'N': 1, 'O': 3, 'S': 0},
        'T': {'C': 4, 'H': 9, 'N': 1, 'O': 3, 'S': 0},
        'V': {'C': 5, 'H': 11, 'N': 1, 'O': 2, 'S': 0},
        'W': {'C': 11, 'H': 12, 'N': 2, 'O': 2, 'S': 0},
        'Y': {'C': 9, 'H': 11, 'N': 1, 'O': 3, 'S': 0},
        'B': {'C': 4, 'H': 7, 'N': 1, 'O': 3, 'S': 0},
        'J': {'C': 6, 'H': 13, 'N': 1, 'O': 2, 'S': 0},
        'O': {'C': 5, 'H': 10, 'N': 2, 'O': 2, 'S': 0},
        'U': {'C': 3, 'H': 7, 'N': 1, 'O': 2, 'S': 1},
        'Z': {'C': 5, 'H': 10, 'N': 2, 'O': 3, 'S': 0}
    }
    composition = Counter()
    for amino_acid in protein_sequence:
        composition += Counter(amino_acid_composition.get(amino_acid, {}))
    molecular_formula = ''.join(f"{element}{count}" for element, count in composition.items())
    return molecular_formula
def calculate_hydrophilicity(protein_sequence):
    hydrophilicity_values = {
        'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
        'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
        'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
        'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
    }    
    protein_sequence = protein_sequence.upper()
    total_hydrophilicity = 0.0    
    for amino_acid in protein_sequence:
        if amino_acid in hydrophilicity_values:
            total_hydrophilicity += hydrophilicity_values[amino_acid]    
    return total_hydrophilicity
def extract_alphabetical_characters(seq_string):
    return ''.join(char for char in seq_string if char.isalpha())
class proteinAnalysis2:
    def __init__(self, myseq) -> None:
        self.myseq = myseq
    def header_info(self):
        return "Sequence\tLength\tMolecular formula\tAA counts\tAA frequencies\tMissing AA\tMost occurring AA\tLess occurring AA\tHydrophobic AA counts\tHydrophilic AA counts\tBasic AA counts\tAcidic AA counts\tModified AA counts\tModified AA frequencies\tMolecular weight\tAliphatic index\tInstability index\tGRAVY\tHydrophilicity\tHydrophobic moment\tNet charge\tIsoelectric point\tStructural class\tMass shift\tAromaticity\tSecondary structure fraction\tMolar extinction coefficient (cysteine|cysteine)".split("\t")
    def all_comp_physP(self):
        seq = self.myseq
        protein_sequence = extract_alphabetical_characters(seq.upper()).replace(",", "").replace("-", "")
        protein = ProteinAnalysis(protein_sequence)
        panObj = proteinAnalysis(protein_sequence)
        peptObj = peptides.Peptide(protein_sequence)
        protein2 = ProteinAnalysis(protein_sequence.replace("B", "").replace("J", "").replace("X", "").replace("U", "").replace("Z", "").replace("O", ""))
        out = ''
        x = ''
        try:
            x = protein2.get_amino_acids_percent()
        except:
            x = "0 [Only special amino acid(s) is/are present in this entry. Please see Special/Modified Amino Acids Frequencies Section (below in the same table)]"
        out+=f"{protein_sequence}\t{len(protein_sequence)}\t{mol_formula(protein_sequence)}\t{protein2.count_amino_acids()}\t{x}\t{panObj.missingResidues()}\t{panObj.mostOccuringResidues()}\t{panObj.lessOccuringResidues()}\t{panObj.hydrophobicAACount()}\t{panObj.hydrophilicAACount()}\t{panObj.basicAACount()}\t{panObj.acidicAACount()}\t{panObj.modifiedAACount()}\t{panObj.modifiedAAFrequency()}\t"
        out+=f"{peptObj.molecular_weight()}\t{peptObj.aliphatic_index()}\t{peptObj.instability_index()}\t{peptObj.hydrophobicity()}\t"
        if len(set(protein_sequence).intersection(set('BJOUXZ'))) == 0:
            out+=f"{calculate_hydrophilicity(protein_sequence)}\t"
        else:
            out+=f"{calculate_hydrophilicity(protein_sequence)} [Note: Unable to consider special amino acids for now]\t"    
        out+=f"{peptObj.hydrophobic_moment()}\t{peptObj.charge()}\t{peptObj.isoelectric_point()}\t{peptObj.structural_class()}\t{peptObj.mass_shift()}\t"
        out+=f"{protein.aromaticity()}\t{protein.secondary_structure_fraction()}\t{protein.molar_extinction_coefficient()}\n"
        _info = out.rstrip().split("\t")
        _header = "Sequence\tLength\tMolecular formula\tAA counts\tAA frequencies\tMissing AA\tMost occurring AA\tLess occurring AA\tHydrophobic AA counts\tHydrophilic AA counts\tBasic AA counts\tAcidic AA counts\tModified AA counts\tModified AA frequencies\tMolecular weight\tAliphatic index\tInstability index\tGRAVY\tHydrophilicity\tHydrophobic moment\tNet charge\tIsoelectric point\tStructural class\tMass shift\tAromaticity\tSecondary structure fraction\tMolar extinction coefficient (cysteine|cysteine)".split("\t")
        _infoDict = {}
        for x, y in zip(_header, _info): _infoDict.update({x: y})
        return _infoDict


if __name__ == "__main__":
    seq = "ZAQALLPIASYAGLTVSAPVFAALVTVYGAYALYRYNIRRRENSYQRIRSDHDSHSCANNRGWCRPTCFSHEYTDWFNNDVCGSYRCCRPGRR"
    pan2obj = proteinAnalysis2(seq)
    # for computing compositional information and physicochemical properties in dictionary format
    print(pan2obj.all_comp_physP())
    # to know the headers only in list format
    print(pan2obj.header_info())
import pdb


class InferredAllele:
    def __init__(self):
        self.inferred_seq = {}
        self.last_allele_index = {}

    def get_inferred_allele(self, sequence: str, allele: str) -> str:
        """Infer allele from the sequence

        Args:
            sequence (str): sequence to infer the allele

        Returns:
            str: inferred allele
        """
        if sequence not in self.inferred_seq:
            return self.set_inferred_allele(sequence, allele)
        return self.inferred_seq[sequence]

    def set_inferred_allele(self, sequence: str, allele: str) -> None:
        """Set the inferred allele for the sequence

        Args:
            sequence (str): sequence to infer the allele
            allele (str): inferred allele
        """
        if allele not in self.last_allele_index:
            self.last_allele_index[allele] = 0
        self.last_allele_index[allele] += 1
        self.inferred_seq[sequence] = self.last_allele_index[allele]
        return self.inferred_seq[sequence]

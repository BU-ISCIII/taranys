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
        inf_value = self.last_allele_index.get(allele, 0) + 1
        self.inferred_seq[sequence] = inf_value
        return self.inferred_seq[sequence]

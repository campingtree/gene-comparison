from itertools import product
from typing import Dict, List

from Bio.Seq import Seq
from config import config
from utils import IO

from core.sequencer import Sequencer


class CodonFrequencyCalculator:
    def __init__(self):
        self.all_codons = self.__get_all_codons()
        self.all_dicodons = self.__get_all_dicodons()

    def get_codon_frequencies(self, seq: Seq, indices: List[List[int]]) -> Dict[str, float]:
        freq_table = dict.fromkeys(self.all_codons, 0)
        seq_str = str(seq)

        IO.info(f'Calculating codon frequencies in {len(indices)} fragments...')

        for fragment_start_id, fragment_end_id in indices:
            fragment = seq_str[fragment_start_id:fragment_end_id]
            for codon_start_id, codon_end_id in Sequencer.get_frame_codons_iter(Seq(fragment)):
                fragment_codon = fragment[codon_start_id:codon_end_id]
                freq_table[fragment_codon] += 1

        total_codon_count = sum(freq_table.values())
        return self.__normalize_frequency_table(freq_table, total_codon_count)

    def get_dicodon_frequencies(self, seq: Seq, indices: List[List[int]]) -> Dict[str, float]:
        freq_table = dict.fromkeys(self.all_dicodons, 0)
        seq_str = str(seq)

        IO.info(f'Calculating dicodon frequencies in {len(indices)} fragments...')

        for fragment_start_id, fragment_end_id in indices:
            fragment = seq_str[fragment_start_id:fragment_end_id]
            for dicodon_start_id, dicodon_end_id in Sequencer.get_frame_dicodons_iter(Seq(fragment)):
                fragment_dicodon = fragment[dicodon_start_id:dicodon_end_id]
                freq_table[fragment_dicodon] += 1
        
        total_dicodon_count = sum(freq_table.values())
        return self.__normalize_frequency_table(freq_table, total_dicodon_count)

    @staticmethod
    def __get_all_codons() -> List[str]:
        codons = [''.join(x) for x in product(config.NUCLEOTIDES, repeat=3)]
        return codons

    @staticmethod
    def __get_all_dicodons() -> List[str]:
        dicodons = [''.join(x) for x in product(config.NUCLEOTIDES, repeat=6)]
        return dicodons

    @staticmethod
    def __normalize_frequency_table(freq_table: Dict[str, int], element_count: int) -> Dict[str, float]:
        for codon in freq_table:
            freq_table[codon] /= element_count

        return freq_table

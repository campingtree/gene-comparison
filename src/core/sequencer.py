from typing import Dict, Iterator, List, Tuple
from Bio.SeqIO.FastaIO import FastaIterator
from Bio.Seq import Seq
from os.path import basename

from utils import IO
from config import config


class Sequencer:
    def __init__(self):
        self.sequences = {}

    # NOTE: galima butu iskelti file atidaryma i main ir cia priimti seq kaip argumenta. Bet jei sito metodo niekur daugiau nereikes kviesti, tai don't bother..
    def find_start_stop_fragments(self, path: str) -> Tuple[Seq, List[List[int]]]:
        start_codons = config.START_CODONS
        stop_codons = config.STOP_CODONS
        base_path_name = basename(path)

        with open(path) as handle:
            for record in FastaIterator(handle):
                nucleotide_sequence = record.seq
                nucleotide_sequence.name = base_path_name.rsplit('.', 1)[0]
                self.sequences[nucleotide_sequence] = []

                IO.info(f'Parsing START-STOP codon fragments in {base_path_name}...')
                for strand in nucleotide_sequence, nucleotide_sequence.reverse_complement():
                    for frame in range(3):
                        strand_start_ids = []
                        for codon_start_id, codon_end_id in self.get_frame_codons_iter(strand, frame):
                            if str(strand[codon_start_id:codon_end_id]) in start_codons:
                                strand_start_ids.append(codon_start_id)
                            if str(strand[codon_start_id:codon_end_id]) in stop_codons:
                                for strand_start_id in strand_start_ids:
                                    self.sequences[nucleotide_sequence].append([strand_start_id, codon_end_id]) # TODO: can prob use tuples here instead, tik tada pakeisti ir return value indikatoriu
                                strand_start_ids.clear()
        return nucleotide_sequence, self.sequences[nucleotide_sequence]
    

    def find_longest_codons_for_common_stop(self, seq: Seq) -> Dict[int, int]:
        if not self.sequences[seq]:
            IO.error('seq not found in sequences. Call find_start_stop_fragments first.')
            return

        unique_stops = {} # stop:largest_length
        for codon_start_id, codon_end_id in self.sequences[seq]:
            if codon_end_id not in unique_stops:
                unique_stops[codon_end_id] = -1
            length = codon_end_id - codon_start_id
            assert length % 3 == 0
            if length > unique_stops[codon_end_id]:
                unique_stops[codon_end_id] = length

        # return { stop-length: stop for stop, length in unique_stops.items() }
        self.sequences[seq] = [[stop-length, stop] for stop, length in unique_stops.items()]
        return self.sequences[seq]

    def filter_out_fragments_shorter_than(self, seq: Seq, n: int) -> None:
        if not self.sequences[seq]:
            IO.error('seq not found in sequences. Call find_start_stop_fragments first.')
            return

        self.sequences[seq] = list(filter(lambda x: x[1] - x[0] >= n, self.sequences[seq]))
        return self.sequences[seq]

    @staticmethod
    def get_frame_codons_iter(seq: Seq, frame=0) -> Iterator[int]:
        seq_length = len(seq)
        length = 3
        for index in range(frame, seq_length, length):
            if index + length > seq_length:
                break

            chunk = index, index+length
            yield chunk

    @staticmethod
    def get_frame_dicodons_iter(seq: Seq, frame=0) -> Iterator[int]:
        # TODO: if this actually works, add explanation why length is 3 here and why length*2 (3 kodonai - 2 dikodonai)
        seq_length = len(seq)
        length = 3
        for index in range(frame, seq_length, length):
            if index + (length+length) > seq_length:
                break

            chunk = index, index+(length+length)
            yield chunk

    # TODO: delete
    @staticmethod
    def __chunk_sequece_iter(seq: Seq, start: int, length: int) -> Iterator[int]:
        seq_length = len(seq)
        for index in range(start, seq_length, length):
            if index + length > seq_length:
                break

            chunk = index, index+length
            yield chunk

from config import config
from core import CodonFrequencyCalculator, PhylipDistanceMatrix, Sequencer
from utils import IO

if __name__ == '__main__':
    frequencyCalculator = CodonFrequencyCalculator()

    codon_frequency_tables = {}
    dicodon_frequency_tables = {}

    sequencer = Sequencer()

    for data_file_path in IO.get_data_paths():
        # 1
        seq, indices = sequencer.find_start_stop_fragments(data_file_path)
        # 2
        longest_fragment_indices = sequencer.find_longest_codons_for_common_stop(seq)
        # 3
        filtered_fragments_indices = sequencer.filter_out_fragments_shorter_than(seq, config.MIN_FRAGMENT_LENGTH)
        # 4 
        codon_frequency_tables[seq] = frequencyCalculator.get_codon_frequencies(seq, filtered_fragments_indices)
        dicodon_frequency_tables[seq] = frequencyCalculator.get_dicodon_frequencies(seq, filtered_fragments_indices)
        print()

    # 5
    IO.info('Calculating Phylip distance matrix for codon frequencies...')
    codon_distance_matrix = PhylipDistanceMatrix(codon_frequency_tables)
    codon_distance_matrix.calculate()
    codon_distance_matrix.print()

    print()

    # 5
    IO.info('Calculating Phylip distance matrix for dicodon frequencies...')
    dicodon_distance_matrix = PhylipDistanceMatrix(dicodon_frequency_tables)
    dicodon_distance_matrix.calculate()
    dicodon_distance_matrix.print()

    print()
    IO.info('FINISHED')

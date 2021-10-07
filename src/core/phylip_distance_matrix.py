import math

from utils import IO


class PhylipDistanceMatrix:
    def __init__(self, freq_tables):
        self.__freq_tables = freq_tables
        self.__matrix = None

    def calculate(self) -> None:
        rows = []
        for first_freq_table_seq in self.__freq_tables:
            row = [first_freq_table_seq.name]
            rows.append(row)
            for second_freq_table_seq in self.__freq_tables:
                s = 0
                freq_items_count = len(self.__freq_tables[first_freq_table_seq])
                assert freq_items_count == 64 or freq_items_count == 4096
                for freq_item_key in self.__freq_tables[first_freq_table_seq]:
                    first_table_frequencies = self.__freq_tables[first_freq_table_seq];
                    second_table_frequencies = self.__freq_tables[second_freq_table_seq];
                    s += self.__MSE(first_table_frequencies[freq_item_key], second_table_frequencies[freq_item_key])
                distance = math.sqrt(s / freq_items_count)
                row.append(f'{format(distance, f".12f")} ')
        self.__matrix = rows

    def print(self) -> None:
        if not self.__matrix:
            IO.error('no calculated matrix found. Call calculate() first.')
            return

        print(self.__get_phylip_header())
        for row in self.__matrix:
            for column in row:
                print(f'{column} ', end='')
            print()

    def __get_phylip_header(self) -> str:
        return f'{len(self.__freq_tables)}'

    @staticmethod
    def __MSE(x: float, y: float) -> float:
        return (x - y) ** 2

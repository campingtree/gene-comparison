from pathlib import Path
from typing import List
from glob import glob
from datetime import datetime

from config import config

class IO:
    __white_escape = '\u001b[37m'
    __red_escape = '\u001b[31m'
    __yellow_escape = '\u001b[33m'
    __reset_escape = '\u001b[0m'

    @staticmethod
    def get_project_root() -> str:
        return Path(__file__).parent.parent.parent

    @staticmethod
    def get_data_paths(extension = 'fasta') -> List[str]:
        data_path = f'{IO.get_project_root()}/data'

        return sorted(glob(f'{data_path}/*.{extension}'))

    @staticmethod
    def __string_with_color(message: str, symbol: str, color: str, with_timestamp=False) -> str:
        timestamp = f'[{datetime.now().strftime(config.TIMESTAMP_FORMAT)}] ' if with_timestamp else ''
        return f'{color}{timestamp}[{symbol}]{IO.__reset_escape} {message}'

    @staticmethod
    def info(message: str, with_timestamp=True) -> None:
        print(IO.__string_with_color(message, '*', IO.__white_escape, with_timestamp))

    @staticmethod
    def warn(message: str, with_timestamp=True) -> None:
        print(IO.__string_with_color(message, '?', IO.__yellow_escape, with_timestamp))

    @staticmethod
    def error(message: str, with_timestamp=True) -> None:
        print(IO.__string_with_color(message, '!', IO.__red_escape, with_timestamp))
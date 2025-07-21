"""
Module for managing the IS database from ISFinder.
"""
from csv import DictReader

from eris import RESOURCES
from eris.io import SeqFile
from eris.seq import Qualifier
from eris.utils import download, is_non_empty_file

# Constants ------------------------------------------------------------------------------------------------------------
_CSV_NAME = 'IS.csv'
_FNA_NAME = 'IS.fna'
_BASE_URL = 'https://raw.githubusercontent.com/thanhleviet/ISfinder-sequences/refs/heads/master/'
_CSV_URL = f'{_BASE_URL}{_CSV_NAME}'
_FNA_URL = f'{_BASE_URL}{_FNA_NAME}'
_CSV_PATH = RESOURCES.data / _CSV_NAME
_FNA_PATH = RESOURCES.data / _FNA_NAME
# _HMMS = ['PF01548', 'PF01609', 'PF01797']
# _CMS = ['RF03014', 'RF03015']


# Classes --------------------------------------------------------------------------------------------------------------
class Database:
    """
    A class to manage the ISFinder database, including downloading and parsing the data.

    Attributes:
        records: dict[str, 'Record'] A dictionary of the ISFinder sequences in the database.
        qualifiers: list[str] A list of metadata qualifiers available for each record.
    """
    def __init__(self):

        if not is_non_empty_file(_CSV_PATH):
            download(_CSV_URL, _CSV_PATH)

        if not is_non_empty_file(_FNA_PATH):
            download(_FNA_URL, _FNA_PATH)

        self.records: dict[str, 'Record'] = {}

        with open(_CSV_PATH, 'rt') as csv_file, SeqFile(_FNA_PATH, 'fasta') as fna_file:
            reader = DictReader(csv_file)
            self.qualifiers: list[str] = reader.fieldnames
            metadata = {i['Name']: i for i in reader}
            for record in fna_file:  # type: Record
                if record_metadata := metadata.get(record.id.partition('_')[0]):
                    record.qualifiers += [Qualifier(k, v) for k, v in record_metadata.items()]
                    self.records[record.id] = record
                # else:
                #     warn(f'{record} does not have metadata', ErisWarning)

    def __repr__(self):
        return f'{self.__class__.__name__}({len(self.records)} records)'

    def __len__(self):
        return len(self.records)

    def __getitem__(self, item):
        return self.records[item]

    def __iter__(self):
        return iter(self.records.values())


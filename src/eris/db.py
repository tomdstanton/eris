"""
Copyright 2025 Tom Stanton (tomdstanton@gmail.com)
https://github.com/tomdstanton/eris

This file is part of eris. eris is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. eris is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with eris.
If not, see <https://www.gnu.org/licenses/>.
"""
from warnings import warn
from csv import DictReader

from eris import RESOURCES, ErisWarning
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
    def __init__(self):

        if not is_non_empty_file(_CSV_PATH):
            download(_CSV_URL, _CSV_PATH)

        if not is_non_empty_file(_FNA_PATH):
            download(_FNA_URL, _FNA_PATH)

        self.records = {}

        with open(_CSV_PATH, 'rt') as csv_file, SeqFile(_FNA_PATH, 'fasta') as fna_file:
            reader = DictReader(csv_file)
            self.qualifiers = reader.fieldnames
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


"""
Module for managing the IS database from ISFinder.
"""
from csv import DictReader
# from re import compile, MULTILINE
# import ssl  # NOT recommended but only way to access ISFinder database
# ssl._create_default_https_context = ssl._create_unverified_context

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
# _HMMS = ['PF01548', 'PF01609', 'PF01797']  # If we ever need hmmer
# _CMS = ['RF03014', 'RF03015']  # If we ever need infernal
# IS Finder URLs -------------------------------------------------------------------------------------------------------
# _ISFINDER_BASE_URL = 'https://www-is.biotoul.fr/scripts/'
# _ISFINDER_SEARCH_URL = 'search-db.php'
# # Regex used to match desired groups (link names, IS name, IS sequence)
# _LINK_REGEX = compile(r'(ficheIS.php\?name=[^\']+)')
# _NUCL_REGEX = compile(
#     r"<div><p class='entete_propriete'>DNA sequence </p>.*?\n*?<div class='seq'>(?P<seq>(?:.|\n)+?)</div>", MULTILINE)
# _PROT_REGEX = compile(
#     r"class='entete_propriete'>ORF sequence : </p>.*?\n*?<div class='seq'>(?P<seq>(?:.|\n)+?)</div>", MULTILINE)
# _STRIP_TAG_REGEX = compile(r'(<.+?>\s?)')
# _ORGANISM_REGEX = compile(r'<td><div class="ascenseurAuto">(.*?)</div></td>')
# _FAMILY_REGEX = compile(r"class='entete_propriete'>Family </span>(.*?)</li>")
# _GROUP_REGEX = compile(r"class='entete_propriete'>Group </span>(.*?)</li>")
# _LENGTH_REGEX = compile(r"class='entete_propriete'>IS Length : </span>(.*?)bp</div>")


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

# Functions ------------------------------------------------------------------------------------------------------------
# def _fetch_data_from_link(link: str):
#     if html := download(link):
#         nucl_seq = _STRIP_TAG_REGEX.sub('', n.group('seq')).strip().upper() if (n := _NUCL_REGEX.search(html)) else ''
#         prot_seq = _STRIP_TAG_REGEX.sub('', n.group('seq')).strip().upper() if (n := _PROT_REGEX.search(html)) else ''
#         group = n.group(1).strip() if (n := _GROUP_REGEX.search(html)) else ''
#         family = n.group(1).strip() if (n := _FAMILY_REGEX.search(html)) else ''
#         length = n.group(1).strip() if (n := _LENGTH_REGEX.search(html)) else ''
#         organism = n.group(1).strip() if (n := _ORGANISM_REGEX.search(html)) else ''
#
#
# def _get_links_from_search(query):
#     """
#     Process search query data and prepare for a POST method request
#     """
#     data = {'name': query, 'namecond': 'contains', 'MGEtype': 'all', 'Onsubmit': 'Submit'}
#     result = download(f'{_ISFINDER_BASE_URL}{_ISFINDER_SEARCH_URL}', data=data, encode_data=True).decode()
#     return _LINK_REGEX.findall(result)  # Find all links in HTML source file
#

***********
Databases
***********

# eris db

The `db` module contains code for parsing, storing and managing eris databases, which represent sets of genes (loci)
each encoding a unique structure (serotype) of a particular antigen (surface-polysaccharide).

One or more databases are curated for a specific organism by a lab group or individual.

Unlike previous versions of eris (v<=3), the `db` module aims to decentralise the eris databases,
shifting the power and responsibility of database management to the curators.

The idea is to use [submodules](https://git-scm.com/book/en/v2/Git-Tools-Submodules) in the `databases` directory
to host databases from different repositories representing a specific organism; such that in the future, you only
need to install databases for
organisms of interest. Database information is kept in `scheme.json` files, and changes are tracked with
`CHANGELOG.md` files and Github releases to generate a version number.

When individual databases are updated, erisDB will pull the latest release and validate any changes, including
passing tests in eris and locus "unique-ness". If the changes are valid, eris release a new version reflecting the
changes.

## Database requirements
- Databases **must** comply to the eris format stipulated [here](https://eris.readthedocs.io/en/latest/Databases.html#format).
- Databases **must** include a `{db}.json` metadata file per database.

The metadata file tracks important information about the database itself and can't currently be embedded into the
database Genbank files (this could change if we use pickles in the future).

For example, the metadata for the _Klebsiella pneumoniae_ Species Complex K-locus database looks like this:

```json
 {
    "organism":  "Klebsiella pneumoniae Species Complex",
    "antigen": "CPS",
    "status": "ok",
    "keywords": ["kpsc_k", "kp_k", "k_k"],
    "n_loci": 161,
    "n_extra": 0,
    "locus_prefix": "KL",
    "type_prefix": "K",
    "id_threshold": 87.5,
    "core_genes": ["galF", "cpsACP", "wza", "wzb", "wzc", "wzi", "wzy", "wzx", "ugd"],
    "version": "0.1.0",
    "last_updated": "21/10/2024",
    "institute": "Monash University",
    "doi": ["https://doi.org/10.1099/mgen.0.000102"],
    "contact": {"Thomas Stanton": "tom.stanton@monash.edu"},
    "url": "https://github.com/tomdstanton/KpSC_erisDB.git"
}
 ```

- Databases **must** include a `CHANGELOG.md` file containing a list of changes made to the database in question.


### The ``Database`` class
This class represents a eris database held in memory. If the first argument is an existing file or keyword
the database will be parsed and loaded. E.g. to load the _Klebsiella pneumoniae_ Species Complex K-locus database:

```python
from eris.db import Database
db = Database('kpsc_k')
```


# `eris.db.glycan`
###### A submodule for parsing, storing and manipulating glycan data

---

`glycan` is a `eris.db` submodule for linking bacterial surface polysaccharide antigens with genetic data.


### The `Glycan` class

The ``Glycan`` class is a submodule of ``Graph``, as glycan molecules can be represented as directed graphs of
``Monosaccharide``s (nodes) connected by ``GlycosidicLinkage``s (subclass of ``Edge``).

Representing glycans as graphs allows for easier comparison, e.g. the graph Smith-Waterman alignment algorithm
implemented in KCaM (KEGG Carbohydrate Matcher). We could also potentially use the graph edit distance method
to compare glycans. This is not currently implemented in ``eris.core.graph`` but it is in NetworkX.

Each ``Glycan`` instance acts as a parser, and accepts a linear glycan representation as input.

# `eris.db.annotate`
###### A submodule for annotating genomes using a eris database as a template

---

`annotate` is a `eris.db` submodule for annotating genomes using a eris database as a template.

### The `Annotator` class

# eris Predict

The `predict` module contains code for predicting antigen strucutres from locus genes using machine learning.

The initial idea is to represent antigen-specific (i.e. non-machinery) locus genes as a directed acyclic graph
(DAG) and the corresponding antigen structure also as a graph. These pairs can be used to train a graph
neural network (GNN), so when we come across a locus with an unknown antigen, we can predict the
structure using the model.

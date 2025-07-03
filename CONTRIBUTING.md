This, like other packages, is built on the `stantlib` framework, which makes heavy use of custom objects
to represent biological components in a hierarchical manner.

Tom is a fussy coder, aided by OCD and being a stubborn Sagittarius.

Tom likes as few dependencies as possible. A classic example is the `Biopython`
library. No one is arguing that `Biopython` is indispensable in bioinformatics,
however, it is a huge library with many modules which are redundant for a tool
performing a specific task. If you are wanting to translate DNA to protein, **write your own function**,
this is ~16 lines of code vs 3 million, and you can have fun making it as efficient as possible; don't
kid yourself, Tom knows he has.

Tom doesn't like classes of a certain hierarchical level referring to their 'parents'; e.g. the `Contig` class
shouldn't have an 'assembly' attribute referring to an `Assembly` object.

Sometimes this expensive to avoid; e.g. with reporting. If 

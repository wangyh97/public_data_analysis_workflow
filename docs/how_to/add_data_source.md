# Add a data source

Create `datasets/<source>/` for the source's retrieval and parsing rules.
Accept accession/dataset ID, raw root and output root on the command line.
Do not include cancer type or one GSE number in the adapter name. Emit raw
identity/checksums and documented metadata fields; hand off normalized
observations to `processing/`. Add a validator using the appropriate data
contract, register the source, and document how to distinguish source IDs
from original study IDs.

When a new accession has unusual metadata, first add a mapping record to
that study's run recipe. Change the adapter only when the file format or
source-wide rule changes. Existing study outputs must remain reproducible.

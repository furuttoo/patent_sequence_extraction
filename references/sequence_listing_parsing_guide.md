# Patent Sequence Listing — Parsing Guide

Read this before touching any patent's sequence listing file. It covers
format identification, the zip wrapper used by newer filings, and which
parser script to hand off to. This is a reference, not a runnable —
the actual parsing logic lives in `parse_sequence_listing_st25.py` and
`parse_sequence_listing_st26.py`.

## The three things you'll encounter, in order of filing date

| Era | Format | Extension | Wrapper |
|---|---|---|---|
| Pre-2022 | ST.25 (flat tags: `<210>`, `<400>`, ...) | `.txt` | none — plain text file |
| Post-2022 | ST.26 (XML: `<SequenceData>`, `<INSDSeq>`, ...) | `.xml` (sometimes `.XML`) | **often zipped** |
| Post-2022, zipped | ST.26 XML, same as above | `.zip` containing one `.xml`/`.XML` | zip archive |

Older ST.25 files come as a bare `.txt` and need no unwrapping. Newer
ST.26 filings frequently arrive zipped — WIPO's own filing tooling
produces a `.zip` per sequence listing rather than a loose XML file, so
treat "the sequence listing is a zip" as the normal case for anything
filed after ~2022, not an edge case.

## Step 1: identify what you actually have

Check extensions in this order — **don't assume the visible filename
tells you the format**, since a `.zip` could in principle wrap either
format (though in practice it's always ST.26 XML so far):

```python
from pathlib import Path

def classify_sequence_file(path):
    path = Path(path)
    ext = path.suffix.lower()
    if ext == ".zip":
        return "zipped"       # needs unzip step before format detection
    if ext == ".xml":
        return "st26"
    if ext == ".txt":
        return "st25"
    return "unknown"
```

## Step 2: if zipped, unzip first — then re-detect

**Never assume the extracted filename matches the patent number or the
zip's own name.** Real example from WO2023019246:

```
WO2023019246_sequence.zip
  └── 121301-16620.XML          <-- internal attorney docket number,
                                     completely unrelated to "WO2023019246",
                                     AND uppercase .XML extension
```

The zip's name follows the patent-number convention your pipeline
expects; the file inside it does not. It's typically named after the
applicant's own internal docket reference (here `121301-16620`, likely
Alnylam's matter number `ALN-...` equivalent under the hood), and the
extension case is inconsistent (`.XML` vs `.xml` across different filers'
export tooling). Handle both:

```python
import zipfile
from pathlib import Path

def unzip_sequence_listing(zip_path, extract_to):
    """Extract a sequence-listing zip and return the path to the sequence
    file inside it — found by extension, not by expecting a matching name."""
    zip_path = Path(zip_path)
    extract_to = Path(extract_to)
    extract_to.mkdir(parents=True, exist_ok=True)

    with zipfile.ZipFile(zip_path) as z:
        names = z.namelist()
        z.extractall(extract_to)

    # find the sequence file by extension, case-insensitively — do not
    # assume it shares any part of the zip's filename
    candidates = [n for n in names if n.lower().endswith((".xml", ".txt"))]
    if not candidates:
        raise ValueError(f"no .xml/.txt found inside {zip_path.name}: {names}")
    if len(candidates) > 1:
        # seen so far: always exactly one file per zip. If a zip ever
        # contains more than one, don't guess — surface it for a human
        # to confirm which is the actual sequence listing.
        raise ValueError(f"multiple candidate files in {zip_path.name}: {candidates}")

    return extract_to / candidates[0]
```

Extract to a working/scratch directory, not next to the original zip —
keeps the patent's source folder clean and makes re-runs idempotent
(delete the scratch dir and re-extract rather than worrying about stale
partial extractions).

## Step 3: route to the matching parser by the *extracted* file's extension

The two parser scripts live alongside this guide:

- **`parse_sequence_listing_st25.py`** — for `.txt` files (`<210>`/`<400>`
  flat-tag format). Exposes `parse_sequence_listing(path)` → `(summary,
  records)`, plus a `build_cache(path, force=False)` helper that writes
  `<stem>_seqcache.json` next to the source and returns
  `(cache_path, "cache_built" | "cache_hit")`.
- **`parse_sequence_listing_st26.py`** — for `.xml`/`.XML` files
  (`<SequenceData>`/`<INSDSeq>` schema). Same two functions, same
  return shapes — the two scripts are interchangeable at the call site.

```python
import importlib.util
from pathlib import Path

def load_parser(fmt):
    """fmt is 'st25' or 'st26', as returned by classify_sequence_file()."""
    script = {
        "st25": "parse_sequence_listing_st25.py",
        "st26": "parse_sequence_listing_st26.py",
    }[fmt]
    spec = importlib.util.spec_from_file_location(fmt, script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

extracted = unzip_sequence_listing(zip_path, scratch_dir)
fmt = classify_sequence_file(extracted)   # re-run classification — now "st25" or "st26"

parser = load_parser(fmt)
cache_path, status = parser.build_cache(extracted)
# or, without writing a cache file: summary, records = parser.parse_sequence_listing(extracted)
```

Equivalently, from the shell rather than importing — both scripts are
runnable directly and this is often the simpler path for a one-off patent:

```bash
python parse_sequence_listing_st26.py path/to/121301-16620.XML
python parse_sequence_listing_st25.py path/to/older_sequence.txt
```

Extension-based routing after extraction is sufficient — no zip
encountered so far has mixed formats, and both parser scripts already
validate their own assumptions (ST.25 via `<210>`/`<400>` cross-check,
ST.26 via XML schema parsing failing loudly on malformed input) so a
misrouted file will fail fast rather than silently producing garbage.

## Caching: key by the *patent's* identifier, not the inner filename

Both parser scripts write `<stem>_seqcache.json` next to the source file
they're given. If you hand them the extracted `121301-16620.XML`, the
cache lands there too, named after the docket number — which is fine
functionally, but disconnects the cache from the patent number your
pipeline actually indexes by. Two reasonable options:

1. **Rename on extraction** — after unzipping, rename the extracted file
   to `<patent_number>_sequence.xml` before handing it to the parser, so
   the cache filename stays consistent with every other patent folder.
2. **Keep the original name, track the mapping separately** — leave
   `121301-16620.XML` as-is and record `{patent_number: extracted_path}`
   in whatever index your pipeline already maintains for the patent
   folder.

Option 1 is simpler and matches the existing convention from the ST.25
side (`<patent>_sequence.txt` → `<patent>_sequence_seqcache.json`), so
default to it unless there's a reason to preserve the original filename
for provenance.

## Known quirks checklist (update this section as new ones surface)

- [x] Inner filename unrelated to outer zip/patent name (docket number,
      not patent number)
- [x] Extension case inconsistency: `.XML` vs `.xml`
- [ ] Zips containing more than one file — not yet seen; current code
      raises rather than guessing if it happens
- [ ] Nested zip-of-zips or zip-in-folder structures — not yet seen
- [ ] ST.25 `.txt` files arriving zipped — not yet seen (only ST.26 XML
      has been observed zipped so far), but the `classify_sequence_file`
      / `unzip_sequence_listing` split above already supports it if it
      does happen, since format detection runs *after* extraction

## Quick reference: which parser for which format

| Detected format | Script | Key functions |
|---|---|---|
| `.txt`, `<210>`/`<400>` tags | `parse_sequence_listing_st25.py` | `parse_sequence_listing(path)`, `build_cache(path)` |
| `.xml`/`.XML`, `<SequenceData>` | `parse_sequence_listing_st26.py` | `parse_sequence_listing(path)`, `build_cache(path)` |
| `.zip` | unzip per this guide, then re-classify the extracted file and call whichever of the above matches | — |

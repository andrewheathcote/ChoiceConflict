#!/bin/sh
# Build the supplement: doc/supplement.docx (Word, native equations) from doc/src/*.md, the generated tables and the figures.
#   sh doc/build.sh        (needs R, pandoc >= 3, poppler's pdftoppm and pdfinfo)
set -e
D=$(cd "$(dirname "$0")" && pwd); A=$(dirname "$D")
Rscript "$D/make_tables.R"
mkdir -p "$D/fig"
for f in "$A"/figures/*.pdf; do [ -e "$f" ] || continue; b=$(basename "$f" .pdf)
  # one PNG per page: <name>.png, or <name>-<page>.png for a multi-page PDF (convergence_grids: one page per task)
  if [ "$(pdfinfo "$f" | awk '/^Pages:/ {print $2}')" -gt 1 ]; then pdftoppm -r 150 -png "$f" "$D/fig/$b"; else pdftoppm -r 150 -png -singlefile "$f" "$D/fig/$b"; fi
done
python3 -I - "$D" <<'PY'
import sys, os, re
d = sys.argv[1]; src = os.path.join(d, "src")
parts = [open(os.path.join(src, f)).read() for f in ["front.md", "s1_s4.md", "s5.md", "s6_s10.md", "figures.md", "references.md"] if os.path.exists(os.path.join(src, f))]
s = "\n\n".join(parts)
s = re.sub(r"<!-- TABLE (\w+) -->", lambda m: open(os.path.join(src, "tables", m.group(1) + ".md")).read(), s)
open(os.path.join(d, "supplement.md"), "w").write(s)
PY
cd "$D"; pandoc supplement.md -o supplement.docx --resource-path="$D" ${REF:+--reference-doc=$REF}
echo "wrote $D/supplement.docx"

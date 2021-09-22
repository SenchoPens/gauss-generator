cat t.tex.head > t.tex
cat /dev/stdin | sed 's/\\leadsto/\\leadsto$$ $$\\leadsto/g' >> t.tex
cat t.tex.tail >> t.tex
# latexmk t.tex && zathura t.pdf

for filename in *.csv; do
        [ -f "$filename" ] || break

        # Trie le fichier selon l’AGI (en excluant le header) et l’enregistre dans un second fichier.
        { head -n 1 $filename && tail -n +2 $filename | sort -k4 -k6n,7n; } > ${filename}_sorted.csv

        # Crée un fichier indexé gz et crée l’index tabix.
        /biotools/htslib/1.9/bin/bgzip ${filename}_sorted.csv && /biotools/htslib/1.9/bin/tabix ${filename}_sorted.csv.gz -S 1 -s 4 -b 6 -e 7

done
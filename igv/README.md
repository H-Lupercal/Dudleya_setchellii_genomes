# IGV Visualization Bundle

These files are prepared for IGV Desktop. The FASTA headers and GFF3 seqids match, so the annotation tracks should appear on the loaded custom genome.

## Recommended Load

1. Open IGV.
2. Select `Genomes > Load Genome from File...`.
3. Load `igv/dudleya_organelles.igv.json`.
4. Use the chromosome dropdown to switch between `chloroplast` and `mitochondria`.

This loads both organelle references and the combined draft annotation track together.

## Chloroplast

1. Open IGV.
2. Select `Genomes > Load Genome from File...`.
3. Load `igv/chloroplast.igv.json`.
4. Use the location box for `chloroplast:1-176964`.

## Mitochondria

1. Select `Genomes > Load Genome from File...`.
2. Load `igv/mitochondria.igv.json`.
3. Use the location box for `mitochondria:1-243359`.

## Files

- `chloroplast.fa`: IGV reference FASTA for the chloroplast assembly.
- `chloroplast.fa.fai`: FASTA index for IGV.
- `chloroplast.gff3`: Best nonredundant draft chloroplast annotations.
- `chloroplast.igv.json`: IGV genome definition for the chloroplast only.
- `dudleya_organelles.fa`: Combined chloroplast and mitochondrial IGV reference FASTA.
- `dudleya_organelles.fa.fai`: Combined FASTA index for IGV.
- `dudleya_organelles.gff3`: Combined draft annotation track.
- `dudleya_organelles.igv.json`: Recommended IGV genome definition for both organelles.
- `mitochondria.fa`: IGV reference FASTA for the mitochondrial assembly.
- `mitochondria.fa.fai`: FASTA index for IGV.
- `mitochondria.gff3`: Best nonredundant draft mitochondrial annotations.
- `mitochondria.igv.json`: IGV genome definition for the mitochondrion only.

These annotations are draft homology-transfer annotations. Use them as evidence-guided feature calls, not as a fully curated organelle annotation.

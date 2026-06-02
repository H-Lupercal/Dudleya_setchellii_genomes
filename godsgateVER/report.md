# Dudleya setchellii candidate organelle genome verification and draft annotation

## Inputs
- chloroplast: Dudleya_hifiasm_purged_manual_chloroplast.fa; header `ptg000216l_1`; length 176,964 bp; GC 37.77%; N=0; ambiguous=0
- mitochondria: Dudleya_hifiasm_purged_manual_mitochondria.fa; header `ptg000317l_1`; length 243,359 bp; GC 45.517%; N=0; ambiguous=0

## Reference searches
- chloroplast: NCBI query `Dudleya[Organism] AND chloroplast[Title] AND complete genome` returned 9 records; fetched 9 records.
  - PX244394.1 (166,371 bp): Dudleya sp. T72600-1 chloroplast, complete genome
  - PX244393.1 (150,419 bp): Dudleya candelabrum chloroplast, complete genome
  - PX244392.1 (152,899 bp): Dudleya greenei chloroplast, complete genome
  - PX244391.1 (150,313 bp): Dudleya lanceolata chloroplast, complete genome
  - PX244390.1 (150,405 bp): Dudleya nesiotica chloroplast, complete genome
  - PX244389.1 (150,377 bp): Dudleya virens chloroplast, complete genome
  - NC_085682.1 (150,780 bp): Dudleya farinosa chloroplast, complete genome
  - OQ076651.1 (150,780 bp): Dudleya farinosa chloroplast, complete genome
  - OL312335.1 (330 bp): Dudleya farinosa voucher MO:Carlsen3301 photosystem II protein D1 (psbA) gene, partial cds; psbA-trnH intergenic spacer, complete sequence; and tRNA-His (trnH-GUG) gene, partial sequence; chloroplast
- mitochondria: NCBI query `Crassulaceae[Organism] AND mitochondrion[Title] AND complete genome` returned 16 records; fetched 16 records.
  - PV256627.1 (242,059 bp): Graptopetalum paraguayense mitochondrion, complete genome
  - PV608516.1 (156,728 bp): Sedum sarmentosum mitochondrion, complete genome
  - CM082116.1 (79,921 bp): Rhodiola kirilowii isolate R_Ki mitochondrion, complete sequence, whole genome shotgun sequence
  - PP024540.1 (259,150 bp): Rhodiola rosea mitochondrion, complete genome
  - NC_082108.1 (202,019 bp): Rhodiola juparensis mitochondrion, complete genome
  - OR188139.1 (202,019 bp): Rhodiola juparensis mitochondrion, complete genome
  - NC_072122.1 (257,378 bp): Rhodiola tangutica mitochondrion, complete genome
  - NC_070303.1 (194,106 bp): Rhodiola crenulata mitochondrion, complete genome
  - NC_069572.1 (212,159 bp): Sedum plumbizincicola voucher KPBK001 mitochondrion, complete genome
  - OP573219.1 (257,378 bp): Rhodiola tangutica mitochondrion, complete genome

## Identity evidence
### chloroplast_vs_chloroplast_refs
- gb|PX244389.1|: query coverage 100.00%, weighted identity 99.565%, reference coverage 99.99%, HSPs 7; Dudleya virens chloroplast, complete genome
- gb|PX244393.1|: query coverage 100.00%, weighted identity 99.388%, reference coverage 99.97%, HSPs 8; Dudleya candelabrum chloroplast, complete genome
- gb|PX244390.1|: query coverage 100.00%, weighted identity 99.368%, reference coverage 99.99%, HSPs 8; Dudleya nesiotica chloroplast, complete genome
- gb|OQ076651.1|: query coverage 99.95%, weighted identity 99.359%, reference coverage 99.83%, HSPs 9; Dudleya farinosa chloroplast, complete genome
- ref|NC_085682.1|: query coverage 99.95%, weighted identity 99.359%, reference coverage 99.83%, HSPs 9; Dudleya farinosa chloroplast, complete genome
- gb|PX244391.1|: query coverage 99.95%, weighted identity 99.342%, reference coverage 99.99%, HSPs 7; Dudleya lanceolata chloroplast, complete genome
- gb|PX244392.1|: query coverage 78.59%, weighted identity 89.802%, reference coverage 77.83%, HSPs 102; Dudleya greenei chloroplast, complete genome
- gb|PX244394.1|: query coverage 76.19%, weighted identity 89.193%, reference coverage 69.47%, HSPs 135; Dudleya sp. T72600-1 chloroplast, complete genome
### mitochondria_vs_mitochondria_refs
- gb|PV256627.1|: query coverage 70.99%, weighted identity 97.69%, reference coverage 42.20%, HSPs 116; Graptopetalum paraguayense mitochondrion, complete genome
- gb|OP588116.1|: query coverage 68.47%, weighted identity 97.288%, reference coverage 46.93%, HSPs 117; Sedum plumbizincicola voucher KPBK001 mitochondrion, complete genome
- ref|NC_069572.1|: query coverage 68.47%, weighted identity 97.288%, reference coverage 46.93%, HSPs 117; Sedum plumbizincicola voucher KPBK001 mitochondrion, complete genome
- gb|PV608516.1|: query coverage 66.01%, weighted identity 97.426%, reference coverage 61.52%, HSPs 122; Sedum sarmentosum mitochondrion, complete genome
- gb|OP573219.1|: query coverage 64.11%, weighted identity 97.485%, reference coverage 36.56%, HSPs 163; Rhodiola tangutica mitochondrion, complete genome
- ref|NC_072122.1|: query coverage 64.11%, weighted identity 97.485%, reference coverage 36.56%, HSPs 163; Rhodiola tangutica mitochondrion, complete genome
- gb|OP312067.1|: query coverage 57.26%, weighted identity 97.436%, reference coverage 43.32%, HSPs 143; Rhodiola crenulata mitochondrion, complete genome
- ref|NC_070303.1|: query coverage 57.26%, weighted identity 97.436%, reference coverage 43.32%, HSPs 143; Rhodiola crenulata mitochondrion, complete genome
### chloroplast_vs_mitochondria_refs
- gb|PV256627.1|: query coverage 18.92%, weighted identity 92.695%, reference coverage 13.82%, HSPs 25; Graptopetalum paraguayense mitochondrion, complete genome
- gb|PP024540.1|: query coverage 15.50%, weighted identity 93.484%, reference coverage 7.42%, HSPs 25; Rhodiola rosea mitochondrion, complete genome
- gb|OP573219.1|: query coverage 12.47%, weighted identity 92.961%, reference coverage 7.57%, HSPs 14; Rhodiola tangutica mitochondrion, complete genome
- ref|NC_072122.1|: query coverage 12.47%, weighted identity 92.961%, reference coverage 7.57%, HSPs 14; Rhodiola tangutica mitochondrion, complete genome
- gb|OR188139.1|: query coverage 9.41%, weighted identity 91.726%, reference coverage 5.55%, HSPs 20; Rhodiola juparensis mitochondrion, complete genome
- ref|NC_082108.1|: query coverage 9.41%, weighted identity 91.726%, reference coverage 5.55%, HSPs 20; Rhodiola juparensis mitochondrion, complete genome
- gb|OP588116.1|: query coverage 7.69%, weighted identity 90.071%, reference coverage 5.39%, HSPs 14; Sedum plumbizincicola voucher KPBK001 mitochondrion, complete genome
- ref|NC_069572.1|: query coverage 7.69%, weighted identity 90.071%, reference coverage 5.39%, HSPs 14; Sedum plumbizincicola voucher KPBK001 mitochondrion, complete genome
### mitochondria_vs_chloroplast_refs
- gb|OQ076651.1|: query coverage 1.24%, weighted identity 77.676%, reference coverage 1.71%, HSPs 7; Dudleya farinosa chloroplast, complete genome
- ref|NC_085682.1|: query coverage 1.24%, weighted identity 77.676%, reference coverage 1.71%, HSPs 7; Dudleya farinosa chloroplast, complete genome
- gb|PX244390.1|: query coverage 1.24%, weighted identity 77.635%, reference coverage 1.71%, HSPs 7; Dudleya nesiotica chloroplast, complete genome
- gb|PX244393.1|: query coverage 1.24%, weighted identity 77.635%, reference coverage 1.71%, HSPs 7; Dudleya candelabrum chloroplast, complete genome
- gb|PX244389.1|: query coverage 1.24%, weighted identity 77.635%, reference coverage 1.71%, HSPs 7; Dudleya virens chloroplast, complete genome
- gb|PX244391.1|: query coverage 1.24%, weighted identity 77.601%, reference coverage 1.72%, HSPs 7; Dudleya lanceolata chloroplast, complete genome
- gb|PX244394.1|: query coverage 1.19%, weighted identity 75.366%, reference coverage 1.42%, HSPs 9; Dudleya sp. T72600-1 chloroplast, complete genome
- gb|PX244392.1|: query coverage 0.95%, weighted identity 75.506%, reference coverage 1.46%, HSPs 7; Dudleya greenei chloroplast, complete genome

## Direct Dudleya setchellii NCBI check
- NCBI nuccore search for `Dudleya setchellii[Organism]` returned 2 records, but they are short marker/partial records, not complete organelle genomes:
  - JX960535.1: Dudleya abramsii subsp. setchellii trnL gene/trnL-trnF spacer, partial chloroplast sequence, 758 bp.
  - JX960458.1: external transcribed spacer and 18S rRNA partial nuclear/ribosomal-region sequence, 573 bp.
- NCBI search for `Dudleya[Organism] AND mitochondrion[Title]` returned 0 records in this run.

## Bottom-line identity call
- Chloroplast FASTA: strongly supported as a Dudleya chloroplast genome. It covers 99.95-100.00% of several complete Dudleya chloroplast references at ~99.34-99.57% weighted nucleotide identity. Cross-comparison to Crassulaceae mitochondrial references is much weaker (best query coverage 18.92%), consistent with normal plastid-derived segments in plant mitochondrial records rather than the query being mitochondrial.
- Mitochondrial FASTA: supported as a Crassulaceae mitochondrial genome, but the taxonomic claim is necessarily weaker than for the chloroplast because no complete Dudleya mitochondrial reference was found. It covers 57.26-70.99% of several Crassulaceae mitochondrial references at ~97.3-97.7% weighted identity. Cross-comparison to Dudleya chloroplast references is only ~1.24% coverage, arguing against the file being a chloroplast assembly.
- I found no evidence from these comparisons that the two files are swapped.

## Draft annotation summary
- chloroplast: 717 raw homology-supported feature hits; best nonredundant draft: 238 features. Raw by type {'gene': 260, 'tRNA': 59, 'CDS': 386, 'rRNA': 12}; raw confidence {'low': 78, 'medium': 114, 'high': 491, 'weak': 34}. Best nonredundant by type {'gene': 119, 'CDS': 85, 'tRNA': 30, 'rRNA': 4}; best confidence {'high': 224, 'medium': 2, 'low': 11, 'weak': 1}.
  - TSV: /home/neil/godsgate_results/dudleya-organelles/chloroplast.annotation.tsv
  - GFF3: /home/neil/godsgate_results/dudleya-organelles/chloroplast.draft.gff3
  - Best nonredundant TSV: /home/neil/godsgate_results/dudleya-organelles/chloroplast.best_nonredundant.annotation.tsv
  - Best nonredundant GFF3: /home/neil/godsgate_results/dudleya-organelles/chloroplast.best_nonredundant.draft.gff3
- mitochondria: 509 raw homology-supported feature hits; best nonredundant draft: 96 features. Raw by type {'CDS': 269, 'gene': 186, 'tRNA': 38, 'rRNA': 16}; raw confidence {'low': 55, 'high': 370, 'medium': 37, 'weak': 47}. Best nonredundant by type {'CDS': 34, 'gene': 47, 'tRNA': 13, 'rRNA': 2}; best confidence {'high': 79, 'medium': 7, 'low': 5, 'weak': 5}.
  - TSV: /home/neil/godsgate_results/dudleya-organelles/mitochondria.annotation.tsv
  - GFF3: /home/neil/godsgate_results/dudleya-organelles/mitochondria.draft.gff3
  - Best nonredundant TSV: /home/neil/godsgate_results/dudleya-organelles/mitochondria.best_nonredundant.annotation.tsv
  - Best nonredundant GFF3: /home/neil/godsgate_results/dudleya-organelles/mitochondria.best_nonredundant.draft.gff3

## Interpretation limits
- This is a homology-transfer draft, not a curated GenBank-submission annotation.
- Complete Dudleya setchellii organelle genomes were not found as public NCBI references here. Two short D. abramsii subsp. setchellii marker records exist, but chloroplast whole-genome evidence uses other complete Dudleya chloroplasts and mitochondrial evidence uses other Crassulaceae mitochondria.
- Plant mitochondrial genomes can be rearranged and can contain plastid-derived segments, so identity is based on cumulative organelle/reference evidence rather than one perfect collinear hit.
- lncRNA/noncoding RNA annotation is limited to annotated ncRNA/misc_RNA features present in related GenBank records and detectable by sequence similarity; this does not discover novel lncRNAs de novo.

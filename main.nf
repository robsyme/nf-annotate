#!/usr/bin/env nextflow

params.masked = "data/masked/*/*.fasta"
params.rnaseq = "data/masked/*/*.bam"
params.refgenome = "Ptt_W11"
params.altpep = "data/altproteins/*.fasta"
params.maxintronlength = 100

altpep = Channel.fromPath(params.altpep)

Channel.fromPath(params.masked)
.map { Path path -> [path.getParent().getBaseName(), path]}
.tap { genomesForUnmasking }
.tap { genomesForCactus }
.tap { genomesForAragorn }
.tap { genomesForExonerate }

process mash {
  tag { id }

  input:
  set id, 'softmasked.fasta' from genomesForUnmasking

  output:
  file "${id}.msh" into sketches1
  file "${id}.msh" into sketches2

  """
awk '/^>/ {print \$0} !/^>/ {print toupper(\$0)}' softmasked.fasta > ${id}.fasta
mash sketch -p ${task.cpus} -o ${id} ${id}.fasta
  """
}

process distance {
  input:
  set "s1.msh", "s2.msh" from sketches1.combine(sketches2)

  output:
  stdout into distances

  "mash dist s1.msh s2.msh | sed 's/\\.fasta//g'"
}

process makeTree {
  input:
  file "distances.txt" from distances.collectFile()

  output:
  file "tree.nwk" into tree
  file "tree.nwk" into tree2

  """
#!/usr/bin/env Rscript
library("phangorn")
library("reshape2")
library("magrittr")
library("dplyr")

read.table("distances.txt", stringsAsFactors=FALSE, sep="\t", col.names=c("ref", "qry", "distance", "pvalue", "mashmatch")) %>%
    select(ref, qry, distance) %>%
    acast(ref ~ qry) %>%
    as.dist(upper=F, diag=F) %>%
    upgma %>%
    write.tree("tree.nwk")
  """
}

process viewTree {
  input:
  file "tree.nwk" from tree2

  output:
  stdout into asciitree

  "nw_display tree.nwk"
}

asciitree.subscribe {
  println("================================================================================")
  println("============================= Estimated  Phylogeny =============================")
  println()
  print(it)
  println("================================================================================")
}

process renamegenomes {
  input:
  set id, "input.fasta" from genomesForCactus

  output:
  file "${id}.fasta" into genomesRenamedForCactus

  "ln -s input.fasta ${id}.fasta"
}

process predict_tRNA {
  tag { id }

  input:
  set id, 'genome.fasta' from genomesForAragorn

  output:
  set id, 'aragorn.gff3' into trnas

  """
aragorn -t genome.fasta > out
grep -E -C2 '(nucleotides|Sequence)' out > 1
aragorn_to_gff3.lua < 1 > 2
gt gff3 -sort -tidy -retainids 2 > aragorn.gff3
  """
}

process make_exonerate_index {
  tag { id }

  input:
  set id, 'genome.fasta' from genomesForExonerate

  output:
  set id, 'genome.fasta', 'index.esi', 'index.esd' into exn_index

  """
fasta2esd --softmask no genome.fasta index.esd
esd2esi index.esd index.esi --translate yes
  """
}

process exonerate {
  tag { id }

  input:
  set id, 'genome.fasta', 'index.esi', 'index.esd', 'altprop.fasta' from exn_index.combine(altpep.splitFasta( by: 200, file: true))

  output:
  set id, 'exn_out' into exonerate_out

  """
exonerate \
 -E false \
 --model p2g \
 --showvulgar no \
 --showalignment no \
 --showquerygff no \
 --showtargetgff yes \
 --percent 80 \
 --geneseed 250 \
 --ryo \"AveragePercentIdentity: %pi\n\" \
 altprop.fasta \
 genome.fasta > exn_out
  """
}

exonerate_out
.collectFile() { id, path -> ["${id}.hints", path.text] }
.map { path -> [path.getBaseName(), path] }
.set { exonerate_hint_input }

process exonerate_make_hints {
  tag { id }

  input:
  set id, 'exn_out' from exonerate_hint_input

  output:
  set id, 'augustus.hints' into exn_hints

  """
exonerate2hints.pl \
 --source=P --maxintronlen=${params.maxintronlength} \
 --in=exn_out \
 --out=augustus.unsorted.hints
sort -k1,1 -k4,4n augustus.unsorted.hints > augustus.hints
  """
}

process progressivecactus {
  scratch true
  cpus 10

  input:
  file "tree.nwk" from tree
  file genomes from genomesRenamedForCactus.toList()

  output:
  file 'genomes.hal' into hal

  """
cp tree.nwk genomes.txt
for fasta in *.fasta; do
  echo \$(basename \$fasta .fasta) \$fasta
done >> genomes.txt
runProgressiveCactus.sh --maxThreads=${task.cpus} genomes.txt outdir genomes.hal
  """
}

process hal2maf {
  input:
  file 'genomes.hal' from hal

  output:
  file 'mafs' into mafs

  "hal2maf_split.pl --halfile genomes.hal --refGenome ${params.refgenome} --cpus 8 --chunksize 50000 --overlap 25000 --outdir mafs"
}

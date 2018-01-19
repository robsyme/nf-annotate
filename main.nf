#!/usr/bin/env nextflow

params.masked = "data/masked/*/*.fasta"
params.rnaseq = "data/masked/*/*.bam"
params.refgenome = "Ptt_W11"

Channel.fromPath(params.masked)
.map { Path path -> [path.getParent().getBaseName(), path]}
.tap { genomesForUnmasking }
.tap { genomesForCactus }
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


process progressivecactus {
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

process predict_tRNA {
  tag { id }

  input:
  set id, 'genome.fasta' from genomesForExonerate

  output:
  set id, 'aragorn.gff3' into trnas

  """
aragorn -t genome.fasta > out
grep -E -C2 '(nucleotides|Sequence)' out > 1
aragorn_to_gff3.lua < 1 > 2
gt gff3 -sort -tidy -retainids 2 > aragorn.gff3
  """
}

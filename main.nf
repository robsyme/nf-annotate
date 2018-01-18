#!/usr/bin/env nextflow

params.masked = "data/*/*.fasta"
params.rnaseq = "data/*/*.bam"

Channel.fromPath(params.masked)
.map { Path path -> [path.getParent().getBaseName(), path]}
.tap { genomesForUnmasking }

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

tree.println()





import "finemap_meta_sub.wdl" as sub

workflow liftover_pheno {

    String zones
    String docker
    File phenolistfile
    String pheno_format

    Array[String] phenos = read_lines(phenolistfile)

    scatter (pheno in phenos) {
        String table = sub(pheno_format, "\\{PHENO\\}", pheno)
        call sub.liftover as liftover {
            input: zones=zones, docker=docker, table=table,
            outname=basename(table, ".txt") + ".liftover.txt"
        }
    }
}

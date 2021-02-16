task liftover {
    File table
    File script
    String reference_genome
    String variant_col
    String outname
    String zones
    String docker
    Int cpu
    Int mem

    command <<<
        # pip3 install cython && pip3 install gnomad
        curl -O https://raw.githubusercontent.com/broadinstitute/gnomad_methods/master/gnomad/utils/liftover.py

        PYSPARK_SUBMIT_ARGS="--conf spark.driver.memory=${mem}g pyspark-shell" \
        python3 ${script} \
        --table ${table} \
        --reference-genome ${reference_genome} \
        --variant-col ${variant_col} \
        --out ${outname}

    >>>

    output {
        File out = outname
    }

    runtime {

        docker: "${docker}"
        cpu: "${cpu}"
        memory: "${mem} GB"
        disks: "local-disk 20 HDD"
        zones: "${zones}"
        preemptible: 2
        noAddress: false
    }
}

task abf {
    File zfile
    String pheno
    String prefix = basename(zfile, ".z")
    Float prior_variance
    Float coverage
    String zones
    String docker
    Int cpu
    Int mem

    command <<<
        run_abf.R \
            --z ${zfile} \
            --snp ${prefix}.abf.snp \
            --cred ${prefix}.abf.cred \
            --log ${prefix}.abf.log \
            --prior-variance ${prior_variance}

    >>>

    output {
        File log = prefix + ".abf.log"
        File snp = prefix + ".abf.snp"
        File cred = prefix + ".abf.cred"

    }

    runtime {

        docker: "${docker}"
        cpu: "${cpu}"
        memory: "${mem} GB"
        disks: "local-disk 5 HDD"
        zones: "${zones}"
        preemptible: 2
        noAddress: false
    }
}

task susie {
    Int n_samples
    Float var_y
    Int n_causal_snps
    File zfile
    String pheno
    String prefix = basename(zfile, ".z")
    String zones
    String docker
    Int cpu
    Int mem
    Float min_cs_corr

    command <<<
        #!/usr/bin/env bash

        run_susieR.R \
            --z ${zfile} \
            -n ${n_samples} \
            --L ${n_causal_snps} \
            --var-y ${var_y} \
            --snp ${prefix}.susie.snp \
            --cred ${prefix}.susie.cred \
            --log ${prefix}.susie.log \
            --susie-obj ${prefix}.susie.rds \
            --save-susie-obj \
            --write-alpha \
            --write-single-effect \
            --min-cs-corr ${min_cs_corr}

    >>>

    output {
        File log = prefix + ".susie.log"
        File snp = prefix + ".susie.snp"
        File cred = prefix + ".susie.cred"
        File rds = prefix + ".susie.rds"

    }

    runtime {

        docker: "${docker}"
        cpu: "${cpu}"
        memory: "${mem} GB"
        disks: "local-disk 5 HDD"
        zones: "${zones}"
        preemptible: 2
        noAddress: false
    }
}


task annotate_gamma {
    File snp
    File sumstats
    File betafile
    File script
    String outname = basename(snp, ".bgz")
    String zones
    String docker
    Int cpu
    Int mem

    command <<<

        Rscript ${script} \
        --snp ${snp} \
        --beta ${betafile} \
        --meta ${sumstats} \
        --out ${outname}

        bgzip -c ${outname} > ${outname}.bgz

    >>>

    output {
        File out = outname + ".bgz"

    }

    runtime {

        docker: "${docker}"
        cpu: "${cpu}"
        memory: "${mem} GB"
        disks: "local-disk 5 HDD"
        zones: "${zones}"
        preemptible: 2
        noAddress: false
    }
}

task combine {
    String pheno
    Int n_causal_snps
    Array[File] abf_snp
    Array[File] abf_cred
    String zones
    String docker
    Int cpu
    Int mem

    command <<<

        cat << "__EOF__" > combine_snp.awk
        BEGIN {
            OFS = "\t"
        }
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                col[$i] = i
            }
            gsub(" ", "\t")
            print "trait", "region", "variant", $0
        }
        FNR == 1 {
            match(FILENAME, /(chr[0-9X]+)\.([0-9]+-[0-9]+)\./, a)
            region = a[1]":"a[2]
        }
        FNR > 1 {
            chrom = substr($col["chromosome"], 4)
            sub(/^0/, "", chrom)
            v = sprintf( \
                "%s:%s:%s:%s", \
                chrom, \
                $col["position"], \
                $col["allele1"], \
                $col["allele2"] \
            )
            gsub(" ", "\t")
            print pheno, region, v, $0 | "sort -V -k2,3"
        }
        __EOF__


        cat << "__EOF__" > combine_cred.awk
        BEGIN {
            OFS = "\t"
        }
        NR == 1 {
            print "trait", "region", $0
        }
        FNR == 1 {
            match(FILENAME, /(chr[0-9X]+)\.([0-9]+-[0-9]+)\./, a)
            region = a[1]":"a[2]
        }
        FNR > 1 {
            print pheno, region, $0
        }
        __EOF__

        # Combine abf .snp files
        awk -f combine_snp.awk -v pheno=${pheno} ${sep=" " abf_snp} | bgzip -c -@ ${cpu} > ${pheno}.ABF.snp.bgz
        tabix -s 5 -b 6 -e 6 -S 1 ${pheno}.ABF.snp.bgz
        # Combine abf .cred files
        awk -f combine_cred.awk -v pheno=${pheno}  ${sep=" " abf_cred} | bgzip -c -@ ${cpu} > ${pheno}.ABF.cred.bgz

    >>>

    output {

        File out_abf_snp = pheno + ".ABF.snp.bgz"
        File out_abf_snp_tbi = pheno + ".ABF.snp.bgz.tbi"
        File out_abf_cred = pheno + ".ABF.cred.bgz"

    }

    runtime {

        docker: "${docker}"
        cpu: "${cpu}"
        memory: "${mem} GB"
        disks: "local-disk 30 HDD"
        zones: "${zones}"
        preemptible: 2
        noAddress: false
    }
}

workflow finemap {

    String zones
    String docker
    String pheno
    Int n_samples
    Float var_y
    Int n_causal_snps
    Array[File] zfiles

    File sumstats
    File betafile

    scatter (zfile in zfiles) {

        call abf {
            input: zones=zones, docker=docker, zfile=zfile, pheno=pheno
        }

        # call annotate as annotate_abf {
        #     input: zones=zones, docker=hail_docker, cpu=hail_cpu, mem=hail_mem,
        #         annotation_script=annotation_script, snp=abf.snp
        # }

    }

    call combine {
        input: zones=zones, docker=docker, pheno=pheno, n_causal_snps=n_causal_snps,
            abf_snp=abf.snp, abf_cred=abf.cred
            # abf_snp=annotate_abf.ld_snp, abf_cred=abf.cred
    }

    call annotate_gamma {
        input: zones=zones, docker=docker, snp=combine.out_abf_snp,
            sumstats=sumstats, betafile=betafile
    }
}

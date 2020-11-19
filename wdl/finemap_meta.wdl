import "finemap_meta_sub.wdl" as sub

task meta_analysis {
    Pair[String, String] config_pheno
    String config = config_pheno.left
    String pheno = config_pheno.right
    File config_file
    String sumstats_format
    String out_format = "meta.pheno%d.config%d.tsv"
    File script

    String zones
    String docker
    Int cpu
    Int mem

    command <<<

        R --slave --vanilla -e 'install.packages("furrr", repos = structure(c(CRAN="https://cloud.r-project.org/")))'

        Rscript ${script} \
        --config ${config} \
        --pheno ${pheno} \
        --config-file ${config_file} \
        --sumstats-format ${sumstats_format} \
        --out-format ${out_format} \
        --threads ${cpu}
    >>>

    output {

        File sumstats = "meta.pheno" + pheno + ".config" + config + ".tsv.gz"
    }

    runtime {
        docker: "${docker}"
        cpu: "${cpu}"
        memory: "52 GB"
        disks: "local-disk 20 HDD"
        zones: "${zones}"
        preemptible: 2
        noAddress: false
    }
}

task preprocess {
    File sumstats
    String pheno = basename(sumstats, ".tsv.gz")
    String zones
    String docker
    Int cpu
    Int mem
    String rsid_col
    String chromosome_col
    String position_col
    String allele1_col
    String allele2_col
    String freq_col
    String beta_col
    String se_col
    String p_col
    String delimiter
    Float p_threshold
    Boolean? grch38

    command <<<

        zcat ${sumstats} | awk '
        BEGIN {
            OFS = "\t"
            n_samples = 0
        }
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                col[$i] = i
            }
        }
        NR > 1 {
            if ($col["n_samples"] > n_samples) {
                n_samples = $col["n_samples"]
            }
        }
        END {
            print n_samples
        }
        ' > n_samples.txt

        make_finemap_inputs.py \
            --sumstats ${sumstats} \
            --rsid-col "${rsid_col}" \
            --chromosome-col "${chromosome_col}" \
            --position-col "${position_col}" \
            --allele1-col "${allele1_col}" \
            --allele2-col "${allele2_col}" \
            --freq-col "${freq_col}" \
            --beta-col "${beta_col}" \
            --se-col "${se_col}" \
            --p-col "${p_col}" \
            --extra-cols "p_het" "n_studies" "n_samples" \
            --delimiter "${delimiter}" \
            ${true='--grch38 ' false=' ' grch38} \
            --no-upload \
            --prefix ${pheno} \
            --out ${pheno} \
            --wdl \
            --p-threshold ${p_threshold}

        res=`cat ${pheno}_had_results`

        if [ "$res" == "False" ]
        then
            touch ${pheno}".z"
            touch ${pheno}".lead_snps.txt"
            touch ${pheno}".bed"
        fi

    >>>

    output {

        String out_pheno = pheno
        Array[File] zfiles = glob("*.z")
        File leadsnps = pheno + ".lead_snps.txt"
        File bed = pheno + ".bed"
        File log = pheno + ".log"
        Boolean had_results = read_boolean("${pheno}_had_results")
        Int n_samples = read_int("n_samples.txt")
        Float var_y = 1
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

workflow finemap_meta {

    String zones
    String docker
    File configlistfile
    File phenolistfile
    String betafile_format

    Array[String] configs = read_lines(configlistfile)
    Array[String] phenos = read_lines(phenolistfile)

    Array[Pair[String, String]] all_pairs = cross(configs, phenos)

    scatter (config_pheno in all_pairs) {

        call meta_analysis {
            input: zones=zones, docker=docker, config_pheno = config_pheno
        }

        call preprocess {
            input: zones=zones, docker=docker, sumstats=meta_analysis.sumstats
        }

        if(preprocess.had_results) {
            File betafile = sub(betafile_format, "\\{PHENO\\}", config_pheno.right)
            call sub.finemap {
                input: zones=zones, docker=docker, pheno=preprocess.out_pheno, zfiles=preprocess.zfiles,
                    n_samples=preprocess.n_samples, var_y=preprocess.var_y,
                    sumstats=meta_analysis.sumstats, betafile=betafile
            }
        }

    }
}

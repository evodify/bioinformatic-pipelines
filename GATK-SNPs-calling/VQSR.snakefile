# call CNVs with GATK 4

REF = 'data/reference/canFam3.fa'
MEM = '12G'

rule all:
    input:
        'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR_SNPs_VQSR.vcf',
        'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR_SNPs_VQSRfilterPASSED.vcf',
        'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR_INDELs_VQSR.vcf',
        'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR_INDELs_VQSRfilterPASSED.vcf',
        'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR_SNPs_VQSR.plots.R',
        'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR_SNPs_VQSR.plots.R.pdf'

rule select_SNPs:
    input:
        ref = REF,
        vcf = 'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR.vcf.gz',
        index = 'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR.vcf.gz.tbi'
    output:
        'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR_SNPs.vcf.gz'
    params: MEM
    shell:
        '''
        gatk --java-options "-Xmx{params} -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" SelectVariants  \
        -R {input.ref} \
        -V {input.vcf}  \
        --select-type-to-include SNP \
        -O {output}
        '''

rule select_INDELs:
    input:
        ref = REF,
        vcf = 'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR.vcf.gz',
        index = 'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR.vcf.gz.tbi'
    output:
        'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR_INDELs.vcf.gz'
    params: MEM
    shell:
        '''
        gatk --java-options "-Xmx{params} -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" SelectVariants  \
        -R {input.ref} \
        -V {input.vcf}  \
        --select-type-to-include INDEL \
        -O {output}
        '''

rule VariantRecalibrator_SNPs:
    input:
        ref = REF,
        vcf = 'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR_SNPs.vcf.gz',
        snps1 = 'data/reference/BQSRreference/HD_array_segregating_canfam3.vcf.gz',
        snps2 = 'data/reference/BQSRreference/Axelsson_2013_SNPs_canfam3.vcf.gz',
        snps3 = 'data/reference/BQSRreference/160DG.99.9.recalibrated_variants.PASS.vcf.gz',
        snps4 = 'data/reference/BQSRreference/00-All_chrAll.vcf.gz',
        snps5 = 'data/reference/BQSRreference/dogs.557.publicSamples.ann.PASS.vcf.gz',
        snps6 = 'data/reference/BQSRreference/wolf.concat.raw.SNPs.filterPASSED.vcf.gz'
    output:
        recal = 'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR_SNPs.recal',
        tranches = 'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR_SNPs.tranches',
        r = 'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR_SNPs_VQSR.plots.R',
        pdf = 'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR_SNPs_VQSR.plots.R.pdf'
    params: MEM
    shell:
        '''
        gatk --java-options "-Xmx{params} -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" VariantRecalibrator \
        -R {input.ref} \
        -V {input.vcf}  \
        --resource:HDarray,known=false,training=true,truth=true,prior=15.0 {input.snps1} \
        --resource:Axelsson_2013,known=false,training=true,truth=false,prior=12.0 {input.snps2} \
        --resource:dogs160,known=false,training=false,truth=false,prior=10.0 {input.snps3} \
        --resource:00-All,known=false,training=false,truth=false,prior=2.0 {input.snps4} \
        --resource:dogs.557,known=false,training=false,truth=false,prior=2.0 {input.snps5} \
        --resource:wolf,known=false,training=true,truth=false,prior=10.0 {input.snps6} \
        -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP -mode SNP \
        -O {output.recal}  \
        --tranches-file {output.tranches} \
        --rscript-file {output.r}
        '''

rule VariantRecalibrator_INDELs:
    input:
        ref = REF,
        vcf = 'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR_INDELs.vcf.gz',
        snps2 = 'data/reference/BQSRreference/Axelsson_2013_indels_canfam3.vcf.gz',
        snps3 = 'data/reference/BQSRreference/160DG.99.9.recalibrated_variants.PASS.vcf.gz',
        snps4 = 'data/reference/BQSRreference/00-All_chrAll.vcf.gz',
        snps5 = 'data/reference/BQSRreference/dogs.557.publicSamples.ann.PASS.vcf.gz',
        snps5 = 'data/reference/BQSRreference/wolf.concat.raw.INDEL.filterPASSED.vcf.gz'
    output:
        recal = 'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR_INDELs.recal',
        tranches = 'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR_INDELs.tranches',
        r = 'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR_INDELs.plots.R'
    params: MEM
    shell:
        '''
        gatk --java-options "-Xmx{params} -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" VariantRecalibrator \
        -R {input.ref} \
        -V {input.vcf}  \
        --resource:Axelsson_2013,known=false,training=false,truth=false,prior=6.0 {input.snps2} \
        --resource:dogs160,known=false,training=true,truth=true,prior=10.0 {input.snps3} \
        --resource:00-All,known=false,training=false,truth=false,prior=2.0 {input.snps4} \
        --resource:dogs.557,known=false,training=false,truth=true,prior=6.0 {input.snps5} \
        --resource:wolf,known=false,training=true,truth=true,prior=6.0
        -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP -mode INDEL \
        -O {output.recal}  \
        --tranches-file {output.tranches} \
        --rscript-file {output.r}
        '''

rule ApplyVQSR_SNPs:
    input:
        ref = REF,
        vcf = 'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR_SNPs.vcf.gz',
        recal = 'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR_SNPs.recal',
        tranches = 'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR_SNPs.tranches'
    output:
        vcf = 'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR_SNPs_VQSR.vcf'
    params: MEM
    shell:
        '''
        gatk --java-options "-Xmx{params} -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" ApplyVQSR \
        -R {input.ref} \
        -V {input.vcf} \
        -O Agrarian_vs_Arctic_DNA_markDupl_BQSR_SNPs_VQSR.vcf \
        --truth-sensitivity-filter-level 99.5 \
        --tranches-file {input.tranches} \
        --recal-file {input.recal} \
        -mode SNP
        '''

rule ApplyVQSR_INDELs:
    input:
        ref = REF,
        vcf = 'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR_INDELs.vcf.gz',
        recal = 'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR_INDELs.recal',
        tranches = 'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR_INDELs.tranches'
    output:
        vcf = 'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR_INDELs_VQSR.vcf'
    params: MEM
    shell:
        '''
        gatk --java-options "-Xmx{params} -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" ApplyVQSR \
        -R {input.ref} \
        -V {input.vcf} \
        -O Agrarian_vs_Arctic_DNA_markDupl_BQSR_INDELs_VQSR.vcf \
        --truth-sensitivity-filter-level 99.5 \
        --tranches-file {input.tranches} \
        --recal-file {input.recal} \
        -mode INDEL
        '''

rule filterVQSR_SNPs:
    input:
        'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR_SNPs_VQSR.vcf'
    output:
        'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR_SNPs_VQSRfilterPASSED.vcf'
    shell:
        '''
        grep -E "CHROM|PASS" {input} > {output}
        '''

rule filterVQSR_INDELs:
    input:
        'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR_INDELs_VQSR.vcf'
    output:
        'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR_INDELs_VQSRfilterPASSED.vcf'
    shell:
        '''
        grep -E "CHROM|PASS" {input} > {output}
        '''

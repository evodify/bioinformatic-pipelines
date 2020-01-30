# map reads with BWA and call variants with GATK

CHR = list(range(1, 39))
CHR.append('X')
CHR.append('M')

# If you have many small scaffolds they can be loaded with:
CHRUN = list()
with open('data/meta/chrUn.list', 'r') as file:
  for line in file:
    CHRUN.append(line.rstrip('\n'))
CHRUNL = ' -L '.join(CHRUN)

SAMPLES = ['sample1', 'sample2', 'sampleN']
REF = '/path/to/reference/canFam3.fa'
CORES = 12
MEM = '12G'
MEMCORE = '2G'

rule all:
    input:
        expand('bam/{sample}_DNA_markDupl_BQSR.bam', sample=SAMPLES),
        expand('qualimap/{sample}_qualimap', sample=SAMPLES),
        'qualimap/mean_genome_coverage.txt',
        expand('VCF/BQSR/{sample}_DNA_markDupl_BQSR.pdf', sample=SAMPLES),
        expand('VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR_chr{j}.vcf.gz', j=CHR),
        'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR_chrUn.vcf.gz',
        'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR.vcf.gz'

rule indexBWA:
    input:
        ref = REF
    output:
        rpac ='/path/to/reference/canFam3.fa.rpac',
        amb = '/path/to/reference/canFam3.fa.amb',
        ann = '/path/to/reference/canFam3.fa.ann',
        pac = '/path/to/reference/canFam3.fa.pac',
        bwt = '/path/to/reference/canFam3.fa.bwt',
        rbwt = '/path/to/reference/canFam3.fa.rbwt',
        rsa = '/path/to/reference/canFam3.fa.rsa',
        sa = '/path/to/reference/canFam3.fa.sa'
    shell:
        '''
        bwa index {input}
        '''

rule indexDict:
    input:
        REF
    output:
        '/path/to/reference/canFam3.dict'
    shell:
        '''
        gatk --java-options "-Xmx{params}" CreateSequenceDictionary -R {input}
        '''

rule indexFai:
    input:
        REF
    output:
        '/path/to/reference/canFam3.fa.fai'
    params: MEM
    shell:
        '''
        samtools faidx {input}
        '''

rule bwa:
    input:
        ref = REF,
        refindex1 = '/path/to/reference/canFam3.fa.bwt',
        refindex2 = '/path/to/reference/canFam3.fa.rpac',
        refindex3 = '/path/to/reference/canFam3.fa.amb',
        refindex4 = '/path/to/reference/canFam3.fa.ann',
        refindex5 = '/path/to/reference/canFam3.fa.pac',
        refindex6 = '/path/to/reference/canFam3.fa.bwt',
        refindex7 = '/path/to/reference/canFam3.fa.rbwt',
        refindex8 = '/path/to/reference/canFam3.fa.bwt',
        refindex9 = '/path/to/reference/canFam3.fa.rsa',
        refindex10 = '/path/to/reference/canFam3.fa.sa',
        R1 = 'data/raw/fastq/{sample}_R1.fastq.gz',
        R2 = 'data/raw/fastq/{sample}_R2.fastq.gz'
    output:
        'bam/{sample}_DNA.bam'
    params:
        mem = MEMCORE,
        samp = '{sample}'
    threads: CORES
    shell:
        '''
        bwa mem -t {threads} -M -R '@RG\\tID:{params.samp}\\tPL:illumina\\tLB:{params.samp}\\tSM:{params.samp}' \
        {input.ref} {input.R1} {input.R2} | samtools sort -O bam -@ {threads} -m {params.mem} > {output}
        '''

rule markDupl:
    input:
        'bam/{sample}_DNA.bam'
    output:
        bam = 'bam/{sample}_DNA_markDupl.bam',
        metrics = 'VCF/metrics/{sample}_DNA_markDupl_metrix.txt'
    params:
        mem = MEM,
        tmp = 'scratch'
    shell:
        '''
        java -Xmx{params.mem} -Djava.io.tmpdir={params.tmp} -jar picard.jar MarkDuplicates \
        VALIDATION_STRINGENCY=LENIENT \
        METRICS_FILE={output.metrics} \
        MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=15000 \
        INPUT={input} \
        OUTPUT={output.bam}
        '''

rule bqsr:
    input:
        ref = REF,
        dic = '/path/to/reference/canFam3.dict',
        fai = '/path/to/reference/canFam3.fa.fai',
        bam = 'bam/{sample}_DNA_markDupl.bam',
        snps1 = '/path/to/reference/BQSRreference/HD_array_segregating_canfam3.vcf.gz',
        snps2 = '/path/to/reference/BQSRreference/Axelsson_2013_SNPs_canfam3.vcf.gz',
        snps3 = '/path/to/reference/BQSRreference/Axelsson_2013_indels_canfam3.vcf.gz',
        snps4 = '/path/to/reference/BQSRreference/160DG.99.9.recalibrated_variants.PASS.vcf.gz',
        snps5 = '/path/to/reference/BQSRreference/00-All_chrAll.vcf.gz',
        snps6 = '/path/to/reference/BQSRreference/dogs.557.publicSamples.ann.PASS.vcf.gz',
        snps7 = '/path/to/reference/BQSRreference/wolf.concat.raw.INDEL.filterPASSED.vcf.gz',
        snps8 = '/path/to/reference/BQSRreference/wolf.concat.raw.SNPs.filterPASSED.vcf.gz'
    output:
        'bam/{sample}_DNA_markDupl_BQSR1.table'
    params: MEM
    shell:
        '''
        gatk --java-options "-Xmx{params}" BaseRecalibrator \
        -R {input.ref} \
        -I {input.bam} \
        --known-sites {input.snps1} \
        --known-sites {input.snps2} \
        --known-sites {input.snps3} \
        --known-sites {input.snps4} \
        --known-sites {input.snps5} \
        --known-sites {input.snps6} \
        --known-sites {input.snps7} \
        --known-sites {input.snps8} \
        -O {output}
        '''

rule apply_bqsr:
    input:
        ref = REF,
        dic = '/path/to/reference/canFam3.dict',
        fai = '/path/to/reference/canFam3.fa.fai',
        bam = 'bam/{sample}_DNA_markDupl.bam',
        table = 'bam/{sample}_DNA_markDupl_BQSR1.table',
    output:
        'bam/{sample}_DNA_markDupl_BQSR.bam'
    params:  MEM
    shell:
        '''
        gatk --java-options "-Xmx{params}" ApplyBQSR \
        -R {input.ref} \
        -I {input.bam} \
        -bqsr {input.table} \
        -O {output}
        '''

rule check_bqsr:
    input:
        ref = REF,
        dic = '/path/to/reference/canFam3.dict',
        fai = '/path/to/reference/canFam3.fa.fai',
        bam = 'bam/{sample}_DNA_markDupl_BQSR.bam',
        snps1 = '/path/to/reference/BQSRreference/HD_array_segregating_canfam3.vcf.gz',
        snps2 = '/path/to/reference/BQSRreference/Axelsson_2013_SNPs_canfam3.vcf.gz',
        snps3 = '/path/to/reference/BQSRreference/Axelsson_2013_indels_canfam3.vcf.gz',
        snps4 = '/path/to/reference/BQSRreference/160DG.99.9.recalibrated_variants.PASS.vcf.gz',
        snps5 = '/path/to/reference/BQSRreference/00-All_chrAll.vcf.gz',
        snps6 = '/path/to/reference/BQSRreference/dogs.557.publicSamples.ann.PASS.vcf.gz',
        snps7 = '/path/to/reference/BQSRreference/wolf.concat.raw.INDEL.filterPASSED.vcf.gz',
        snps8 = '/path/to/reference/BQSRreference/wolf.concat.raw.SNPs.filterPASSED.vcf.gz'
    output:
        'bam/{sample}_DNA_markDupl_BQSR2.table'
    params: MEM
    shell:
        '''
        gatk --java-options "-Xmx{params}" BaseRecalibrator \
        -R {input.ref} \
        -I {input.bam} \
        --known-sites {input.snps1} \
        --known-sites {input.snps2} \
        --known-sites {input.snps3} \
        --known-sites {input.snps4} \
        --known-sites {input.snps5} \
        --known-sites {input.snps6} \
        --known-sites {input.snps7} \
        --known-sites {input.snps8} \
        -O {output}
        '''

rule plot_bqsr:
    input:
        table1 = 'bam/{sample}_DNA_markDupl_BQSR1.table',
        table2 = 'bam/{sample}_DNA_markDupl_BQSR2.table',
    output:
        'VCF/BQSR/{sample}_DNA_markDupl_BQSR.pdf'
    params: MEM
    shell:
        '''
        gatk --java-options "-Xmx{params}"  AnalyzeCovariates \
        -before {input.table1} \
        -after {input.table2} \
        -plots {output}
        '''

rule qualimap:
    input:
        bam = 'bam/{sample}_DNA_markDupl_BQSR.bam',
    output:
        dir = 'qualimap/{sample}_qualimap',
        res = 'qualimap/{sample}_qualimap/genome_results.txt'
    params: MEM
    threads: CORES
    shell:
        '''
        qualimap bamqc -nt {threads} --java-mem-size={params} -bam {input} -outdir {output.dir}
        '''

rule getMeanCoverage:
    input:
        expand('qualimap/{sample}_qualimap/genome_results.txt', sample=SAMPLES),
    output:
        'qualimap/mean_genome_coverage.txt'
    shell:
        '''
        grep "mean coverageData" {input} > {output}
        '''

rule gVCF:
    input:
        ref = REF,
        dic = '/path/to/reference/canFam3.dict',
        fai = '/path/to/reference/canFam3.fa.fai',
        bam = 'bam/{sample}_DNA_markDupl_BQSR.bam',
    output:
        vcf = 'gVCF/{sample}_DNA_markDupl_BQSR.g.vcf.gz',
        index = 'gVCF/{sample}_DNA_markDupl_BQSR.g.vcf.gz.tbi'
    params: MEM
    shell:
        '''
        gatk --java-options "-Xmx{params}" HaplotypeCaller \
        -R {input.ref} \
        -ERC GVCF \
        -I {input.bam} \
        -O {output.vcf}
        '''

rule combineGVCF:
    input:
        ref = REF,
        dic = '/path/to/reference/canFam3.dict',
        fai = '/path/to/reference/canFam3.fa.fai',
        vcf = expand('gVCF/{sample}_DNA_markDupl_BQSR.g.vcf.gz', sample=SAMPLES),
    output:
        vcf = 'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR.g.vcf.gz',
        index = 'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR.g.vcf.gz.tbi'
    params:
        mem = MEM,
        vcfs = lambda wildcards, input: ' -V '.join(input.vcf), 
    shell:
        '''
        gatk --java-options "-Xmx{params.mem} -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" CombineGVCFs \
        -R {input.ref} \
        -V {params.vcfs} \
        -O {output.vcf}
        '''

# rule genotypeGVCFs:
#     input:
#         ref = REF,
#         dic = '/path/to/reference/canFam3.dict',
#         fai = '/path/to/reference/canFam3.fa.fai',
#         vcf = 'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR.g.vcf.gz',
#     output:
#         vcf = 'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR.vcf.gz',
#         index = 'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR.vcf.gz.tbi'
#     params: MEM
#     shell:
#         '''
#         gatk --java-options "-Xmx{params} -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" GenotypeGVCFs \
#         -R {input.ref} \
#         -V {input.vcf} \
#         -O {output.vcf}
#         '''
#
# NOTE! Due to large size of the dog genome it is more efficient to genotype each
# chromosome separately and merge the results later.

rule genotypeGVCFs:
    input:
        ref = REF,
        dic = '/path/to/reference/canFam3.dict',
        fai = '/path/to/reference/canFam3.fa.fai',
        vcf = 'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR.g.vcf.gz',
    output:
        vcf = 'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR_chr{j}.vcf.gz',
        index = 'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR_chr{j}.vcf.gz.tbi'
    params:
        mem = MEM,
        chrN = 'chr{j}'
    shell:
        '''
        gatk --java-options "-Xmx{params.mem} -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" GenotypeGVCFs \
        -R {input.ref} \
        -L {params.chrN} \
        -V {input.vcf} \
        -O {output.vcf}
        '''

rule genotypeGVCFs_chrUn:
    input:
        ref = REF,
        dic = '/path/to/reference/canFam3.dict',
        fai = '/path/to/reference/canFam3.fa.fai',
        vcf = 'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR.g.vcf.gz',
    output:
        vcf = 'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR_chrUn.vcf.gz',
        index = 'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR_chrUn.vcf.gz.tbi'
    params:
        mem = MEM,
        chrN = CHRUNL
    shell:
        '''
        gatk --java-options "-Xmx{params.mem} -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" GenotypeGVCFs \
        -R {input.ref} \
        -L {params.chrN} \
        -V {input.vcf} \
        -O {output.vcf}
        '''

rule mergeVCFs:
    input:
        vcf = expand('VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR_chr{j}.vcf.gz', j=CHR),
        vcfUn = 'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR_chrUn.vcf.gz'
    output:
        'VCF/Agrarian_vs_Arctic_DNA_markDupl_BQSR.vcf.gz'
    params:
        mem = MEM,
        temp = '$SNIC_TMP',
        vcfs = lambda wildcards, input: ' I='.join(input.vcf)
    shell:
        '''
        java -Xmx{params.mem} -Djava.io.tmpdir={params.temp} -jar picard.jar MergeVcfs \
        I={params.vcfs} \
        I={input.vcfUn} \
        O={output}
        '''

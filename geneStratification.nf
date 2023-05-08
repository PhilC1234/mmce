nextflow.enable.dsl=2

process tabix {
    input:
        tuple val(meta), path(vcf)

    output:
        tuple val(meta), path("*vcf.gz"), emit: vcf
        path "*vcf.gz*", emit: publishFiles

    executor = 'float'
    image = 'tabix'
    cpu = 2
    mem = 2
    publishDir path: "${params.outputFolder}", mode: 'move', overwrite: 'true'
    
    script:
    """
        bgzip ${vcf}
        tabix ${vcf}.gz
    """
}

workflow METADATA {
    
    take:
        f
        
    main:
        fileTuple = f.map{ [ [ id: it.simpleName, baseName: it.baseName, size: it.size(), folder: it.parent, fileType: it.extension, fullFile: it ], it ] }

    emit:
        fileTuple

}

process sam_index {
    input:
        tuple val(sample), path(bam)

    output:
        path "${bam}*" 

    executor = 'float'
    image = 'bcftools'
    cpu = 6
    mem = 2

    publishDir path: "${params.outputFolder}", mode: 'move', overwrite: 'true'

    script:
        """
            samtools index "${bam}" -@ ${task.cpus} "${bam}.bai" 
        """

}

process intersectBam {

    input:
        tuple val(sample), path(bam), val(bed), path(bedFile)

    output:
        path "${sample.id}----${bed}.bam"

    tag "${bam}----${bed}"

    executor = 'float'
    image = 'bedtools'
    cpu = 2
    mem = 2

    script:
    """
        intersectBed -a ${bam} -b ${bedFile} > "${sample.id}----${bed}.bam"
    """
}

process intersectVcf {

    input:
        tuple val(sample), path(vcf), val(bed), path(bedFile)

    output:
        path "${sample.id}----${bed}.vcf"

    executor = 'float'
    image = 'bedtools'
    //extra = '--migratePolicy [enable=true]'
    extra = '--migratePolicy [cpu.upperBoundRatio=80,cpu.lowerBoundRatio=10,cpu.upperBoundDuration=3s,cpu.lowerBoundDuration=30s,cpu.step=50,mem.upperBoundRatio=80,mem.lowerBoundRatio=10,mem.upperBoundDuration=3s,mem.lowerBoundDuration=30s,mem.step=50]'
    cpu = 2
    mem = 4
//    errorStrategy 'retry' 
    tag "${vcf}----${bed}"
    
    script:
    """
        intersectBed -a ${vcf} -b ${bedFile} -header > "${sample.id}----${bed}.vcf"
    """
}

workflow {
    beds  = Channel.fromPath(params.bedTSV).splitCsv(header: false, sep: "\t")
    beds.map { [ bed: it[0], bedFile: it[1] ] }
//    bams = Channel.fromPath(params.bams) | METADATA
    vcfs = Channel.fromPath(params.vcfs) | METADATA
//    bams = bams.combine(beds).map { [ sample: it[0], bam: it[1], bed: it[2], bedFile: 's3://bdtest2003/tgen/GRCh38/' + it[3] ] }
    vcfs = vcfs.combine(beds).map { [ sample: it[0], vcf: it[1], bed: it[2], bedFile: '/home/ec2-user/mnt/s3/GRCh38/' + it[3] ] }

//    sam_index(intersectBam(bams) | METADATA)
    tabix(intersectVcf(vcfs) | METADATA)
}

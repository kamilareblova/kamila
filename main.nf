include { ALIGN; VARDICT; NORMALIZACE; ANOTACE; COVERAGE } from "${params.projectDirectory}/modules"


workflow {
rawfastq = Channel.fromPath("${params.homeDir}/samplesheet.csv")
    . splitCsv( header:true )
    . map { row ->
        def meta = [name:row.name, run:row.run]
        [meta.name, meta, [
            file("${params.inputDirectory}/${meta.name}_R1.fastq.gz"),
            file("${params.inputDirectory}/${meta.name}_R2.fastq.gz"),
        ]]
    }
     . view()

aligned	= ALIGN(rawfastq)
varcalling = VARDICT(aligned)
normalizovany = NORMALIZACE(varcalling)
anotovany = ANOTACE(normalizovany)
coverage = COVERAGE(aligned)
}

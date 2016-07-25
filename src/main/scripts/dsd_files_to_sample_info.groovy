#!/usr/bin/env groovy 
// vim: ts=4:sw=4:expandtab:cindent
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// HaloPlex Batch Preparation Script
//
// This script performs a number of tasks to go from FASTQ files to a Cpipe batch ready for analysis.
//
// These include:
//
//  - Rename any samples where the sample name starts with a number - this prevents complications down the line
//  - Create a samples.txt file with sample id, fastq and sex inferred from the file names
//  - Create a samples.ped file with ids, family structure and sex (karyotype) inferred from file names
//  - Copy the variant database so that there is a dedicated copy for the batch
//  - Create a target_regions.txt file that sets the DSD target region, and the correct update and annotation
//    variant databases
// 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Author: Simon Sadedin, MCRI
// 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

System.err.println "=" * 100
System.err.println "HaloPlex DSD Batch PreparationTool"
System.err.println "=" * 100

msg = {
    System.err.println ""
    System.err.println((" "+it+" ").center(100,"="))
    System.err.println ""
}

Cli cli = new Cli(header: """Creates a samples.txt and optionally a samples.ped file from data 
                             files that contain sample information in their names
                          """)
cli.with {
    batch "batch to which samples belong", args:1, required: false
    design "target region design id", args:1, required: true
    mask "regex mask to apply to file names", args:1
    ped "generate PED file by parsing karyotype and relationships from sample names", args:1
    meta "name of sample meta data file to write", args:1, required: true
}

opts = cli.parse(args)
if(!opts)
        System.exit(0)

batch = opts.batch 

// If batch is not provided but we are run within a batch directory, use that
if(!batch) {
    if(new File(".").canonicalFile.parentFile.name == "batches") {
        batch = new File(".").canonicalFile.name
    }
    else {
        System.err.println "ERROR: please run this command inside your batch directory OR provide -batch"
    }
}


// For CRAM file support
System.properties.reference="/home/simon/work/hg19/gatk.ucsc.hg19.fasta"

// println "Reference = " + System.properties.reference

mask = opts.mask

fastqFiles = new File("data").listFiles().grep { it.name.endsWith("fastq.gz") }.collect { "data/" +  it.name }

if(!fastqFiles) {
    System.err.println """
        ERROR: No FASTQ files were found in the data directory. 
        
        Please copy the files to the data directory before running this script
    """
    System.exit(1)
}

samples = SampleInfo.fromFiles(fastqFiles.grep { !(new File(it).name.startsWith("Undetermined"))}, mask ?:null)

def r = System.in.newReader()

// Start by picking out any samples that start with numeric identifiers: we will rename these to start with X
samples.grep { it.key ==~ '^[0-9].*$' }.each { e ->
    def (id,sample) = [e.key, e.value]
    System.err.println "Sample $id starts with a numeric value. Do you want to rename to start with X? (y/n)"
    def answer = r.readLine().trim()
    if(answer == "y") {
        sample.files.fastq.each {  fastq ->
            def fastqFile = new File(fastq)
            def newName = new File(fastqFile.parentFile, "X" + fastqFile.name) 
            System.err.println "Rename: $fastq -> $newName.absolutePath" 
            fastqFile.renameTo(newName.absolutePath)
        }
    }
}

samples = SampleInfo.fromFiles(fastqFiles.grep { !(new File(it).name.startsWith("Undetermined"))}, mask ?:null)

fathers = []
mothers = []
children = []

samples.each {  id, s ->
    s.batch = batch
    s.target = opts.design 
    if(opts.ped) {
        def parts = (s.sample =~ '(-[MF]){0,1}-(4[0-9][XY]{1,5})_(S|ML)[0-9]{1,2}')
        if(parts) {
            s.sex = Sex.OTHER
            if(parts[0][2] == '46XY')
                s.sex = Sex.MALE
            else
            if(parts[0][2] == '46XX')
                s.sex = Sex.FEMALE

            if(parts[0][1] == '-F') 
                fathers << s
            else
            if(parts[0][1] == '-M') 
                mothers << s
            else
                children << s
        }
        else {
            System.err.println "WARNING: cannot parse info from sample id $s.sample"
        }
    }
}

sexes = ["UNKNOWN" : 0, "MALE" : 1, "FEMALE" : 2, "OTHER" : 0 ]

simplify_sample_id = { sampleId ->
	sampleId.replaceAll('(-4[0-9][XY][XY]).*$','$1')
}

if(opts.ped)  {
    msg("Creating PED file")

    System.err.println "Please check the following information is correct: \n"

    System.err.println "Fathers: " + fathers*.sample.join("\n")
    System.err.println "\nMothers: " + mothers*.sample.join("\n")
    System.err.println "\nChildren: " + children*.sample.join("\n")

    System.err.println ""

    new File(opts.ped).withWriter() { w ->
        for(child in children) {

            def family = child.sample.replaceAll('(-[MF]){0,1}-(4[0-9][XY]{1,5})_(S|ML)[0-9]{1,2}.*$','')
            def father = samples.find { it.key =~ (family + '-F-') }?.value
            def mother = samples.find { it.key =~ (family + '-M-') }?.value

            child.sample = simplify_sample_id(child.sample)

           if(father == null) {
                System.err.println "No father for $child.sample"
            }
            else {
                System.err.println "Father of $child.sample is $father.sample"
                father.sample = simplify_sample_id(father.sample)
                w.println([ family, father.sample, 0, 0, sexes[father.sex.name()], 1 ].join("\t"))
            }

            if(mother == null) {
                System.err.println "No mother for $child.sample"
            }
            else {
                System.err.println "Mother of $child.sample is $mother.sample"
                mother.sample = simplify_sample_id(mother.sample)
                w.println([ family, mother.sample, 0, 0, sexes[mother.sex.name()], 1 ].join("\t"))
            }

            w.println([ family, child.sample, father?.sample?:0, mother?.sample?:0, sexes[child.sex.name()], 2 ].join("\t"))
        }
    }
    println ""
}

msg "Write sample meta data file (${opts.meta})"
new File(opts.meta).withWriter { w ->
    w.println(
        samples*.value*.toTsv().join("\n")
    )
}

msg "Copying Variant Database ..." 

ant = new AntBuilder()
ant.saveStreams = false
ant.copy(file:"../../tools/haloplex_dsd/variants.dsd.db", todir:".")

msg("Writing target_regions.txt")
targetRegions = """
EXOME_TARGET="\$BASE/designs/$opts.design/${opts.design}.bed"
VARIANT_DB="\$BASE/tools/haloplex_dsd/variants.dsd.db"
UPDATE_VARIANT_DB="\$BASE/tools/haloplex_dsd/variants.dsd.db"
ANNOTATION_VARIANT_DB="\$BASE/tools/haloplex_dsd/variants.dsd.db"
"""

// Create a target_regions.txt file, if we are confident about where to do it
if(new File(".").canonicalFile.parentFile.name == "batches") {
    new File("target_regions.txt").text = targetRegions
}
else {
    System.err.println "WARNING: because you ran this command outside of the batch directory, a target_regions.txt file was not created. Please create this file and add the following contents: \n\n" + targetRegions + "\n"
}

msg("Creating analysis directory")
ant.mkdir(dir:"analysis")

if(opts.ped) {
    System.err.println """

    NOTE: it is assumed that all singletons / children are probands / affected. Please check 
          and correct this if it is not the case
    """
}

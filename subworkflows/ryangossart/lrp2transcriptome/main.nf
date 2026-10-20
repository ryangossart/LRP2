/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SQANTI_LRP2QC              } from '../../../modules/ryangossart/sqanti/lrp2qc/main'
include { GTF_LRP2GENERATEHASHIDS    } from '../../../modules/ryangossart/gtf/lrp2generatehashids/main'
include { GTF_LRP2FILTERTRANSCRIPTOME } from '../../../modules/ryangossart/gtf/lrp2filtertranscriptome/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW DEFINITION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow TRANSCRIPTOME {
    take:
    ch_samples                  // channel: [meta, [collapsed_gff, flnc_counts]]
    reference_gtf               // path: reference GENCODE GTF
    reference_fasta             // path: reference GENCODE FASTA
    sample_metadata             // path: sample metadata CSV

    main:
    ch_versions = channel.empty()

    //
    // MODULE: Run SQANTI3 QC on collapsed transcripts
    //
    SQANTI_LRP2QC (
        ch_samples,
        reference_gtf,
        reference_fasta
    )
    ch_versions = ch_versions.mix(SQANTI_LRP2QC.out.versions)

    //
    // MODULE: Generate hashids for all transcripts
    //
    GTF_LRP2GENERATEHASHIDS (
        SQANTI_LRP2QC.out.corrected_gtf
            .join(SQANTI_LRP2QC.out.classification, by: 0),
        reference_gtf
    )
    ch_versions = ch_versions.mix(GTF_LRP2GENERATEHASHIDS.out.versions)

    //
    // MODULE: Filter SQANTI3 output (requires hashids_mapping from GTF_LRP2GENERATEHASHIDS)
    //
    GTF_LRP2FILTERTRANSCRIPTOME (
        SQANTI_LRP2QC.out.classification
            .join(SQANTI_LRP2QC.out.corrected_gtf, by: 0)
            .join(SQANTI_LRP2QC.out.corrected_fasta, by: 0)
            .join(GTF_LRP2GENERATEHASHIDS.out.hashids_mapping, by: 0)
            .join(GTF_LRP2GENERATEHASHIDS.out.corrected_psl, by: 0)
            .map { meta, classification, gtf, fasta, hashids, psl ->
                [meta, classification, gtf, fasta, hashids, psl] },
        reference_gtf,
        reference_fasta,
        sample_metadata
    )
    ch_versions = ch_versions.mix(GTF_LRP2FILTERTRANSCRIPTOME.out.versions)

    emit:
    // SQANTI_LRP2QC outputs
    classification          = SQANTI_LRP2QC.out.classification                    // [meta, *_classification.txt]
    corrected_gtf           = SQANTI_LRP2QC.out.corrected_gtf                     // [meta, *_corrected.gtf]
    corrected_fasta         = SQANTI_LRP2QC.out.corrected_fasta                   // [meta, *_corrected.fasta]
//    corrected_genepred      = SQANTI_LRP2QC.out.corrected_genepred                // [meta, *_corrected.genePred]
//    corrected_cds_gff       = SQANTI_LRP2QC.out.corrected_cds_gff                 // [meta, *_corrected.gtf.cds.gff]
//    isoforms_gtf            = SQANTI_LRP2QC.out.isoforms_gtf                      // [meta, *.isoforms.gtf]
    junctions               = SQANTI_LRP2QC.out.junctions                         // [meta, *_junctions.txt]
//    params                  = SQANTI_LRP2QC.out.params                            // [meta, *.params.txt]
//    refannotation_genepred  = SQANTI_LRP2QC.out.refannotation_genepred            // [meta, refAnnotation_*.genePred]

    // GTF_LRP2GENERATEHASHIDS outputs
    hashids_mapping         = GTF_LRP2GENERATEHASHIDS.out.hashids_mapping               // [meta, *_transcriptome_hashids_mapping.txt]

    // GTF_LRP2FILTERTRANSCRIPTOME outputs
//    transcriptome_classification = GTF_LRP2FILTERTRANSCRIPTOME.out.classification           // [meta, *_transcriptome_classification.txt]
    classification_filtered = GTF_LRP2FILTERTRANSCRIPTOME.out.classification_filtered       // [meta, *_classification_filtered.txt]
//    transcriptome_corrected_gtf = GTF_LRP2FILTERTRANSCRIPTOME.out.corrected_gtf             // [meta, *_transcriptome_corrected.gtf]
    corrected_gtf_filtered  = GTF_LRP2FILTERTRANSCRIPTOME.out.corrected_gtf_filtered        // [meta, *_corrected_filtered.gtf]
    corrected_bed_filtered  = GTF_LRP2FILTERTRANSCRIPTOME.out.corrected_bed_filtered        // [meta, *_corrected_filtered.bed]
//    transcriptome_corrected_fasta = GTF_LRP2FILTERTRANSCRIPTOME.out.corrected_fasta         // [meta, *_transcriptome_corrected.fasta]
    corrected_fasta_filtered= GTF_LRP2FILTERTRANSCRIPTOME.out.corrected_fasta_filtered      // [meta, *_corrected_filtered.fasta]
//    corrected_psl           = GTF_LRP2FILTERTRANSCRIPTOME.out.corrected_psl                 // [meta, *_transcriptome_corrected.psl]
//    hashids_mapping         = GTF_LRP2FILTERTRANSCRIPTOME.out.hashids_mapping               // [meta, *_transcriptome_hashids_mapping.txt]
    hashids_filtered        = GTF_LRP2FILTERTRANSCRIPTOME.out.hashids_filtered              // [meta, *_hashids_with_cpm_filtered.txt]
    hashids_all             = GTF_LRP2FILTERTRANSCRIPTOME.out.hashids_all                   // [meta, *_all_hashids_with_cpm.txt]

    versions                = ch_versions.unique().collectFile(name: 'versions.yml')
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

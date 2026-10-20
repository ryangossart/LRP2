/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CPAT_CPAT                     } from '../../../modules/ryangossart/cpat/cpat/main'
include { CPAT_FILTER                   } from '../../../modules/ryangossart/cpat/filter/main'
include { SQANTI_PROTEIN                } from '../../../modules/ryangossart/sqanti/protein/main'
include { SQANTI_PROTEINCLASSIFICATION  } from '../../../modules/ryangossart/sqanti/proteinclassification/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW DEFINITION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow PREDICTED_PROTEOME {
    take:
    ch_samples                  // channel: [meta, [corrected_fasta, corrected_gtf, classification_filtered, hashids_filtered, hashids_mapping]]
    reference_gtf               // path: reference GENCODE GTF
    hexamer_file                // path: CPAT hexamer file (species-specific)
    logit_model                 // path: CPAT logit model (species-specific)

    main:
    ch_versions = channel.empty()

    //
    // MODULE: Predict ORFs using CPAT
    //
    CPAT_CPAT (
        ch_samples.map { meta, fasta, gtf, classification, hashids, mapping ->
            [meta, fasta] },
        hexamer_file,
        logit_model
    )
    ch_versions = ch_versions.mix(CPAT_CPAT.out.versions)

    //
    // MODULE: Filter CPAT output and call best ORFs
    //
    CPAT_FILTER (
        CPAT_CPAT.out.orf_prob
            .join(CPAT_CPAT.out.orf_seqs, by: 0)
            .join(ch_samples, by: 0)
            .map { meta, orf_prob, orf_seqs, fasta, gtf, _classification, hashids, mapping ->
                [meta, orf_prob, orf_seqs, fasta, gtf, mapping] },
        reference_gtf
    )
    ch_versions = ch_versions.mix(CPAT_FILTER.out.versions)

    //
    // MODULE: Run SQANTI3 protein classification
    //
    SQANTI_PROTEIN (
        CPAT_FILTER.out.cds_gtf,
        reference_gtf
    )
    ch_versions = ch_versions.mix(SQANTI_PROTEIN.out.versions)

    //
    // MODULE: Custom protein and 5'UTR classification
    //
    SQANTI_PROTEINCLASSIFICATION (
        SQANTI_PROTEIN.out.protein_classification
            .join(CPAT_FILTER.out.cds_gtf, by: 0)
            .join(ch_samples, by: 0)
            .join(CPAT_FILTER.out.all_orfs_mapped, by: 0)
            .map { meta, protein_class, cds_gtf, fasta, gtf, classification, hashids, mapping, all_orfs ->
                [meta, protein_class, cds_gtf, fasta, hashids, all_orfs] },
        reference_gtf
    )
    ch_versions = ch_versions.mix(SQANTI_PROTEINCLASSIFICATION.out.versions)

    emit:
    // CPAT outputs
    orf_prob                    = CPAT_CPAT.out.orf_prob                                     // [meta, *_cpat.ORF_prob.tsv]
    orf_prob_best               = CPAT_CPAT.out.orf_prob_best                                // [meta, *_cpat.ORF_prob.best.tsv]
    orf_seqs                    = CPAT_CPAT.out.orf_seqs                                     // [meta, *_cpat.ORF_seqs.fa]
    no_orf                      = CPAT_CPAT.out.no_orf                                       // [meta, *_cpat.no_ORF.txt]
    cpat_error                  = CPAT_CPAT.out.error_log                                    // [meta, *_cpat.error]

    // CPAT_FILTER outputs
    all_orfs_mapped             = CPAT_FILTER.out.all_orfs_mapped                           // [meta, *_all_orfs_mapped.tsv]
    cds_gtf                     = CPAT_FILTER.out.cds_gtf                                   // [meta, *_CDS.gtf]

    // SQANTI_PROTEIN outputs
    protein_classification      = SQANTI_PROTEIN.out.protein_classification  // [meta, *.sqanti_protein_classification.tsv]

    // SQANTI_PROTEINCLASSIFICATION outputs
//    protein_all_isoforms         = SQANTI_PROTEINCLASSIFICATION.out.protein_all_isoforms       // [meta, protein_sqanti/*_protein_classification_all_isoforms.txt]
    protein_all_orfs_fasta       = SQANTI_PROTEINCLASSIFICATION.out.protein_all_orfs_fasta     // [meta, protein_sqanti/*.predicted.proteome.all_best_orfs.fa]
    hashids_orf                  = SQANTI_PROTEINCLASSIFICATION.out.hashids_orf                // [meta, protein_sqanti/*.predicted.proteome.high_confidence_ORF_cpm.txt]
    protein_gtf                  = SQANTI_PROTEINCLASSIFICATION.out.protein_gtf                // [meta, protein_sqanti/*.predicted.proteome.high_confidence_ORF.gtf]
    protein_bed                  = SQANTI_PROTEINCLASSIFICATION.out.protein_bed                // [meta, protein_sqanti/*.predicted.proteome.high_confidence_ORF.bed]
    versions                     = ch_versions.unique().collectFile(name: 'versions.yml')
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

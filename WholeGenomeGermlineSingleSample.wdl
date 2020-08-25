version 1.0

## Copyright Broad Institute, 2018
##
## This WDL pipeline implements data pre-processing and initial variant calling (GVCF
## generation) according to the GATK Best Practices (June 2016) for germline SNP and
## Indel discovery in human whole-genome data.
##
## Requirements/expectations :
## - Human whole-genome pair-end sequencing data in unmapped BAM (uBAM) format
## - One or more read groups, one per uBAM file, all belonging to a single sample (SM)
## - Input uBAM files must additionally comply with the following requirements:
## - - filenames all have the same suffix (we use ".unmapped.bam")
## - - files must pass validation by ValidateSamFile
## - - reads are provided in query-sorted order
## - - all reads must have an RG tag
## - GVCF output names must end in ".g.vcf.gz"
## - Reference genome must be Hg38 with ALT contigs
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

#import "./tasks/UnmappedBamToAlignedBam.wdl" as ToBam
#import "./tasks/AggregatedBamQC.wdl" as AggregatedQC
#import "./tasks/Qc.wdl" as QC
#import "./tasks/BamToCram.wdl" as ToCram
#import "./tasks/VariantCalling.wdl" as ToGvcf
#import "./structs/GermlineStructs.wdl"

import "https://raw.githubusercontent.com/microsoft/gatk4-genome-processing-pipeline-azure/om_bwaSplitHardReduceNoMerge/tasks/UnmappedBamToAlignedBam.wdl" as ToBam
import "https://raw.githubusercontent.com/microsoft/gatk4-genome-processing-pipeline-azure/om_bwaSplitHardReduceNoMerge/tasks/AggregatedBamQC.wdl" as AggregatedQC
import "https://raw.githubusercontent.com/microsoft/gatk4-genome-processing-pipeline-azure/om_bwaSplitHardReduceNoMerge/tasks/Qc.wdl" as QC
import "https://raw.githubusercontent.com/microsoft/gatk4-genome-processing-pipeline-azure/om_bwaSplitHardReduceNoMerge/tasks/BamToCram.wdl" as ToCram
import "https://raw.githubusercontent.com/microsoft/gatk4-genome-processing-pipeline-azure/om_bwaSplitHardReduceNoMerge/tasks/VariantCalling.wdl" as ToGvcf
import "https://raw.githubusercontent.com/microsoft/gatk4-genome-processing-pipeline-azure/om_bwaSplitHardReduceNoMerge/structs/GermlineStructs.wdl"

# WORKFLOW DEFINITION
workflow WholeGenomeGermlineSingleSample {

  String pipeline_version = "1.4"

  input {
    SampleAndUnmappedBams sample_and_unmapped_bams
    GermlineSingleSampleReferences references
    PapiSettings papi_settings
    File wgs_coverage_interval_list

    File? haplotype_database_file
    Boolean provide_bam_output = false
    Boolean use_gatk3_haplotype_caller = false

    String bwa_commandline = "bwa mem -K 100000000 -p -v 3 -t 16 -Y $bash_ref_fasta"
    String bwa_cores = "16"
    String bwa_docker = "bwaoptimized.azurecr.io/genomes-in-the-cloud:2.4.3-1564508330-msopt"


    Float large_bam_in_gb = 20.0


  }

  # Not overridable:
  Int read_length = 250
  Float lod_threshold = -20.0
  String cross_check_fingerprints_by = "READGROUP"
  String recalibrated_bam_basename = sample_and_unmapped_bams.base_file_name + ".aligned.duplicates_marked.recalibrated"

  call ToBam.UnmappedBamToAlignedBam {
    input:
      sample_and_unmapped_bams    = sample_and_unmapped_bams,
      references                  = references,
      papi_settings               = papi_settings,

      contamination_sites_ud = references.contamination_sites_ud,
      contamination_sites_bed = references.contamination_sites_bed,
      contamination_sites_mu = references.contamination_sites_mu,

      cross_check_fingerprints_by = cross_check_fingerprints_by,
      haplotype_database_file     = haplotype_database_file,
      lod_threshold               = lod_threshold,
      recalibrated_bam_basename   = recalibrated_bam_basename,

      bwa_commandline = bwa_commandline,
      bwa_cores = bwa_cores,
      bwa_docker = bwa_docker,

      cutoff_for_large_rg_in_gb = large_bam_in_gb
  }


  call ToGvcf.VariantCalling as BamToGvcf {
    input:
      calling_interval_list = references.calling_interval_list,
      evaluation_interval_list = references.evaluation_interval_list,
      haplotype_scatter_count = references.haplotype_scatter_count,
      break_bands_at_multiples_of = references.break_bands_at_multiples_of,
#      contamination = UnmappedBamToAlignedBam.contamination,
      input_bam = UnmappedBamToAlignedBam.output_bam,
      input_bam_index = UnmappedBamToAlignedBam.output_bam_index,
      ref_fasta = references.reference_fasta.ref_fasta,
      ref_fasta_index = references.reference_fasta.ref_fasta_index,
      ref_dict = references.reference_fasta.ref_dict,
      dbsnp_vcf = references.dbsnp_vcf,
      dbsnp_vcf_index = references.dbsnp_vcf_index,
      base_file_name = sample_and_unmapped_bams.base_file_name,
      final_vcf_base_name = sample_and_unmapped_bams.final_gvcf_base_name,
      agg_preemptible_tries = papi_settings.agg_preemptible_tries,
      use_gatk3_haplotype_caller = use_gatk3_haplotype_caller
  }

  if (provide_bam_output) {
    File provided_output_bam = UnmappedBamToAlignedBam.output_bam
    File provided_output_bam_index = UnmappedBamToAlignedBam.output_bam_index
  }

  # Outputs that will be retained when execution is complete 
  output {
    File duplicate_metrics = UnmappedBamToAlignedBam.duplicate_metrics
    File output_bqsr_reports = UnmappedBamToAlignedBam.output_bqsr_reports

    File? output_bam = provided_output_bam
    File? output_bam_index = provided_output_bam_index

    File output_vcf = BamToGvcf.output_vcf
    File output_vcf_index = BamToGvcf.output_vcf_index
  }
}

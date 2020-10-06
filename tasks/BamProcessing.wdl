version 1.0

## Copyright Broad Institute, 2018
##
## This WDL defines tasks used for BAM file processing of human whole-genome or exome sequencing data.
##
## Runtime parameters are often optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

# Sort BAM file by coordinate order
task SortSam {
  input {
    File input_bam
    String output_bam_basename
    Int preemptible_tries
    Int compression_level
  }
  # SortSam spills to disk a lot more because we are only store 300000 records in RAM now because its faster for our data so it needs
  # more disk space.  Also it spills to disk in an uncompressed format so we need to account for that with a larger multiplier
  Float sort_sam_disk_multiplier = 3.25
  Int disk_size = ceil(sort_sam_disk_multiplier * size(input_bam, "GB")) + 20

  command {
    java -Dsamjdk.compression_level=~{compression_level} -Xms4000m -jar /usr/gitc/picard.jar \
      SortSam \
      INPUT=~{input_bam} \
      OUTPUT=~{output_bam_basename}.bam \
      SORT_ORDER="coordinate" \
      CREATE_INDEX=true \
      CREATE_MD5_FILE=true \
      MAX_RECORDS_IN_RAM=300000

  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564508330"
    disk: disk_size + " GB"
    cpu: "1"
    memory: "5000 MB"
    preemptible: true
    maxRetries: preemptible_tries
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
    File output_bam_index = "~{output_bam_basename}.bai"
    File output_bam_md5 = "~{output_bam_basename}.bam.md5"
  }
}

# Mark duplicate reads to avoid counting non-independent observations
task MarkDuplicates {
  input {
    Array[File] input_bams
    String output_bam_basename
    String metrics_filename
    Float total_input_size
    Int compression_level
    Int preemptible_tries

    # The program default for READ_NAME_REGEX is appropriate in nearly every case.
    # Sometimes we wish to supply "null" in order to turn off optical duplicate detection
    # This can be desirable if you don't mind the estimated library size being wrong and optical duplicate detection is taking >7 days and failing
    String? read_name_regex
    Int memory_multiplier = 1
    Int additional_disk = 20							
  }

  # The merged bam will be smaller than the sum of the parts so we need to account for the unmerged inputs and the merged output.
  # Mark Duplicates takes in as input readgroup bams and outputs a slightly smaller aggregated bam. Giving .25 as wiggleroom
  Float md_disk_multiplier = 3
  Int disk_size = ceil(md_disk_multiplier * total_input_size) + additional_disk

  Float memory_size = 7.5 * memory_multiplier
  Int java_memory_size = (ceil(memory_size) - 2)

  # Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly
  # This works because the output of BWA is query-grouped and therefore, so is the output of MergeBamAlignment.
  # While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"

  command {
    java -Dsamjdk.compression_level=~{compression_level} -Xms~{java_memory_size}g -jar /usr/gitc/picard.jar \
      MarkDuplicates \
      INPUT=~{sep=' INPUT=' input_bams} \
      OUTPUT=~{output_bam_basename}.bam \
      METRICS_FILE=~{metrics_filename} \
      VALIDATION_STRINGENCY=SILENT \
      ~{"READ_NAME_REGEX=" + read_name_regex} \
      OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
      ASSUME_SORT_ORDER="queryname" \
      CLEAR_DT="false" \
      ADD_PG_TAG_TO_READS=false
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564508330"
    preemptible: true
    maxRetries: preemptible_tries
    memory: "~{memory_size} GB"
    disk: disk_size + " GB"
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
    File duplicate_metrics = "~{metrics_filename}"
  }
}

# Generate Base Quality Score Recalibration (BQSR) model
task BaseRecalibrator {
  input {
    File input_bam
    File input_bam_index
    String recalibration_report_filename
    Array[String] sequence_group_interval
    File dbsnp_vcf
    File dbsnp_vcf_index
    Array[File] known_indels_sites_vcfs
    Array[File] known_indels_sites_indices
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Int bqsr_scatter
    Int preemptible_tries
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.0.10.1"
  }

  Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB")
  Float dbsnp_size = size(dbsnp_vcf, "GB")
  Int disk_size = ceil(size(input_bam, "GB") + size(input_bam_index, "GB") + ref_size + dbsnp_size) + 20

  parameter_meta {
    input_bam: {
      localization_optional: true
    }
  }

  command {
    gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
      -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
      -Xloggc:gc_log.log -Xms5g" \
      BaseRecalibrator \
      -R ~{ref_fasta} \
      -I ~{input_bam} \
      --use-original-qualities \
      -O ~{recalibration_report_filename} \
      --known-sites ~{dbsnp_vcf} \
      --known-sites ~{sep=" -known-sites " known_indels_sites_vcfs} \
      -L ~{sep=" -L " sequence_group_interval}
  }
  runtime {
    docker: gatk_docker
    preemptible: true
    maxRetries: preemptible_tries
    memory: "6 GB"
    disk: disk_size + " GB"
  }
  output {
    File recalibration_report = "~{recalibration_report_filename}"
  }
}

# Apply Base Quality Score Recalibration (BQSR) model
task ApplyBQSR {
  input {
    File input_bam
    File input_bam_index
    String output_bam_basename
    File recalibration_report
    Array[String] sequence_group_interval
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Int compression_level
    Int bqsr_scatter
    Int preemptible_tries
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.0.10.1"
    Int memory_multiplier = 1
    Int additional_disk = 20							
  }

  Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB")
  Int disk_size = ceil(size(input_bam, "GB") * 3 + size(input_bam_index, "GB") * 3 + ref_size) + additional_disk

  Int memory_size = ceil(3500 * memory_multiplier)

  parameter_meta {
    input_bam: {
      localization_optional: true
    }
  }

  command {
    gatk --java-options "-XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps \
      -XX:+PrintGCDetails -Xloggc:gc_log.log \
      -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Dsamjdk.compression_level=~{compression_level} -Xms3000m" \
      ApplyBQSR \
      --create-output-bam-md5 \
      --add-output-sam-program-record \
      -R ~{ref_fasta} \
      -I ~{input_bam} \
      --use-original-qualities \
      -O ~{output_bam_basename}.bam \
      -bqsr ~{recalibration_report} \
      --static-quantized-quals 10 \
      --static-quantized-quals 20 \
      --static-quantized-quals 30 \
      -L ~{sep=" -L " sequence_group_interval}
  }
  runtime {
    docker: gatk_docker
    preemptible: true
    maxRetries: preemptible_tries
    memory: "~{memory_size} MB"
    disk: disk_size + " GB"
  }
  output {
    File recalibrated_bam = "~{output_bam_basename}.bam"
    File recalibrated_bam_checksum = "~{output_bam_basename}.bam.md5"
  }
}

# Combine multiple recalibration tables from scattered BaseRecalibrator runs
task GatherBqsrReports {
  input {
    Array[File] input_bqsr_reports
    String output_report_filename
    Int preemptible_tries
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.0.10.1"
  }

  command {
    gatk --java-options "-Xms3000m" \
      GatherBQSRReports \
      -I ~{sep=' -I ' input_bqsr_reports} \
      -O ~{output_report_filename}
    }
  runtime {
    docker: gatk_docker
    preemptible: true
    maxRetries: preemptible_tries
    memory: "3500 MB"
    disk: "20 GB"
  }
  output {
    File output_bqsr_report = "~{output_report_filename}"
  }
}

# Combine multiple *sorted* BAM files
task GatherSortedBamFiles {
  input {
    Array[File] input_bams
    String output_bam_basename
    Float total_input_size
    Int compression_level
    Int preemptible_tries
  }

  # Multiply the input bam size by two to account for the input and output
  Int disk_size = ceil(2 * total_input_size) + 20

  command {
    java -Dsamjdk.compression_level=~{compression_level} -Xms2000m -jar /usr/gitc/picard.jar \
      GatherBamFiles \
      INPUT=~{sep=' INPUT=' input_bams} \
      OUTPUT=~{output_bam_basename}.bam \
      CREATE_INDEX=true \
      CREATE_MD5_FILE=true
    }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564508330"
    preemptible: true
    maxRetries: preemptible_tries
    memory: "3 GB"
    disk: disk_size + " GB"
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
    File output_bam_index = "~{output_bam_basename}.bai"
    File output_bam_md5 = "~{output_bam_basename}.bam.md5"
  }
}

# Combine multiple *unsorted* BAM files
# Note that if/when WDL supports optional outputs, we should merge this task with the sorted version
task GatherUnsortedBamFiles {
  input {
    Array[File] input_bams
    String output_bam_basename
    Float total_input_size
    Int compression_level
    Int preemptible_tries
  }

  # Multiply the input bam size by two to account for the input and output
  Int disk_size = ceil(2 * total_input_size) + 20

  command {
    java -Dsamjdk.compression_level=~{compression_level} -Xms2000m -jar /usr/gitc/picard.jar \
      GatherBamFiles \
      INPUT=~{sep=' INPUT=' input_bams} \
      OUTPUT=~{output_bam_basename}.bam \
      CREATE_INDEX=false \
      CREATE_MD5_FILE=false
    }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564508330"
    preemptible: true
    maxRetries: preemptible_tries
    memory: "3 GB"
    disk: disk_size + " GB"
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
  }
}


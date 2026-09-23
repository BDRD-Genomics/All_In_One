#!/bin/bash

fastq_dir=$1

rsync -av ${fastq_dir}/*.fastq.gz /export/sequencing_dir/

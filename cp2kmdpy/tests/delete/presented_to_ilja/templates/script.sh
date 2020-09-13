#!/bin/sh
{% extends base_script %}
{% block project_header %}
#PBS -l walltime=00:10:00,mem=30gb,nodes=1:ppn=1

module purge
module load mkl
module load fftw
module load intel/cluster/2018

date >> execution.log
{{ super() }}
{% endblock %}

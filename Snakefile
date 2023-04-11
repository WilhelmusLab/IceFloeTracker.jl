configfile: "./hpc/snakemake-config.yaml"
from snakemake.utils import validate
from os import *

envvars:
  "SPACEUSER",
  "SPACEPSWD"

validate(config, "./hpc/snakemake-config.yaml") ## requires a schema

rule all:
  input: "runall.txt", "soit.txt", "preprocess.txt"

rule fetchdata:
  output: parent = directory(config["fetchdata-outdir"]), t = touch("runall.txt"), truedir = directory(config["truecolor-outdir"]), refdir = directory(config["reflectance-outdir"])      
  params:
    help = "-h",
    start = config["startdate"],
    end = config["enddate"],
    bb = config["bounding-box"],
  shell: "./scripts/fetchdata.sh -o {output.parent} -s {params.start} -e {params.end} {params.bb}"

rule soit:
  input: rules.fetchdata.output.t
  output: outdir = directory(config["soit-outdir"]), s = touch("soit.txt")
  shell: """
          mkdir -p {output.outdir}
          python3 ./scripts/pass_time_snakemake.py
         """
# run delta_time script

rule landmask:
  input: rules.soit.output.s
  params: lmdir = rules.fetchdata.output.parent
  output: outdir = directory(config["landmask-outdir"]), outfile = "output/landmasks/generated_landmask.jls"
  shell: """
          mkdir -p {output.outdir}
          ./scripts/ice-floe-tracker.jl landmask {params.lmdir} {output.outdir}
         """

rule preprocess: 
  input: rules.landmask.output.outfile
  output: outdir = directory(config["preprocess-outdir"]), p = touch("preprocess.txt")
  shell: """
          mkdir -p {output.outdir}
          julia -t auto ./scripts/ice-floe-tracker.jl preprocess -t {rules.fetchdata.output.truedir} -r {rules.fetchdata.output.refdir} -l {rules.landmask.output.outdir} -o {output.outdir}
         """

rule extractfeatures:
  input: rules.preprocess.output.p
  params: minarea = config["minfloearea"], maxarea = config["maxfloearea"]
  output: directory(config["features-outdir"])
  shell: """
          mkdir -p {output}
          julia -t auto ./scripts/ice-floe-tracker.jl extractfeatures -i {rules.preprocess.output.outdir} -o {output} --minarea {params.minarea} --maxarea {params.maxarea}
         """

rule cleanup:
  input: rules.preprocess.output.p
  shell: """
          rm runall.txt
          rm soit.txt
          rm preprocess.txt
         
         """
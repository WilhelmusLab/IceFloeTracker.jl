configfile: "./snakemake-config.yaml"
from snakemake.utils import validate
from os import *

# make wildcard for truecolor images

envvars:
  "SPACEUSER",
  "SPACEPSWD"

validate(config, "./snakemake-config.yaml") ## requires a schema

rule all:
  input: "file.done"
  
#   input: 
#     directory(config["landmask-outdir"])
#lmdir = directory(config["landmask-outdir"]),

rule mkdir:
  output:
          lmdir = directory(config["landmask-outdir"]),
          lmgen = directory(config["lm-generated"]),
          soitdir = directory(config["soit-outdir"]),
  shell: """
          mkdir -p {output.soitdir}
          mkdir -p {output.lmdir}
          mkdir -p {output.lmgen}
         """

rule fetchdata:
  output: parent = directory(config["fetchdata-outdir"]), t = touch("file.done")      
  params:
    help = "-h",
    start = config["startdate"],
    end = config["enddate"],
    bb = config["bounding-box"],
  shell: "./scripts/fetchdata.sh -o {output.parent} -s {params.start} -e {params.end} {params.bb}"

rule soit:
  output: directory(config["soit-outdir"])
  shell: "python3 ./scripts/pass_time_snakemake.py"
# run delta_time script

# rule mvlm:
#   params: infile = "images/landmask.tiff",
#           outdir = directory("landmasks")
#   shell: "mv {params.infile} {params.outdir}"

rule landmask:
  input: "file.done"
  params: lmdir = "landmasks",
          outdir = "landmasks/generated"
  shell: """
          mv images/landmask.tiff landmasks
          ./scripts/ice-floe-tracker.jl landmask {params.lmdir} {params.outdir}
        """

# rule preprocess:

# rule segmentation:

# rule feature_extraction:

#output: directory(rules.mkdir.output.soitdir)
#output: directory(config["landmask-outdir"])
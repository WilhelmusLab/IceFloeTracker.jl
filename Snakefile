configfile: "./snakemake-config.yaml"
from snakemake.utils import validate

envvars:
  "SPACEUSER",
  "SPACEPSWD"

validate(config, "./snakemake-config.yaml") ## requires a schema

# rule all:
#   input: 
#     directory(config["landmask-outdir"])

rule mkdir:
  params: lmdir = directory(config["landmask-outdir"])
  shell: "mkdir -p {params.lmdir}"

rule mkex:
  shell: "chmod a+x ./scripts/ice-floe-tracker.jl"

rule fetchdata:
  output: directory(config["fetchdata-outdir"])
  params:
    help = "-h",
    start = config["startdate"],
    end = config["enddate"],
    bb = config["bounding-box"],
  shell: "./scripts/fetchdata.sh -o {output} -s {params.start} -e {params.end} {params.bb}"

# rule soit:
# #   input:
# #   output:
#   params:
#     spacetrackuser=os.environ["SPACEUSER"],
#     spacetrackpassword=os.environ["SPACEPSWD"]
  #shell: 
#         run soit 
#          run delta_time script

rule landmask:
  input: {rules.fetchdata.output}
  params: lmdir = {rules.fetchdata.output}
  script: "./scripts/ice-floe-tracker.jl landmask {input} {params.lmdir}"

# rule preprocess:

# rule segmentation:

# rule feature_extraction:
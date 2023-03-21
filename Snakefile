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
  output: main = directory(config["main-output"]),
          lmdir = directory(config["landmask-outdir"]),
          soitdir = directory(config["soit-outdir"])
  shell: """
          mkdir -p {output.main}
          mkdir -p {output.lmdir}
          mkdir -p {output.soitdir}
         """

rule fetchdata:
  output: parent = directory("output/data")      
  params:
    help = "-h",
    start = config["startdate"],
    end = config["enddate"],
    bb = config["bounding-box"],
  shell: "./scripts/fetchdata.sh -o {output.parent} -s {params.start} -e {params.end} {params.bb}"

rule soit:
  output: directory(rules.mkdir.output.soitdir)
  shell: "python3 ./scripts/pass_time_snakemake.py"
#          run delta_time script

rule landmask:
  #input: "output/data/truecolor"
  params: indir = "output/data/truecolor"
  output: directory(rules.mkdir.output.lmdir)
  script: "./scripts/ice-floe-tracker.jl landmask {params.indir} {output}"

# rule preprocess:

# rule segmentation:

# rule feature_extraction:
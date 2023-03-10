from datetime import datetime
current_time = datetime.now().strftime("%Y%m%d%H%M%S")
from pathlib import Path

# rule all:
#   input: "data/output/landmask", "data/output/cloudmask"

rule fetchdata:
  output:
    out = directory("data")
  params:
    help = "-h",
    startdate = "2023-02-08",
    enddate = "-e 2023-02-14",
    outdir = lambda wildcards, output: Path(output.out),
    bb = "81 -22 79 12",
    
  shell: "./scripts/fetchdata.sh -o {output}_{current_time} -s {params.startdate} {params.enddate} {params.bb}"

# rule soit:
#   input:
#   output:
#   shell: run soit 
#          run delta_time script

# rule landmask:
#   input: 
#     metadata="data/input/metadata.json",
#     images="data/input"
#   output: directory("data/output/landmask")
#   shell: "./scripts/ice-floe-tracker.jl landmask {input.metadata} {input.images} {output}"

# rule cloudmask:
#   input: 
#     metadata="data/input/metadata.json",
#     images="data/input/reflectance"
#   output: directory("data/output/cloudmask")
#   shell: "./scripts/ice-floe-tracker.jl cloudmask {input.metadata} {input.images} {output}"


# rule preprocess:

# rule segmentation:

# rule feature_extraction:
rule all:
  input: "data/output/landmask", "data/output/cloudmask"

rule fetchdata:
  output: 
    directory("data/input"),
    "data/input/landmask.tiff", 
    directory("data/input/truecolor"),
    directory("data/input/reflectance"),
    "data/input/metadata.json"
  shell: "./scripts/ice-floe-tracker.jl fetchdata {output[0]}"

rule landmask:
  input: 
    metadata="data/input/metadata.json",
    images="data/input"
  output: directory("data/output/landmask")
  shell: "./scripts/ice-floe-tracker.jl landmask {input.metadata} {input.images} {output}"

rule cloudmask:
  input: 
    metadata="data/input/metadata.json",
    images="data/input/reflectance"
  output: directory("data/output/cloudmask")
  shell: "./scripts/ice-floe-tracker.jl cloudmask {input.metadata} {input.images} {output}"

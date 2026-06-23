# Workflow Management

IceFloeTracker.jl can be used with workflow management tools for batch processing.

- The example in this directory uses [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html).
- The main file specifying the workflow is the [snakefile](./snakefile). 
- The [config.yaml](./config.yaml) file sets environment variables.

## Install dependencies

Install the workflow and its dependencies by calling:
```
# Install the snakemake runner in an isolated environment
pipx install snakemake  

# Add pyproj and pandas to the snakemake environment
pipx inject snakemake pyproj pandas  
```

## Run the workflow

To run the workflow, call:
```bash
snakemake -c 1
```

This will use one processor core (`-c 1`) 
to run the workflow to calculate all the tracking results for the cases in [case.csv](./case.csv)
which uses the region definitions from [region.csv](./region.csv). 
The results will be written to the `./results` directory which is not checked into the repository. 

(This command is equivalent to `snakemake -c 1 region_case_results_by_filetype` 
which is the default target of the snakemake workflow.)

You can also specify the output files which are desired, 
and Snakemake should download or create any (missing) prerequisites 
and process them in the correct order.

For example, the command:
```bash
snakemake track/beaufort_sea-100km.250m.2019-03-22.2019-03-23.LopezAcosta2019Tiling.tracked.csv 
```
- Will use the number of cores from the [default profile](./profiles/default/config.yaml)
- to run tracking (`track/...tracked.csv`)
- on a 100km x 100km region in the Beaufort Sea, 
- using 250m scale images, 
- from 22nd to 23rd March 2019,
- with the `Lopez2019Tiling.Segment` algorithm.

It reports that 29 tasks of different types are planned before executing the first:
```
Assuming unrestricted shared filesystem usage.
host: c486485cb3cb
Building DAG of jobs...
Using shell: /usr/bin/bash
Provided cores: 1 (use --cores to define parallelism)
Rules claiming more threads will be scaled down.
Job stats:
job                                count
-------------------------------  -------
get_region_month_overpass_times        1
get_single_overpass_time               4
load_case_coastal_buffer_mask          4
load_case_landmask                     4
load_coastal_buffer_mask               1
load_falsecolor                        4
load_landmask                          1
load_raw_landmask                      1
load_truecolor                         4
segment_lopez_tiling                   4
tracking                               1
total                                 29
```

Intermediate files can be made in the same way, e.g.:
```bash
snakemake -c 4 scene/hudson_bay-1500km.250m.2023-03-{22..25}.terra/falsecolor.tiff
```
- will use four cores (`-c 4`)
- to load falsecolor images (`.../falsecolor.tiff`)
- of the 1500km x 1500km image of Hudson Bay
- at the 250m scale
- for the 22nd through 25th March 2023.

To keep going if a single function in the pipeline fails, 
use the `--keep-going` flag to continue to continue to run any independent tasks.
However, any tasks which depend on failed jobs will still fail.

### Configuration

By default, the [default configuration file](./configs/default/config.yaml) will be used. 

It specifies:
- the [default case list](./configs/default/case.csv)
- the [default region list](./configs/default/region.csv)
- which Julia environment is used to run IceFloeTracker.jl commands
- and other commands used in the workflow.

If an additional configuration is specified on the command line using `--configfile`
it will be combined with the default configuration
and can override variables from the default configuration.

Example configurations are included:
- [with regions and timepoints in a single CSV file](./configs/case-region-file/config.yaml), 
  call
  ```bash
  snakemake --configfile workflow/configs/case-region-file/config.yaml
  ```
- [a large validation dataset](https://github.com/danielmwatkins/ice-floe-validation-dataset/blob/main/data/validation_dataset/validation_dataset.csv) 
  call
  ```bash
  snakemake --configfile workflow/configs/validation-cases/config.yaml
  ```
- [use the version of IceFloeTracker.jl from this repo](./configs/ift-from-this-repo/config.yaml)
  ```bash
  snakemake --configfile workflow/configs/ift-from-this-repo/config.yaml
  ```

> [!TIP]
> Multiple config files can be used at once, e.g.:
> ```bash
> snakemake --configfile workflow/configs/ift-from-this-repo/config.yaml --configfile workflow/configs/validation-cases/config.yaml
> ```

### Files produced in the workflow

Files produced in the workflow include:

- Observation specific files, organized by observation directory `scene/{region}.{scale}m.{date_}.{satellite}/`, e.g. `beaufort_sea.250m.2019-03-23.aqua/`:
  - `truecolor.tiff` – the input truecolor image (MODIS bands 143),
  - `falsecolor.tiff` – the input falsecolor image (MODIS bands 721),
  - `cloud.tiff` – the input cloud mask from the MODIS observations,
  - `landmask.tiff` – the landmask for the region,
  - `coastal_buffer_mask.tiff` – the dilated landmask,
- Observation and pipeline results files, organized by observation and pipeline directory `scene/{region}.{scale}m.{date_}.{satellite}/{pipeline}/`, e.g. `scene/beaufort_sea.250m.2019-03-23.aqua/LopezAcosta2019/`:
  - `segmentation.hdf5` – the combined results file with segmentation results, floe masks, floe properties etc.,
  - `labels_map.tiff` – the floe labels from the segmentation,
  - `segment_mean_truecolor.tiff` – each region in `labels_map.tiff` with the mean color of that region from `truecolor.tiff`,
  - `segment_mean_falsecolor.tiff` – each region in `labels_map.tiff` with the mean color of that region from `falsecolor.tiff`,
  - `cloud_mask.tiff` – regions where cloud is detected in the image,
  - `ice_mask.tiff` – regions where ice is detected in the image,
  - `config.txt` – listing of the parameters used in the segmentation,


## Satellite Overpass Identification Tool Concurrency Limit

The Satellite Overpass Identification Tool depends on the rate limits of space-track.org. 
As of March 6, 2026, the rate limits on API calls are 30 per minute _and_ 300 per hour.[^1]

A conservative global resource limit for `soit_api_calls` 
is set in [profiles/default/config.yaml](./profiles/default/config.yaml).
The two rules `get_region_overpass_times` and `get_region_month_overpass_times` each make two API calls,
so a resource limit of `soit_api_calls=4` means that only two instances of those rules can run concurrently.

You can specify a different limit in the snakemake call:
```bash
snakemake <other arguments> --resources soit_api_calls=16
```

[^1]: https://www.space-track.org/documentation

## Specifying new regions

Regions like `hudson_bay-1500km` can be specified by adding them to the [region.csv](./region.csv) file.
Each region is identified by a name, and specified by:
- its center in a particular Coordinate Reference System (CRS), usually [EPSG:4326](https://epsg.io/4326) (latitude and longitude in decimal degrees),
- the target (output) CRS, usually [EPSG:3413](https://epsg.io/3413) (NSIDC Sea Ice Polar Stereographic North)
- and its extent in the target CRS, which for EPSG:3413 is in metres. 

## Batch processing

For batch processing, [case.csv](./case.csv) specifies which regions and which cases will be run in the batch. 
Each row specifies: 
- the region name (from [region.csv](./region.csv)), 
- the start and end dates (in ISO 8601 YYYY-MM-DD format)
- the image scale in metres
- the `pipeline` value, corresponding to the Julia pipeline/module used by the workflow.

The supported batch processing rules are:
- `region_case_results`: organize the segmentation results from `case.csv` in directories named `{region}.{scale}.{date}.{satellite}`,
- `region_case_results_by_filetype`: reorganize the segmentation results from `case.csv` in directories by `{filetype}`. 

The files are organized as follows:

|Type|`region_case_results` path|`region_case_results_by_filetype` path|
|----|--------|---|
|Input file template|`scene/{scene}/{filetype}.{extension}`|`filetype/{filetype}/{scene}.{extension}`|
|Input file example|`scene/beaufort_sea.250m.2019-03-23.aqua/falsecolor.tiff`|`filetype/falsecolor/beaufort_sea.250m.2019-03-23.aqua.tiff`|
|Output file template|`scene/{scene}/{pipeline}/{filetype}.{extension}`|`filetype/{pipeline}/{filetype}/{scene}.{extension}`|
|Output file example|`scene/beaufort_sea.250m.2019-03-23.aqua/LopezAcosta2019/segment_mean_falsecolor.tiff`|`filetype/LopezAcosta2019/segment_mean_falsecolor/beaufort_sea.250m.2019-03-23.aqua.tiff`|


To invoke the batch processing rule, call `snakemake <other arguments> <batch_processing_rule_name>`

In all cases, the tracking results are stored in the top-level of the `results` directory, 
named like `{region}.{scale}.{start date}.{end date}.{pipeline}.tracked.csv`.

### `region_case_results`

Organizes the segmentation results from `case.csv` in directories named 
- `scene/{region}.{scale}.{date}.{satellite}` for input files, and
- `scene/{region}.{scale}.{date}.{satellite}/{pipeline}` for output files. 

Example invocation using 4 processing cores:
```shell
snakemake -c 4 region_case_results
```

### `region_case_results_by_filetype`

Reorganize the segmentation results from `case.csv` in directories by 
- `{filetype}` for input files, and
- `{pipeline}/{filetype}` for output files. 

Example invocation using 4 processing cores:
```shell
snakemake -c 4 region_case_results_by_filetype
```

The reorganized files are hard-linked (i.e., not copies), so they don't take up any additional space. 

## Anatomy of a Snakemake Rule

The workflow relies on being able to call the IceFloeTracker.jl from the command line.
In this workflow example, we make each rule based on a Julia function, a one-line or very short script.
A non-trivial example is shown here:

```
rule segment_lopez:
    output:
        labels_map = "{dir}/lopez.labels_map.tiff",
        cloud_mask = "{dir}/lopez/cloud_mask.tiff",
    input:
        truecolor = "{dir}/truecolor.tiff",
        falsecolor = "{dir}/falsecolor.tiff",
        landmask = "{dir}/landmask.tiff",
        coastal_buffer_mask = "{dir}/coastal_buffer_mask.tiff",
    shell:
        """
        {config[IFT]} -e 'using IceFloeTracker, Images; LopezAcosta2019.Segment()(
            load("{input.truecolor}"),
            load("{input.falsecolor}"),
            load("{input.landmask}") |> binarize_mask,
            load("{input.coastal_buffer_mask}") |> binarize_mask;
            intermediate_results_callback=call_kwargs(;
                labels_map=l -> l .|> UInt16 |> save("{output.labels_map}"),
                cloud_mask=save("{output.cloud_mask}"),
            ),
        )'
        """
```

It is comprised of the following parts:

- Output and input filenames containing a wildcard `{dir}` – these are where all the outputs generated by the script will be written.
    ```
    rule segment_lopez:
        output:
            labels_map = "{dir}/lopez.labels_map.tiff",
            cloud_mask = "{dir}/lopez/cloud_mask.tiff",
        input:
            truecolor = "{dir}/truecolor.tiff",
            falsecolor = "{dir}/falsecolor.tiff",
            landmask = "{dir}/landmask.tiff",
            coastal_buffer_mask = "{dir}/coastal_buffer_mask.tiff",
        ...
    ```
- Shell command:
    ```
    rule segment_lopez:
        ...
        shell:
            """
            {config[IFT]} -e 'using IceFloeTracker, Images; 
            LopezAcosta2019.Segment(
                diffusion_algorithm = PeronaMalikDiffusion()
            )(
                load("{input.truecolor}"),
                load("{input.falsecolor}"),
                load("{input.landmask}") |> binarize_mask,
                load("{input.coastal_buffer_mask}") |> binarize_mask;
                intermediate_results_callback=call_kwargs(;
                    cloud_mask=save("{output.cloud_mask}"),
                    labels_map=l -> l .|> UInt16 |> save("{output.labels_map}"),
                ),
            )'
            """
    ```
    with the following subcomponents:
    - `{config[IFT]}`: the name of the Julia environment, like `"julia --project=/path/to/instantiated/IceFloeTracker.jl` configured in the [config.yaml](config.yaml) file.
    - `-e` flag to call julia in "eval" mode, which will evaluate the following string as a (series of) Julia command(s).
    - `'using IceFloeTracker...'` the Julia command to run.

Each of the high-level algorithms defined by IceFloeTracker, like `LopezAcosta2019.Segment(;kwargs...)(images...)`, 
has a structure which supports its evaluation on the command line. 

Each function is a "functor", which accepts keyword arguments to define how it behaves, 
like `.Segment(diffusion_algorithm = PeronaMalikDiffusion())(...` 
which sets the diffusion algorithm used in the preprocessing step.

The instantiated functor can be called,
and will accept a series of arguments to run the actual calculation.
The arguments are often images which have to be loaded. 
- Images which need no conversion can be loaded simply like `load("{input.truecolor}")`,
- Images which need conversion, for instance converting to a binary mask, 
  can be piped to a converting function, e.g. `load("{input.landmask}") |> binarize_mask`

For saving outputs, the high-level algorithms accept a callback function in the   
`intermediate_results_callback` keyword argument. 
The `call_kwargs(; kwargs...)` function passes each matching result, 
like `labels_map` (the overall output of a segmentation algorithm),
or `cloud_mask` (an intermediate result),
into the function specified. 
In this example, two functions are evaluated:
- `cloud_mask |> save("{output.cloud_mask}")`, 
  which saves the cloud mask 
  to the file specified in the output list, 
  `"{dir}/lopez/cloud_mask.tiff"`
- and `labels_map .|> UInt16 |> save("{output.labels_map}")`, 
  which converts the array of labels 
  from the default 64-bit integer to an unsigned 16-bit integer, 
  and then saves the result to the file specified in the output list,
  `"{dir}/lopez.labels_map.tiff"`


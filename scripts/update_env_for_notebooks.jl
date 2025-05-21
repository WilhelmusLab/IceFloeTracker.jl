using TOML
using Pkg

# Parse the Project.toml file
project = TOML.parsefile("Project.toml")

# Get the list of packages under the notebooks target
notebooks_targets = get(project["targets"], "notebooks", String[])

# Add each package
Pkg.add.(notebooks_targets);

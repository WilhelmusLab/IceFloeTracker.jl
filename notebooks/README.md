# Example Notebooks

## Recommened Settings - VS Code
1. Install Julia
2. Install VS Code
3. Add the following extensions for VS Code: Julia, Jupyter
4. (Optional) Ensure Julia is set up to work with threading enabled for shared memory multiprocessing.
   a. Press Cmd+Shift+P (mac) or Ctrl+Shift+P (windows) to display the command palette and search for `settings.json` and choose the option that appears
   
    ![Alt text](imgs/settings.png)
    
   b. In the `settings.json` file, find the ` "julia.NumThreads"` option (near the bottom). Change its value to `"auto"` and save the file
   
   ![Alt text](imgs/set-threads-auto.png)

5. Open the `./notebooks/preprocessing-workflow/preprocessing-workflow.ipynb` file within VS Code

## Recommended Settings - Jupyter/JupyterLab
1. Install Julia
2. Install [miniconda](https://www.anaconda.com/docs/getting-started/miniconda/install)
3. Make a new miniconda environment with Jupyter and activate it:
   ```@bash
   conda create --name ift-env jupyter
   conda activate ift-env
   ```
4. Open Julia in the terminal
5. Activate a Julia environment
   ```@julia
   using Pkg
   Pkg.activate("ift-demo")
   ```
6. Install `IJulia` and `IceFloeTracker`
   ```@julia
   Pkg.add("IJulia")
   Pkg.add("IceFloeTracker")
   ```
7. Install optional packages if desired
   ```@julia
   Pkg.add("Plots")
   ```
8. Resolve and instantiate:
   ```@julia
   Pkg.resolve()
   Pkg.instantiate()
   ```
9. You can then exit Julia via `exit()`, and open Jupyter lab:
   ```
   jupyter lab
   ```
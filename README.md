# Miscellaneous scripts
A collection of (mostly Julia) scripts I've been at least a bit of effort into. All scripts are contained within their own subdirectory and exist as their own Julia project (where applicable).

Generally, `Manifest.toml` files **are** included, to allow for reproducibility and slightly quicker plug-and-play. If you want the most up-to-date versions of plugins, make sure you delete these files before running the corresponding script.
## Contents
### [alphaMover](alphaMover/)
The `alphaMover.jl` script extracts AlphaFold2-generated top ranked model files (`.pdb` and `.pkl`) from multi-prediction subdirectories. This can be useful when batch predicting many models and when you only want to keep the best model and its metadata.

Requires single argument, an AlphaFold2 output directory (e.g. `outputDir`).

```sh
julia ./alphaMover/alphaMover.jl outputDir
```

### [pickledJSON](pickledJSON/)

> [!CAUTION]
> Pickling is not secure, so make sure you trust the origin of your pickle files. [Better explained here](https://docs.python.org/3/library/pickle.html).

The `pickledJSON.jl` script converts AlphaFold2-generated pickle files into JSON, keeping only pTM, pLDDT, PAE, and max PAE metrics. This can be useful when needing to minimise storage space of a lot of files, without needing to keep all the metadata contained within the pickle files. Output JSON files are compatible with the [PAE Viewer webserver](https://pae-viewer.uni-goettingen.de/).

Requires single argument, a target directory of pickle files (e.g. `pklDir`).

```sh
julia ./pickledJSON/pickledJSON.jl pklDir
```

**Note:** I'm yet to add ipTM to the list of metrics to keep, so for now this works for monomer predictions only.


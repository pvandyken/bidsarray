# bidsarray

Compile a tabular output of input and output files from a bids dataset. Output can be fed into [gnu parallel](https://www.gnu.org/software/parallel/sphinx.html) to run arbitrary commands on an entire bids dataset with parallelization!


## Command syntax

### Prelude

The first part of a `bidsarray` call consists of prelude arguments that initialize the input and output dirs:

```bash
bidsarray INPUT_DIR OUTPUT_DIR [--derivatives] \
    [--participant-label LABEL ...] [--exclude-participant-label LABEL ...] \
    [--pybidsdb-dir DIR] [--pybidsdb-reset]
```

`--derivatives` enables indexing of derivative datasets (in the `derivatives/` folder). `--participant-label` allows the specification of one or more subject labels (the bit after sub: `sub-LABEL`). Only files from these subjects will be produced. `--exclude-participant-label` does the same but excludes subjects. `--pybidsdb-dir` can be used to specify a pybids database created via the `database_path` argument in `bids.BIDSLayout`. This can speed up indexing for large datasets. `--pybidsdb-reset` forces reindexing of the database.

### Components

Components are specified after the prelude, each seperated by `:::` (just like in gnu parallel!). So a complete command looks like this:

```bash
bidsarray <PRELUDE...> ::: <COMPONENT 1> ::: <COMPONENT 2> ::: ...
```

Each component may be specified as an input or an output. Inputs are read from an existing dataset, outputs are created based on the provided inputs.

#### Inputs

Input components are specified as follows:

```bash
::: --input [label] [--groupby ENTITY ...] [--aggregate ENTITY ...] [--filter ENTITY[:METHOD]=VALUE ...]
```

The `label` is optional: if provided, it will add a header row to the top of the tabular output. If a label is provided to any one component, all components must receive a label!

`--filter` narrows the set of selected paths to those containing the selected entity-values. For example, to select diffusion images, one might use:

```bash
::: --input --filter suffix=dwi datatype=dwi extension=.nii.gz
```

`METHOD` tells the filter how to do the selection. By default, it looks for an exact match, but regex can also be used using `:match` and `:search`:

```bash
::: --input --filter 'suffix:match=[Tt][12]w?'
```

Overall, `--filter` tends to REDUCE the number of rows in the output.

`--groupby` and `--aggregate` both mark variable parts of the path. Often, they'll be entities such as `subject`, `session`, and `run`. These variable entities are referred to generically as wildcards.

`--groupby` will create a seperate row per path matched by the filters and wildcards. This is used if you want to run a tool seperately on each file. For instance, if you have a group of images, and you want to apply a smoothing function to each one seperately, you may use:

```bash
::: --input --filter suffix=T1w extension=.nii.gz --groupby subject session
```

`--aggregate` joins multiple paths into the same row. This is useful if you want to perform an aggregation, such as an average or standard deviation. For example, if you want to get an average fractional anisotropy map from all subjects, you may use:

```bash
::: --input --filter desc=FA suffix=mdp extension=.nii.gz --aggregate subject session
```

These two flags can be combined. So to get an average FA map for each subject across all sessions, you could use:

```bash
::: --input --filter desc=FA suffix=mdp extension=.nii.gz --groupby subject --aggregate session
```

#### Outputs

Outputs are generated using the bids function from [`snakebids`](https://snakebids.readthedocs.io/en/stable/api/paths.html#snakebids.bids). The syntax is:

```bash
::: --output [label] --entities [ENTITY=VALUE ...]
```

Each `ENTITY_VALUE` pair provided to `--entities` specifies a static value that should be applied to all wildcards. For instance, continuing our smoothing example, you may use:

```bash
::: --input --filter suffix=T1w extension=.nii.gz --groupby subject session ::: --output --entities datatype=anat suffix=T1w desc=smoothed
```

Importantly, any wildcards specified using `--groupby` are AUTOMATICALLY provided to each output. So in the above example, our table will include a row for each subject-session combination, with a correctly formatted output path.

## Feeding outputs to `parallel`

### Basics

While the tabular output could be used for any number of purposes, it really shines in combination with gnu parallel. Parallel is a powerful and complicated tool; its full usage can be found on its [documentation](https://www.gnu.org/software/parallel/sphinx.html). Here we show just the basics using it with bidsarray.

```bash
bidsarray <PRELUDE...> ::: <COMPONENT 1> ::: <COMPONENT 2> | parallel --colsep '\t' echo 1={1} 2={2}
```

The `--colsep` argument to `parallel` allows it to read from an incoming columnar data. We use `\t` as the argument to read `bidsarray` tabular data.

### Examples

In this example, we apply a transform to all of the fractional anisotropy maps of a preprocessed diffusion dataset:

```bash
bidsarray . derivatives/template --derivatives \
    ::: --input --filter suffix=mdp desc=FA extension=.nii.gz space=participant --groupby subject session \
    ::: --input --filter suffix=xfm from=participant to=MNI6 extension=.nii.gz --groupby subject session \
    ::: --input --filter suffix=T1w extension=.nii.gz space=MNI6 \
    ::: --output --entities desc=FA datatype=dwi space=MNI6 suffix=mdp.nii.gz |
    parallel --bar --colsep '\t' antsApplyTransform -d3 -i {1} -o {4} -r {3} -t {2}
```

### Using labels

The above examples use numeric ids for argument substitution, which may be difficult when handling many components. It's possible to use labels instead:

```bash
bidsarray . derivatives/template --derivatives \
    ::: --input image --filter suffix=mdp desc=FA extension=.nii.gz space=participant --groupby subject session \
    ::: --input transform --filter suffix=xfm from=participant to=MNI6 extension=.nii.gz --groupby subject session \
    ::: --input reference --filter suffix=T1w extension=.nii.gz space=MNI6 \
    ::: --output out --entities desc=FA datatype=dwi space=MNI6 suffix=mdp.nii.gz |
    parallel --bar --colsep '\t' --header : antsApplyTransform -d3 -i {image} -o {out} -r {reference} -t {transform}
```

Note that the use of labels in `parallel` is not compatible with all features.

### Aggregation commands with `parallel`

Parallel automatically shell escapes each column when using `--colsep`, including spaces, but `bidsarray` uses spaces to seperate files that should be aggregated. So if we try to get the average of set of files, `parallel` will read the the filenames as one giant filename. So we need to use a trick to get this to work. There are two basic approaches (both complements of this [SO Q/A](https://superuser.com/questions/1135733/gnu-parallel-remove-escape-before-space-characters-in-command)):

The simplest is to prepend the command with `eval` to remove all escapes:

```bash
bidsarray . derivatives/average \
    ::: --input --filter suffix=T1w extension=.nii.gz --groupby subject session --aggregate run \
    ::: --output --entities suffix=T1w.nii.gz datatype=anat \
    parallel --bar --colsep '\t' eval mrcalc {1} -add {2}
```

This may not work with complex commands, as eval will aggressively strip away quotes. For a more surgical approach, you can use the `uq` function:

```bash
bidsarray . derivatives/average \
    ::: --input --filter suffix=T1w extension=.nii.gz --groupby subject session --aggregate run \
    ::: --output --entities suffix=T1w.nii.gz datatype=anat \
    parallel --bar --colsep '\t' mrcalc {=1 uq=} -mean {2}
```

`uq` will only apply to the targeted variable. This approach is not compatible with labels.

### Creating output folders

Note that for most commands (especially any involving `--groupby`), you will likely need to create the ouput folders as part of the command. Use a command like this:

```bash
bidsarray . derivatives/template --derivatives \
    ::: --input --filter suffix=mdp desc=FA extension=.nii.gz space=participant --groupby subject session \
    ::: --input --filter suffix=xfm from=participant to=MNI6 extension=.nii.gz --groupby subject session \
    ::: --input --filter suffix=T1w extension=.nii.gz space=MNI6 \
    ::: --output --entities desc=FA datatype=dwi space=MNI6 suffix=mdp.nii.gz |
    parallel --bar --colsep '\t' mkdir -p {4} \&\& antsApplyTransform -d3 -i {1} -o {4} -r {3} -t {2}
```


# Motivation and Design

`bidsarray` was built both to be a useful tool and a demonstration of the `bidsapp` module in [`snakebids`](https://github.com/khanlab/snakebids). The parsing, components, and CLI are all provided by `snakebids`, so `bidsarray` can organize a tabular output with just a few small files.

`snakebids.bidsapp` uses a system of hooks and plugins to build an app. In <bidsarray/run.py>, three hooks can be seen:

* `get_argv`: Retrieve the provided CLI arguments and split them at `:::`. This will allow the app to seperately parse the command prelude and each of the components. The prelude arguments are returned to be parsed by the bidsapp, and the components are saved into the `config` to be parsed later.
* `finalize_config`: The prelude has now been parsed by bidsapp. We retrieve the components arguments we saved earlier and parse them (using functions in <bidsarray/component.py>). The results are saved back into `config`
* `run`: We use configured components we calculated earlier and call `generate_inputs`, which uses pybids to parse the input dataset and create a [`BidsDataset`](https://snakebids.readthedocs.io/en/stable/api/structures.html#snakebids.BidsDataset). The methods on this dataset can be used to retrieve and organize the indexed paths.

[`snakebids.bidsapp.app`](https://snakebids.readthedocs.io/en/stable/api/app.html#snakebids.bidsapp.app) is used to create the app with several plugins, providing the basic functionality. The last plugin: `sys.modules[__name__]`, loads the hooks defined in the `run.py` file so that the app will work. Finally, the entrypoint `app.run()` is called behind an `if __name__ == "__main__"` block. It's also specified in the `pyproject.toml` file as a script.

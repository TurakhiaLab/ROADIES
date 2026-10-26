# Input and Output

This directory converts user-facing files and command-line arguments into the validated structures consumed by `MlipperSession`, and serializes final trees and placements.

## Read This Folder in This Order

1. `parse_file.hpp/.cpp`: the small in-memory alignment and model structures.
2. `input_validation.hpp/.cpp`: reusable validation and user-facing errors.
3. `CLI.hpp/.cpp`: workflow-specific command-line parsing and config assembly.
4. `tree_newick.hpp/.cpp` and `jplace.hpp/.cpp`: output formatting.

## Data Structures

- `parse::Alignment`: parallel `names` and `sequences` vectors plus `sites`.
- `parse::ModelConfig`: supported substitution model, rates, frequencies, Gamma categories, alpha, invariant proportion, and scaling mode.
- `cli::SmallTipConfig`: all validated inputs/options for placement or commit.
- `cli::DivideAndConquerConfig`: all validated inputs/options for D&C.

These config objects own parsed values. They do not own CUDA state.

```text
files + CLI text -> parsed values -> validation -> workflow config
                                                       |
                                                       v
                                                 MlipperSession
```

Nothing in this directory should allocate a tree on the GPU.

## Key Functions

- `parse::read_alignment_file()`: performs the [files-to-parsed-values step](#data-structures) and returns an `Alignment`.
- `cli::selectWorkflow()`: detect `--divide-and-conquer`; `main.cpp` uses this before invoking either full parser.
- `cli::loadSmallTipConfigFromCommandLine()`: builds the small-tip [validated workflow config](#data-structures): parse flags, load the tree/alignment/model, preprocess repeated sites, validate combinations, and return `SmallTipParseResult`.
- `cli::loadDivideAndConquerConfigFromCommandLine()`: perform the corresponding D&C setup without accepting a backbone tree. It retains the original full alignment for DIPPER; `main.cpp` preprocesses repeated sites after DIPPER returns the starting tree.
- `input::normalize_cli_path()` and `validate_output_path()` implement the [error contract](#error-contract) for paths: resolve paths relative to the invocation directory and validate output destinations.
- `input::validate_alignment_names()`, `validate_alignment_symbols()`, and `validate_query_reference_name_overlap()`: enforce alignment contracts before GPU initialization.
- `input::validate_model_inputs()`: reject unsupported model shapes rather than silently approximating them.
- `treeio::write_tree_to_newick_string()` / `write_tree_to_newick_file()` serialize the CPU source of truth, [`TreeBuildResult`](../tree/README.md#treebuildresult).
- `jplaceio::write_jplace()`: write ranked placements and edge annotations in EPA-ng-compatible jplace form.

## Error Contract

Use `input::CliError` for invalid user input. `ValidationError` prefixes an option name so `main.cpp` can report a concise failure. Parsing should finish before constructing expensive GPU state.

## Adding a CLI Option

1. Add storage to the appropriate temporary CLI arguments/config structure.
2. Register the option in only the workflow(s) that support it.
3. Validate ranges and mutually exclusive flags in `CLI.cpp`.
4. Copy the validated value into `SmallTipConfig` or `DivideAndConquerConfig`.
5. Consume it from `main.cpp` or a workflow-level session option, not from a low-level kernel.
6. Update help text, the repository README, and `tests/input_validation_test.cpp`.

Tree parsing itself is in `tree/tree_generation.cpp` because it constructs the algorithm's `TreeBuildResult`; this folder owns serialization and CLI I/O.

# Currently planned feature implementations/fixes

## Metaboverse

#### Bugs
- Metaboverse building still tries to prematurely read `.mvrs` output file before it has been fully loaded. Need to address sync/async issues.

#### Features
- scRNA-seq implmentation (#90) -- likely just need a walkthrough demonstrating how to use this effectively
- Interactive data prep and name mapping check (#74, #86)
    - see https://datatables.net/examples/data_sources/js_array.html 
    - Cross-check labels with MetaboAnalyst API and Metaboverse naming dictionaries
    - Show which ones likely won't map (a way to modify output `.mvrs` directly?)

## metaboverse-cli
#### Bugs
- Update tests broken by previous updates

#### Features
- Option to allow users to display user-provided labels (#87)
    - New node key (`user_label`)
    - During viz, if `user_label` is not None, use that 
- Flux metabolomics integration (#21)
- Give user option to choose inferred complex values and stats or vanilla
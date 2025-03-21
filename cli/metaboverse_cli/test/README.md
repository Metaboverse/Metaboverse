# Metaboverse Redox Pair Mapping Tests

This directory contains unit tests for the redox pair mapping functionality in Metaboverse.

## Overview

The tests verify that the redox pair mapping functionality correctly handles metabolites that exist in different oxidation states, ensuring they are treated as distinct entities during the mapping process. This is particularly important for metabolomics studies focusing on redox biology, energy metabolism, and oxidative stress.

## Tests Included

1. **Test gather_synonyms function**:
   - Verifies that synonyms are not shared between oxidized and reduced forms of metabolites
   - Tests NAD+/NADH, NADP+/NADPH, FAD/FADH2, and GSH/GSSG pairs

2. **Test map_attributes function**:
   - Verifies that redox pairs are correctly identified during the mapping process
   - Checks that each member of a redox pair is mapped to its correct value

3. **Test REDOX_PAIRS dictionary**:
   - Verifies that all expected redox pairs are defined in the dictionary
   - Checks that there is no overlap between synonyms of different redox pairs

## Redox Pairs Covered

The tests cover the following redox pairs:
- NAD+/NADH
- NADP+/NADPH
- FAD/FADH2
- GSH/GSSG (Glutathione)
- Thioredoxin (oxidized/reduced)
- Coenzyme Q/Ubiquinol
- Cytochrome C (Fe3+/Fe2+)
- Ferredoxin (oxidized/reduced)
- Ascorbic acid/Dehydroascorbic acid
- Lipoic acid/Dihydrolipoic acid
- FMN/FMNH2
- Methylene-THF/THF

## Running the Tests

To run the tests, navigate to the `cli/metaboverse_cli/test` directory and run:

```bash
python run_tests.py
```

Alternatively, you can run the tests using the unittest module directly:

```bash
python -m unittest unit_tests.py
```

## Expected Output

If all tests pass, you should see output similar to:

```
test_gather_synonyms_fad (__main__.TestRedoxPairMapping) ... ok
test_gather_synonyms_glutathione (__main__.TestRedoxPairMapping) ... ok
test_gather_synonyms_nad (__main__.TestRedoxPairMapping) ... ok
test_gather_synonyms_nadp (__main__.TestRedoxPairMapping) ... ok
test_map_attributes_redox_pairs (__main__.TestRedoxPairMapping) ... ok
test_redox_pairs_completeness (__main__.TestRedoxPairMapping) ... ok
test_redox_pairs_no_overlap (__main__.TestRedoxPairMapping) ... ok

=== Test Summary ===
Ran 7 tests
Failures: 0
Errors: 0
```

## Adding New Tests

To add tests for additional redox pairs:

1. Add the pair to the `REDOX_PAIRS` dictionary in `cli/metaboverse_cli/analyze/model.py`
2. Add test cases for the new pair in `unit_tests.py`
3. Run the tests to verify that the new pair is correctly handled 
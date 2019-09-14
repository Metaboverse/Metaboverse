# Curation
- Compile total analyte dictionary (chebi + uniprot + etc)
id: {
  id: analyte_id,
  name: analyte_name,
  type: ensemble, chebi, etc,
  other_info: ...
}

# Viz
- Fix coloring issues
- Include IntAct database?
- Prototype expression overlay with some missing data points
  - Compile all data types and normalize scoring together?
  - Have individual scales and color node edge different color
  - Layout organization
  - Add genes if analyte in complex database
  - Make elastic layout with D3?
  - for calculating difference, per reaction, find min and max product and reactant and see if more than an x fold change and return that as a hit if so
    - How to calculate differences
      - `for n, nbrsdict in G.adjacency():`
      - will give you neighboring pairs
      - or
      - cycle through reactions and collect min/max from reactants and products then calculate largest difference between reactants and products
  - on Hover, show name, p value if available
  - auto-scale network based on connectivity and number of hubs
  - figure out how to handle participating complexes and genes in a concise manner
  - figure out when something should have HSA vs ALL
  - would be nice if nodes stick to where they are replaced so user gets a starting configuration, but can then be re-arranged until its formatted nicely and then they can export. Edges stay the same
  - reactome already has a tool that does something similar -- how is this better?
    - time analysis
    - ID breakpoints
    - looks better and more organized
    - no stats developed
    - multi-omics?
    - in-house or piped automated normalization
    - all analysis is pathway level for hits, nothing at individual points
    - no predictive model for other omics measurements

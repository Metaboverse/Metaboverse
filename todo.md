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
- Prototype expression overlay with some missing data points
  - Compile all data types and normalize scoring together?
  - Have individual scales and color node edge different color
  - Layout organization
  - Add genes if analyte in complex database
  - Make elastic layout with D3?
  - for calculating difference, per reaction, find min and max product and reactant and see if more than an x fold change and return that as a hit if so
  - on Hover, show name, p value if available
  - auto-scale network based on connectivity and number of hubs
  - figure out how to handle participating complexes and genes in a concise manner
  - figure out when something should have HSA vs ALL

###############
Updates
###############

=================================
v0.5.0-beta (in progress)
=================================
| **Major**
| - Addresses `issue #66 <https://github.com/Metaboverse/Metaboverse/issues/66>`_ , by hosting curated reference and template files for each organism per version of Metaboverse. Also provided user argument options to directly include already downloaded or curated files (:data:`--organism_curation_file`, :data:`--neighbor_dictionary_file`, :data:`--graph_template_file`). Using pre-downloaded files, this reduced processing time for curating data on the human network from ~30-40 min to ~2 min. These pre-curated files will be processed with each release of Metaboverse and are hosted on SourceForge currently.
|
| **Minor**
| - Fixes an issue where an empty unmapped dataframe would cause an error (fixed by  `#32e9283 <https://github.com/Metaboverse/metaboverse-cli/commit/32e9283363bb9ce8c4ef2325184ad01d102f4680>`_ )
| - Fixes an issue the working path would be appended to the organism ID (fixed by  `#91a490d <https://github.com/Metaboverse/metaboverse-cli/commit/91a490dec409c7a27d1b2cc0207ded5dd0fa60c1>`_ )
| - Addresses `issue #67 <https://github.com/Metaboverse/Metaboverse/issues/67>`_ , where experiment name inputs with spaces would cause an error.
| - Bump required version of Electron to >=9.4.0 (see `pull request #68 <https://github.com/Metaboverse/Metaboverse/pull/68>`_ ).
| - Removed some unused user arguments from command-line interface.
| - Fixed an issue where the backend argument parse checker would try to append a file path to the organism ID.
| - Fixes internal warning for UI when CLI did not output blocklist or labels.
| - Updated copyright info.
| - Removed deprecated files.
| - Migrated from Travis-CI to GitHub Actions.

=================================
v0.4.0-beta
=================================
| **Major**
| - Partial collapse: Addresses  `#51 <https://github.com/Metaboverse/Metaboverse/issues/51>`_  , which introduces partial collapsing to the reaction collapsing utility within Metaboverse. Previously, perfect matches were required between two reactions to collapse the reactions. However, this can be overly stringent in key metabolic pathways where a metabolite that is output by one reaction may not be required for the subsequent reaction (perhaps ATP is produced by reaction A but is not required for reaction B). To perform a partial collapse, Metaboverse operates by largely the same scheme as before, but now checks for a perfect match from each neighboring reaction, and if a perfect match is not available, checks for partial matches by filtering out high-degree nodes (quartile 98 of all non-reaction node degrees) and then checking if at least 30% of the nodes match with its neighbor.
| - Improvements to nearest neighbor searches where all iterations of a species are included in the graphing. Previously, it would only use the literally selected node to search for neighbors, but Reactome provides separate species IDs for a metabolite's different organelle-localizations, which was complicating these searches.

| **Minor**
| - Displays a preview of the user-selected reaction in an interactive format on the Pattern Search Analysis page until the user selects a Pathway to visualize. If a reaction is collapsed and spans two pathways, no pathways will be shown and instead a note that the reaction spans two pathways is displayed.
| - The Pattern Search Analysis page now allows users to filter out collapsed reactions from the search results. By default, collapsed reactions will be displayed until the checkbox is unchecked by the user.
| - Minor updates to Pattern Search Analysis page to make better usage of blank space
| - Fixes  `#60 <https://github.com/Metaboverse/Metaboverse/issues/60>`_  , where the :data:`.mvrs` file extension would not be automatically added to the user-provided output file name in Linux.
| - Addresses  `#62 <https://github.com/Metaboverse/Metaboverse/issues/62>`_  , where the some time-course/multi-condition slider bars would be improperly formatted.

=================================
v0.3.3-beta
=================================
| **Minor**
| - Closes `#63 <https://github.com/Metaboverse/Metaboverse/issues/63>`_ by applying :data:`safestr()`` function to all user input encodings to make sure no errors arise.

=================================
v0.3.2-beta
=================================
| **Minor**
| - Closes  `#59 <https://github.com/Metaboverse/Metaboverse/issues/59>`_  where non-ascii characters in reaction names would break the info extraction. Added a safestring conversion utility to prevent ascii-character issues.

=================================
v0.3.1-beta
=================================
| **Minor**
| - Fixes path separator for motif page name identification to allow for including modifiers in motif ID and exclusion of hubs ( `#55 <https://github.com/Metaboverse/Metaboverse/issues/55>`_ )
| - Fixes CHEBI mapping so that CHEBI IDs provided as input data are more reliably used as mapping IDs if it cannot match the metabolite by name ( `#58 <https://github.com/Metaboverse/Metaboverse/issues/58>`_ )
| - Fixes issue that arose in :data:`v0.3.0b` where some motif stamps could not be clicked on for viewing for timecourse/multi-condition data where it could not identify the shape for an unknown component type ( `#54 <https://github.com/Metaboverse/Metaboverse/issues/54>`_ )
| - Addresses  `#59 <https://github.com/Metaboverse/Metaboverse/issues/59>`_  where non-ascii characters in reaction names would break the info extraction. Was not able to recapitulate the error, but this fix, where relevant reaction metadata is forced to a string data-type, should allow for some flexibility here.
| - Updates walkthroughs and documentation to address ( `#31 <https://github.com/Metaboverse/Metaboverse/issues/31>`_ ) and update formatting

===========
v0.3.0-beta
===========
| **Major**
| - Allows for more flexible gene/protein mapping with Reactome-formatted node names. For example, Reactome will label a gene or protein with its isomer coordinates. Metaboverse now ignores those coordinates during attribute mapping of the user's data.
| - New naming of modules: :data:`Motif Search` is now called :data:`Pattern Analyis`, :data:`Visualize` is now called :data:`Explore`, and :data:`Connectivity` is now called :data:`Perturbation Networks`. Changed to be more descriptive and accessible to all users from broader backgrounds

| **Minor**
| - Fixed nearest neighbors capabilities in Perturbation Network visualization. In a previous release, a change had interfered with its function.
| - Use of outdated version will now direct user to the download page for the most recent version
| - The :data:`Back` button from any of the analysis modules will now redirect back to the index page
| - Removed compartment visualization from the :data:`Perturbation Networks` page as these often would just clutter the visualization and would not actually be helpful
| - Fixed reaction node formatting to turn off motif symbols/highlighting when timepoint/condition changed as this had been disabled by a previous change
| - Fixed collapsed reactions to ensure they were included in all reaction and motif formatting events
| - General formatting changes
| - Updated documentation

===========
v0.2.0-beta
===========
| **Major**
| - Fixes issues with missing metabolites during the network mapping stage ( `#37 <https://github.com/Metaboverse/Metaboverse/issues/37>`_ ). Addressed by re-working the metabolite synonym scheme to provide the same coverage of synonym look-up as before, but with more robustness so that some that were difficult to map would start mapping (i.e., Fructose)
| - Added dynamic line-plots of a selected reaction motif for time-course and multi-condition data ( `#15 <https://github.com/Metaboverse/Metaboverse/issues/15>`_). When exploring motifs on the Motif page for time-course and multi-condition experiments, a new panel appears at the bottom of the page which, for a selected motif, will show those reaction motif's component's behavior across all the time-points or conditions.
| - Added the option to exclude reaction motifs at a given time-point or condition that appear in another selected time-point or condition. ( `#16 <https://github.com/Metaboverse/Metaboverse/issues/16>`_ ). For example, if a user has selected to view motifs at a terminal time-point, but they want to know which reactions are motifs at this time-point but not at the initial time-point, they can exclude the motifs that show up at both time-points using the appropriate drop-down menu on the motif page.

| **Minor**
| - Metaboverse now outputs a table of unmapped metabolites ( `#35 <https://github.com/Metaboverse/Metaboverse/issues/35>`_ ).
| - Exploration pages now have pop-out bubbles with all information for compartments and node/link types `7d17d34 <https://github.com/Metaboverse/Metaboverse/commit/7d17d34aca5e900c307e266a07b4d82bd19a222d>`_.
| - Metaboverse new remembers and provides session info for experiment name, experiment type, labels, etc. and automatically fills those out for the user if returning to a page within the session `172d21a <https://github.com/Metaboverse/Metaboverse/commit/172d21a719bbc855fd46d4d8da223140c512a18f>`_.
| - Updated minor page formatting to make display more stable between Windows/Linux/Mac `52a100d <https://github.com/Metaboverse/Metaboverse/commit/52a100da0958af75c489165bc2f7c9eaf80294e8>`_.
| - Added test cases to CI for new/updated features
| - Updated package dependency information
| - Updated docs and FAQs

===========
v0.1.4-beta
===========
| - Fixes `#26 <https://github.com/Metaboverse/Metaboverse/issues/26>`_, where an error log is output if build fails
| - Removes direct Matplotlib imports in metaboverse-cli modules to prevent unnecessary errors and incompatibilities

===========
v0.1.3-beta
===========
| - Fixes bug where user paths with spaces were unable to be used ( `#26 <https://github.com/Metaboverse/Metaboverse/issues/26>`_ )

===========
v0.1.2-beta
===========
| - Fixes bug that prevented the curation from running without a blocklist ( `#19 <https://github.com/Metaboverse/Metaboverse/issues/19>`_ )
| - Fixes bug during data mapping that caused protein or gene values to occasionally map to metabolites ( `#20 <https://github.com/Metaboverse/Metaboverse/issues/20>`_ )

===========
v0.1.1-beta
===========
| - Fixes minor run-time issues with the Metaboverse interactive app
| - Fixes version alert to let users know if there is a newer version of Metaboverse available

===========
v0.1.0-beta
===========
| Initial beta release

===========
v0.0.1-beta
===========
| Demo pre-release with included human network data file for network visualization and exploration. Currently only available for MacOS.
|
| How to run:
|
| 1. Download attached :data:`.zip` demo file.
| 2. Double-click on :data:`.zip` file to uncompress Metaboverse and the accompanying test file
| 3. Within the uncompressed folder, right-click on Metaboverse to launch the app
| 4. Drag and drop the :data:`.json` file to the appropriate load icon and click the Visualize button.
|

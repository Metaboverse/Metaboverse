###############
Updates
###############

.. note::
    This page is no longer being updated. Please refer to the release notes for each release on the Metaboverse GitHub page: https://github.com/Metaboverse/Metaboverse/releases

=================================
v0.11.3
=================================
| - Fixes issues accessing pre-built files on Sourceforge and changing to Reactome database file keys (by @j-berg in `#159 <https://github.com/Metaboverse/Metaboverse/issues/159>`_)
| - Removes unnecessary warning on app load/exit (by @j-berg in `#160 <https://github.com/Metaboverse/Metaboverse/issues/160>`_)

=================================
v0.11.2
=================================
| - Entity mapping improvements so user can select between standard display names and user-provided names for metabolites (by @j-berg in `#153 <https://github.com/Metaboverse/Metaboverse/issues/153>`_)
| - Fixes to error on window close (by @j-berg in `c3f153a <https://github.com/Metaboverse/Metaboverse/commit/c3f153aec78d354ace178cf4b8bf6403a6bc9c60>`_)

=================================
v0.11.1
=================================
| - Fixes issues with HMDB ID mapping in back-end (`#144 <https://github.com/Metaboverse/Metaboverse/issues/144>`_, `#146 <https://github.com/Metaboverse/Metaboverse/issues/146>`_)
| - Prints note about using HMDB IDs with MetaboAnalyst API (`#144 <https://github.com/Metaboverse/Metaboverse/issues/144>`_)
| - Updates MetaboAnalyst API address for metabolite name conversion (`#144 <https://github.com/Metaboverse/Metaboverse/issues/144>`_)
| - Reverts to using Reactome metabolite names for visualization (will give the option in the future - `#151 <https://github.com/Metaboverse/Metaboverse/issues/151>`_, but need to handle cases where HMDB IDs are used as these are not very useful for visualization)
| - Updated documentation specifying that Ensembl transcript IDs or names should be used as Ensembl gene IDs do not work well. (`#146 <https://github.com/Metaboverse/Metaboverse/issues/146>`_)

=================================
v0.11.0
=================================
| - `#132 <https://github.com/Metaboverse/Metaboverse/issues/132>`_: Fixes issue with p-value/FDR calculation in interactive datatable formatting module
| - `#134 <https://github.com/Metaboverse/Metaboverse/issues/134>`_, `#136 <https://github.com/Metaboverse/Metaboverse/issues/136>`_: Moves curated files to Sourceforge, adds automated scripts for release building
| - `#128 <https://github.com/Metaboverse/Metaboverse/issues/128>`_, `#129 <https://github.com/Metaboverse/Metaboverse/issues/129>`_, `#130 <https://github.com/Metaboverse/Metaboverse/issues/130>`_, `#131 <https://github.com/Metaboverse/Metaboverse/issues/131>`_: Returns more information for values not being mapped to networks
| - Updates to front-end calls to improve performance

=================================
v0.10.1
=================================
| - Updated license to MIT
| - Fixed version check when launching GUI to make sure it sorts version numbers correctly
| - Updated documentation (expanded walkthrough, added dev notes)
| - Added more test datasets distributed with each release 
| - Added sample distribution display for Data Formatting tool 
| - Fixes issue where collapsed reactions in other timepoints were not being highlighted

=================================
v0.10.0
=================================
| - Introduces **Co-factor Selection**, a drop-down menu available in Pattern Analysis. Selected a specific factor will limit the Pattern Analysis results to only the reactions containing the selected entity.
| - Implements a Pattern Analysis **Table Export** utility, whereby the user can export the Reaction Patterns of a given type.
| - Fixes an issue where some non-collapsed reactions were not showing up in Pattern results.

=================================
v0.9.0
=================================
.. note::
    Files curated and analyzed using Metaboverse v0.9.0 will not be backwards compatible with files generated using earlier versions of Metaboverse.

| **Major**
| - Integration of confidence values during data formatting, curation, and analysis.
| - Fixed issues during curation arising from updates to the formatting of Reactome source files.
|
| **Minor**
| - Streamlined option buttons available during :data:`Pattern Analysis`.
| - Force directed networks now have sticky nodes once a user a dragged them to a new position.
| - Added statistical thresholding visualization to :data:`Pattern Analysis` module.
| - Output session file with each successful build - will have experiment name appended to file.
| - Sort :data:`Average` reaction pattern by best statistic on each side of reaction
| - Fix reaction filtering to more strictly and flexibly remove transport reactions from reaction pattern categories, and to only display transport reactions in the :data:`TransReg` reaction pattern category.
| - Node labels will default to the user-provided names if not the same as the curated default name.

=================================
v0.8.0
=================================
| **Major**
| - Added data formatting aid page (`see issue #86 <https://github.com/Metaboverse/Metaboverse/issues/86>`_).
| - Added ability to cross-reference metabolite names with MetaboAnalyst and find more compatible names (`see issue #74 <https://github.com/Metaboverse/Metaboverse/issues/74>`_).
|
| **Minor**
| - Fixed issue where Metaboverse would not populate selected organism from drop-down if previously selected and the user returns to the page (`see commit <https://github.com/Metaboverse/Metaboverse/commit/80d6ba995a71a1306d490cda768b2ed16174cf2a>`_).
| - Fixed issue where Metaboverse would not show the :data:`Continue` button on the :data:`Curation` page if necessary inputs were previously selected and the user returns to the page (`see commit <https://github.com/Metaboverse/Metaboverse/commit/80d6ba995a71a1306d490cda768b2ed16174cf2a>`_).


=================================
v0.7.1
=================================
| **Important Note**
| Many of the changes introduced in :data:`v0.7.0` to session and intermediate file metadata will likely be incompatible with previous versions of Metaboverse.
|
| **Minor**
| - Fix Session Data page to format variables, file paths better (`see commit <https://github.com/Metaboverse/Metaboverse/commit/07962e2a5d70a47a8acd341860237c1fcc16cafa>`_)
| - More flexible blocklist to capture all components with the same name, even if they have different species IDs (see `commit1 <https://github.com/Metaboverse/Metaboverse/commit/8975d24a556d31b2aa6e8013659bb80f22ff6a2a>`_ ; `commit2 <https://github.com/Metaboverse/Metaboverse/commit/1273b94acf1c1ee8fd4f60b175e61cf1bd506774>`_)
| - Find largest change possible for modifier regulation patterns (`see commit <https://github.com/Metaboverse/Metaboverse/commit/de1148b35d415cfa20ad3e68e47a3cbb3d729d25>`_)
| - Sort by best p-value (previously had taken a more conservative approach by using the worst p-value on each side of the reaction) (`see commit <https://github.com/Metaboverse/Metaboverse/commit/8975d24a556d31b2aa6e8013659bb80f22ff6a2a>`_)
| - Add button and capabilities to switch between inferred complex values or to compare each complex component individually within the reaction pattern (`commit1 <https://github.com/Metaboverse/Metaboverse/commit/31ece06c7476cc8d568bdd67f46dbceae2193d65>`_ ; `commit2 <https://github.com/Metaboverse/Metaboverse/commit/de1148b35d415cfa20ad3e68e47a3cbb3d729d25>`_)
|   - Will still display the complex as inferred value, but evaluated as each individual component during reaction pattern search
| - Protein complex inference/aggregation
|   - mean -> median for generating aggregate magnitude value for protein complex from component parts (`see commit <https://github.com/Metaboverse/metaboverse-cli/commit/e6755ca67322745dc40af89fdd67b894f5732fc8>`_)
|   - Aggregate statistic calculated using :data:`e * gmean(p-array)` (`see commit <https://github.com/Metaboverse/metaboverse-cli/commit/e6755ca67322745dc40af89fdd67b894f5732fc8>`_)
|   - Max aggregate p-value set to 1 (`see commit <https://github.com/Metaboverse/metaboverse-cli/commit/ce4ccad650f3e1bf51635e3415ca5759ab513f78>`_)
| - Allow exporting line plots for timecourse and multi-condition datasets (`see issue #89 <https://github.com/Metaboverse/Metaboverse/issues/89>`_)
| - Use user-provided names in labeling (`see issue #87 <https://github.com/Metaboverse/Metaboverse/issues/87>`_)
| - Toggle analyte labels on by default (`see commit <https://github.com/Metaboverse/Metaboverse/commit/1f79661240c196cdffd0114f91dcae51ed4e4ee1>`_)
| - Allow flexibility for input data where commas are used in place of decimals (`see issue #92 <https://github.com/Metaboverse/Metaboverse/issues/92#issuecomment-854090294>`_)
| - Remove duplicate rows from input data (interactive input will warn about these) (`see commit <https://github.com/Metaboverse/metaboverse-cli/commit/a2fc6642168adb3fc7bcc4e10e4b21aff4e272e3>`_)
| - Fix **Sustained** reaction pattern to not identify is input and output value being compared are identical (`see commit <https://github.com/Metaboverse/Metaboverse/commit/1273b94acf1c1ee8fd4f60b175e61cf1bd506774>`_)
| - Fix issue with :data:`parseComponents()` function where usage of modifiers in pattern determination was pushing all modifiers (catalysts and inhibitors) to reactants list (`see commit <https://github.com/Metaboverse/Metaboverse/commit/de1148b35d415cfa20ad3e68e47a3cbb3d729d25>`_)
| - For upregulated sustained reactions, get max of inputs and outputs (previously was getting min) (`see commit <https://github.com/Metaboverse/Metaboverse/commit/de1148b35d415cfa20ad3e68e47a3cbb3d729d25>`_)
| - Fixed global motif searching for pathway and perturbation visualization to search non-collapsed reaction dictionary too. (`see commit <https://github.com/Metaboverse/Metaboverse/commit/54a2e44d4913e1d4f903271bdae8af3617f0f33c>`_)
| - Added reaction pattern tooltips on button to show a graphical example of each (`see commit <https://github.com/Metaboverse/Metaboverse/commit/66d7ecc210c224451370772b4de3749af055aa69>`_)
| - Move some shared utilities to the `motif-utils.js` file (`see commit <https://github.com/Metaboverse/Metaboverse/commit/1273b94acf1c1ee8fd4f60b175e61cf1bd506774>`_)


=================================
v0.7.0
=================================
| **Important Note**
| Many of the changes introduced in :data:`v0.7.0` to session and intermediate file metadata will likely be incompatible with previous versions of Metaboverse.
|
| **Major**
| - :data:`Enzyme`/:data:`Metabolite` reaction patterns added: The :data:`Enzyme` reaction pattern evaluates for two neighboring reactions both with perturbed enzymes matching the given threshold. This will allow for better pattern identification, especially in RNA-seq/proteomics-only datasets. The :data:`Metabolite` reaction pattern looks for neighboring reactions both with perturbed metabolites matching the given threshold (see `issue #81 <https://github.com/Metaboverse/Metaboverse/issues/81>`_).
|
| **Minor**
| - Removed pathway-specific pattern detection: In our testing, these seemed to be minimally helpful.
| - Option added for users to define percentage of matching nodes between two reaction to allow for a collapse (see `issue #82 <https://github.com/Metaboverse/Metaboverse/issues/82>`_).
| - SVG export option (Full support for Inkscape, partial support for Illustrator) (see `issue #83 <https://github.com/Metaboverse/Metaboverse/issues/83>`_).
| - Improved and more explicit session data for all intermediate files (see `issue #78 <https://github.com/Metaboverse/Metaboverse/issues/78>`_).
| - Migrated source files to rutter.chpc.utah.edu/Metaboverse/source. This change should allow for faster downloads of pre-curated intermediate source files (see `issue #80 <https://github.com/Metaboverse/Metaboverse/issues/80>`_).
| - Loading icon in reaction pattern page to let user know patterns are loading, especially in cases where many reaction patterns are discovered and the software may take some time to display them all.



=================================
Previous versions
=================================

---------------------------------
v0.6.0
---------------------------------
| **Major**
| - New database integration: First supported release with the ability to overlay data on BiGG and BioModels network models and enable reaction pattern searching across a wider array of organisms. Note: Network models from these sources can be less uniform as Reactome sources, so users should exercise some caution when using these capabilities and perform some sanity checks (see `issue #73 <https://github.com/Metaboverse/Metaboverse/issues/73>`_).
| - kNN visualization improved to allow for more stable NN building without error (see `commit 2395cd6 <https://github.com/Metaboverse/Metaboverse/commit/2395cd6fe44167def52ae991b8db5f9559a9eba9>`_).
| - Neighbors dictionary backend curation is simplified and sped up (see `commit 355abd4 <https://github.com/Metaboverse/metaboverse-cli/commit/355abd4a6c5196bf6b4e46304eb1984d22597d7c>`_).
| - Improved security policies. Specifically, external websites are opened in an isolated browser window and explicitly are context isolated and unable to access node integration. Enabled GitHub and Reactome URLs (:data:`connect-src`) are more specific (see `commit 96b1c9f <https://github.com/Metaboverse/Metaboverse/commit/96b1c9fa3135cbe2aea97e4a132e57063acbcf38>`_).
|
| **Minor**
| - Progress bar during network build now update with more incremental steps for longer processes (see `issue #77 <https://github.com/Metaboverse/Metaboverse/issues/77>`_).
| - New variables for more unified backend processing. Metaboverse v0.6.0 and later will not be compatible with files curated using Metaboverse v0.5.0b or earlier.

---------------------------------
v0.5.0-beta
---------------------------------
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


---------------------------------
v0.4.0-beta
---------------------------------
| **Major**
| - Partial collapse: Addresses  `#51 <https://github.com/Metaboverse/Metaboverse/issues/51>`_  , which introduces partial collapsing to the reaction collapsing utility within Metaboverse. Previously, perfect matches were required between two reactions to collapse the reactions. However, this can be overly stringent in key metabolic pathways where a metabolite that is output by one reaction may not be required for the subsequent reaction (perhaps ATP is produced by reaction A but is not required for reaction B). To perform a partial collapse, Metaboverse operates by largely the same scheme as before, but now checks for a perfect match from each neighboring reaction, and if a perfect match is not available, checks for partial matches by filtering out high-degree nodes (quartile 98 of all non-reaction node degrees) and then checking if at least 30% of the nodes match with its neighbor.
| - Improvements to nearest neighbor searches where all iterations of a species are included in the graphing. Previously, it would only use the literally selected node to search for neighbors, but Reactome provides separate species IDs for a metabolite's different organelle-localizations, which was complicating these searches.

| **Minor**
| - Displays a preview of the user-selected reaction in an interactive format on the Pattern Search Analysis page until the user selects a Pathway to visualize. If a reaction is collapsed and spans two pathways, no pathways will be shown and instead a note that the reaction spans two pathways is displayed.
| - The Pattern Search Analysis page now allows users to filter out collapsed reactions from the search results. By default, collapsed reactions will be displayed until the checkbox is unchecked by the user.
| - Minor updates to Pattern Search Analysis page to make better usage of blank space
| - Fixes  `#60 <https://github.com/Metaboverse/Metaboverse/issues/60>`_  , where the :data:`.mvrs` file extension would not be automatically added to the user-provided output file name in Linux.
| - Addresses  `#62 <https://github.com/Metaboverse/Metaboverse/issues/62>`_  , where the some time-course/multi-condition slider bars would be improperly formatted.


---------------------------------
v0.3.3-beta
---------------------------------
| **Minor**
| - Closes `#63 <https://github.com/Metaboverse/Metaboverse/issues/63>`_ by applying :data:`safestr()`` function to all user input encodings to make sure no errors arise.


---------------------------------
v0.3.2-beta
---------------------------------
| **Minor**
| - Closes  `#59 <https://github.com/Metaboverse/Metaboverse/issues/59>`_  where non-ascii characters in reaction names would break the info extraction. Added a safestring conversion utility to prevent ascii-character issues.


---------------------------------
v0.3.1-beta
---------------------------------
| **Minor**
| - Fixes path separator for motif page name identification to allow for including modifiers in motif ID and exclusion of hubs ( `#55 <https://github.com/Metaboverse/Metaboverse/issues/55>`_ )
| - Fixes CHEBI mapping so that CHEBI IDs provided as input data are more reliably used as mapping IDs if it cannot match the metabolite by name ( `#58 <https://github.com/Metaboverse/Metaboverse/issues/58>`_ )
| - Fixes issue that arose in :data:`v0.3.0b` where some motif stamps could not be clicked on for viewing for timecourse/multi-condition data where it could not identify the shape for an unknown component type ( `#54 <https://github.com/Metaboverse/Metaboverse/issues/54>`_ )
| - Addresses  `#59 <https://github.com/Metaboverse/Metaboverse/issues/59>`_  where non-ascii characters in reaction names would break the info extraction. Was not able to recapitulate the error, but this fix, where relevant reaction metadata is forced to a string data-type, should allow for some flexibility here.
| - Updates walkthroughs and documentation to address ( `#31 <https://github.com/Metaboverse/Metaboverse/issues/31>`_ ) and update formatting


---------------------------------
v0.3.0-beta
---------------------------------
| **Major**
| - Allows for more flexible gene/protein mapping with Reactome-formatted node names. For example, Reactome will label a gene or protein with its isomer coordinates. Metaboverse now ignores those coordinates during attribute mapping of the user's data.
| - New naming of modules: :data:`Motif Search` is now called :data:`Pattern Analyis`, :data:`Visualize` is now called :data:`Explore`, and :data:`Connectivity` is now called :data:`Perturbation Networks`. Changed to be more descriptive and accessible to all users from broader backgrounds
|
| **Minor**
| - Fixed nearest neighbors capabilities in Perturbation Network visualization. In a previous release, a change had interfered with its function.
| - Use of outdated version will now direct user to the download page for the most recent version
| - The :data:`Back` button from any of the analysis modules will now redirect back to the index page
| - Removed compartment visualization from the :data:`Perturbation Networks` page as these often would just clutter the visualization and would not actually be helpful
| - Fixed reaction node formatting to turn off motif symbols/highlighting when timepoint/condition changed as this had been disabled by a previous change
| - Fixed collapsed reactions to ensure they were included in all reaction and motif formatting events
| - General formatting changes
| - Updated documentation


---------------------------------
v0.2.0-beta
---------------------------------
| **Major**
| - Fixes issues with missing metabolites during the network mapping stage ( `#37 <https://github.com/Metaboverse/Metaboverse/issues/37>`_ ). Addressed by re-working the metabolite synonym scheme to provide the same coverage of synonym look-up as before, but with more robustness so that some that were difficult to map would start mapping (i.e., Fructose)
| - Added dynamic line-plots of a selected reaction motif for time-course and multi-condition data ( `#15 <https://github.com/Metaboverse/Metaboverse/issues/15>`_). When exploring motifs on the Motif page for time-course and multi-condition experiments, a new panel appears at the bottom of the page which, for a selected motif, will show those reaction motif's component's behavior across all the time-points or conditions.
| - Added the option to exclude reaction motifs at a given time-point or condition that appear in another selected time-point or condition. ( `#16 <https://github.com/Metaboverse/Metaboverse/issues/16>`_ ). For example, if a user has selected to view motifs at a terminal time-point, but they want to know which reactions are motifs at this time-point but not at the initial time-point, they can exclude the motifs that show up at both time-points using the appropriate drop-down menu on the motif page.
|
| **Minor**
| - Metaboverse now outputs a table of unmapped metabolites ( `#35 <https://github.com/Metaboverse/Metaboverse/issues/35>`_ ).
| - Exploration pages now have pop-out bubbles with all information for compartments and node/link types `7d17d34 <https://github.com/Metaboverse/Metaboverse/commit/7d17d34aca5e900c307e266a07b4d82bd19a222d>`_.
| - Metaboverse new remembers and provides session info for experiment name, experiment type, labels, etc. and automatically fills those out for the user if returning to a page within the session `172d21a <https://github.com/Metaboverse/Metaboverse/commit/172d21a719bbc855fd46d4d8da223140c512a18f>`_.
| - Updated minor page formatting to make display more stable between Windows/Linux/Mac `52a100d <https://github.com/Metaboverse/Metaboverse/commit/52a100da0958af75c489165bc2f7c9eaf80294e8>`_.
| - Added test cases to CI for new/updated features
| - Updated package dependency information
| - Updated docs and FAQs


---------------------------------
v0.1.4-beta
---------------------------------
| - Fixes `#26 <https://github.com/Metaboverse/Metaboverse/issues/26>`_, where an error log is output if build fails
| - Removes direct Matplotlib imports in metaboverse-cli modules to prevent unnecessary errors and incompatibilities


---------------------------------
v0.1.3-beta
---------------------------------
| - Fixes bug where user paths with spaces were unable to be used ( `#26 <https://github.com/Metaboverse/Metaboverse/issues/26>`_ )


---------------------------------
v0.1.2-beta
---------------------------------
| - Fixes bug that prevented the curation from running without a blocklist ( `#19 <https://github.com/Metaboverse/Metaboverse/issues/19>`_ )
| - Fixes bug during data mapping that caused protein or gene values to occasionally map to metabolites ( `#20 <https://github.com/Metaboverse/Metaboverse/issues/20>`_ )


---------------------------------
v0.1.1-beta
---------------------------------
| - Fixes minor run-time issues with the Metaboverse interactive app
| - Fixes version alert to let users know if there is a newer version of Metaboverse available


---------------------------------
v0.1.0-beta
---------------------------------
| Initial beta release


---------------------------------
v0.0.1-beta
---------------------------------
| Demo pre-release with included human network data file for network visualization and exploration. Currently only available for MacOS.
|
| How to run:
|
| 1. Download attached :data:`.zip` demo file.
| 2. Double-click on :data:`.zip` file to uncompress Metaboverse and the accompanying test file
| 3. Within the uncompressed folder, right-click on Metaboverse to launch the app
| 4. Drag and drop the :data:`.json` file to the appropriate load icon and click the Visualize button.

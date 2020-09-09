###############
Updates
###############

===========
v0.3.0-beta
===========
| **Major**
| - Allows for more flexible gene/protein mapping with Reactome-formatted node names. For example, Reactome will label a gene or protein with its isomer coordinates. Metaboverse now ignores those coordinates during attribute mapping of the user's data.

| **Minor**
| - 

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

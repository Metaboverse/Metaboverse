##############
XPRESSpipe
##############
|build-status| |docs| |Docker|

=================
About
=================
| `XPRESSpipe <https://github.com/XPRESSyourself/XPRESSpipe>`_ is a part of the `XPRESSyourself <https://github.com/XPRESSyourself/>`_ suite of sequencing tools. `XPRESSpipe <https://github.com/XPRESSyourself/XPRESSpipe>`_ is an automated, efficient, and flexible pipeline for end-to-end processing of ribosome profiling data.
|
| `XPRESSpipe <https://github.com/XPRESSyourself/XPRESSpipe>`_ is currently capable of handling single-end (SE), paired-end (PE), and ribosome profiling data.
|
| If you have limited or no computational experience, please see our :ref:`beginners_link`.
|
| Please refer to the :ref:`overview_link` page for more details regarding functionality.
|
| For each sequencing type pipeline, all relevant quality control analyses are performed for that sequencing type, as well as gene coverage analysis of a housekeeping gene.
|
| Other analyses can be performed by `XPRESSplot <https://github.com/XPRESSyourself/XPRESSplot>`_. Please read the relevant documentation for more information.

=================
Table of contents
=================
.. toctree::
   :hidden:

   self

.. toctree::
   :maxdepth: 1

   content/overview
   content/quickstart
   content/beginner
   content/installation
   content/general-usage
   content/reference-building
   content/seRNAseq
   content/peRNAseq
   content/riboseq
   content/quality-control
   content/analysis
   content/trim
   content/align
   content/count
   content/normalization
   content/other-features
   content/faqs
   content/updates

=======
License
=======
| `XPRESSpipe <https://github.com/XPRESSyourself/XPRESSpipe>`_ and the `XPRESSyourself <https://github.com/XPRESSyourself/>`_ suite is developed and maintained by Jordan Berg in the `Rutter Lab <https://biochem.utah.edu/rutter/index.html>`_ @ the `University of Utah <https://www.utah.edu/>`_, along with other collaborators. We welcome pull requests if you would like to contribute to the project.
|
| `XPRESSpipe <https://github.com/XPRESSyourself/XPRESSpipe>`_ is perpetually open source under a GNU General Public License (v3.0).

==========
Questions?
==========
| If you have questions, requests, or bugs to report, please use the `XPRESSpipe issues forum <https://github.com/XPRESSyourself/XPRESSpipe/issues>`_.



.. |build-status| image:: https://travis-ci.org/XPRESSyourself/XPRESSpipe.svg?branch=master
    :target: https://travis-ci.org/XPRESSyourself/XPRESSpipe
    :alt: Build Status

.. |codecov| image:: https://codecov.io/gh/XPRESSyourself/XPRESSpipe/XPRESSpipe.svg?branch=master
    :target: https://codecov.io/gh/XPRESSyourself/XPRESSpipe
    :alt: Code Coverage

.. |docs| image:: https://readthedocs.org/projects/xpresspipe/badge/?version=latest
    :target: https://xpresspipe.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. |Docker| image:: https://img.shields.io/static/v1.svg?label=docker&message=dowload&color=informational
    :target: https://cloud.docker.com/repository/docker/jordanberg/xpresspipe/general
    :alt: Docker

####################
Backend Development
####################

.. note::
    This page includes information about the underlying scripts of Metaboverse and Metaboverse-cli and how to contribute to and modify the software.

===============================
Expectations
===============================
| It is essential before getting started that you read through the `Contributing <https://github.com/Metaboverse/Metaboverse/blob/main/CONTRIBUTING.md>`_ and `Code of Conduct <https://github.com/Metaboverse/Metaboverse/blob/main/CODE_OF_CONDUCT.md>`_ documents. Harassment, bullying, or other inappropriate behavior of any form will absolutely not be tolerated. Those violating these standards will be banned from contributing to the Metaboverse project.

===============================
Communication
===============================
| We highly recommend making use of the `Issues <https://github.com/Metaboverse/Metaboverse/issues>`_, `Discussions <https://github.com/Metaboverse/Metaboverse/discussions>`_, and `Projects <https://github.com/Metaboverse/Metaboverse/projects>`_ pages to communicate developments, issues, questions and more.

===============================
First Steps
===============================
| The code base for Metaboverse is currently housed in two primary Github repositories. `metaboverse-cli <https://github.com/Metaboverse/metaboverse-cli>`_ handles the backend processing of metabolic networks and preparing user data for analysis, and `Metaboverse <https://github.com/Metaboverse/Metaboverse>`_ handles front-end processing and visualization of user data.

| If you would like to modify either of these repositories, you should be begin by forking a branch of the repository. You should reference the documentation: `https://docs.github.com/en/get-started/quickstart/fork-a-repo <https://docs.github.com/en/get-started/quickstart/fork-a-repo>`_.

Once you have made the desired changes, you can submit a Pull Request, which will be reviewed and incorporated into the code base for future deployment if the changes meet the required standards. You should reference the documentation: `https://docs.github.com/en/pull-requests <https://docs.github.com/en/pull-requests>`_.

===============================
:data:`metaboverse-cli`
===============================
| :data:`metaboverse-cli` handles a large portion of the backend processing, primarily to curate the metabolic network for the organism of interest and to layer the user's data onto said network. The source code is generally organized as follows:

.. code-block:: shell

    metaboverse-cli
    ├── setup.py:                The command line tool setup file 
    ├── requirements.txt:        The depencency file used for distribution 
    ├── metaboverse-cli.spec:    The Pyinstaller file to compile the distributable executable
    └── /metaboverse_cli          
        ├── arguments.py:        A file to set and parse user-provided command line arguments
        ├── __main__.py:         The main execution script for the backend
        ├── /mapper:             This submodule generates the metabolite synonym mapper that assists 
        │                          in mapping user data to the metabolic network
        ├── /curate:             This submodule builds the organism specific reaction network and 
        │                          stores any relevant metadata about the reaction network component
        └── /analyze:            This submodule handles layering the user's data onto the metabolic 
                                   network. This includes implementing the metabolite synonym mappper.
                                   This submodule also handles generating the reaction network with 
                                   reaction collapsing


===============================
:data:`Metaboverse`
===============================
| :data:`Metaboverse` handles the front-end user interface and data visualization, and performed reaction pattern recognition in real-time. The source code is generally organized as follows:

.. code-block:: shell

    Metaboverse      
    ├── build.sh:                A helpful script for compiling Metaboverse for distribution 
    ├── /docs                    Source files for metaboverse.readthedocs.io 
    └── /app          
        ├── /css:                CSS style files, including both custom and distributed files
        ├── /data:               Test data that is distributed with each release of Metaboverse, 
        │                          along with icons and images used within the GUI
        ├── /html:               HTML files for each of the pages used with in the Metaboverse GUI
        ├── /js:                 Javascript files for the Metaboverse GUI. Selected important files
        │   │                      are listed below:
        │   ├── index.js:        Handles the home screen interactions
        │   ├── curate.js:       Handles the Curate screen interactions, where users provide 
        │   │                      information about their model organism and more 
        │   ├── variables.js:    Handles the Variables and Data screen interactions, where users 
        │   │                      provide input data and other experimental data 
        │   ├── build.js:        Handles the Build screen interactions, where Metaboverse compiles 
        │   │                      the information provided by the user and performs back-end 
        │   │                      processing with metaboverse-cli 
        │   ├── visualize.js:    Handles network visualization 
        │   ├── timecourse.js:   Handles formatting and toggling between time-course or multi-condition
        │   │                      samples during visualization
        │   ├── js-colormaps.js: A Javascript implementation of the "seismic" colormap from Matplotlib
        │   ├── motif-script.js: Front handler for the Reaction Pattern screen 
        │   ├── motif-graph.js:  Handles interactions and visualization of Reaction Patterns 
        │   ├── motifs.js:       Handles the actual, real-time Reaction Pattern searching across all 
        │   │                      reactions in the network
        │   ├── motif-global.js: Searches for every possible Reaction Pattern during Pathway visualization
        │   ├── perturbations.js:Handles the Perturbation screen interactions and network generation of 
        │   │                      a network consisting of all perturbed reactions within a given pathway
        │   └── datatable.js     Handles the Format Data screen interactions          
        ├── main.js:             This handles initializing the Electron app window
        ├── package.json         This handles initializing dependencies and other settings for an 
        │                          an Electron app
        └── __version__.txt:     This supplies the current software version displayed within the 
                                   software

===============================
Documentation and Testing
===============================
| Any changes made should be supported by intuitive and descriptive comments throughout the code, documentation for analysis features for `metaboverse.readthedocs.io <metaboverse.readthedocs.io>`_, and with test cases.

===============================
Distribution
===============================
| Distribution of a new release of Metaboverse requires several steps, outlined below. 
| 1) Compile :data:`metaboverse-cli`
|       a) Update the version number in :data:`metaboverse_cli/__init__.py`
|       b) You will then need to do the following for each operating system you are distibuting Metaboverse on.

.. code-block:: shell

    $ conda create --name pyinstaller
    $ conda activate pyinstaller
    $ conda config --add channels conda-forge
    $ conda install python=3.8

| On Windows:

.. code-block:: shell

    $ conda install pyinstaller 
    $ conda install --file requirements.txt

| On Mac/Linux:

.. code-block:: shell

    $ pip install pyinstaller 
    $ pip install -r requirements.txt  

| Then, on each operating system

.. code-block:: shell

    $ cd /path-to/metaboverse-cli 
    $ pyinstaller metaboverse-cli.spec

| From here, you can archive each operating system executable on the Github release page for :data:`metaboverse-cli`.

| 2) Compile :data:`Metaboverse`
|       a) You will need to download the :data:`metaboverse-cli` executables for each operating system. These should be stored in the directory, :data:`/path-to/Metaboverse/app/python/`
|       b) Update version numbers in :data:`/path-to/Metaboverse/docs/conf.py`, :data:`/path-to/Metaboverse/build.sh`, :data:`/path-to/Metaboverse/app/__version__.txt`, :data:`/path-to/Metaboverse/app/package.json`, and :data:`/path-to/Metaboverse/CITATION.cff`
|       c) Compile the Metaboverse distributions:

.. code-block:: shell

    $ cd /path-to/Metaboverse
    $ bash build.sh


.. warning:: 
    You will likely need to compile the Metaboverse package for Mac on a Mac. If the Mac-version :data:`metaboverse-cli` executable is packaged into Metaboverse via Electron on another operating system, we have consistently received an error in executing :data:`metaboverse-cli`.

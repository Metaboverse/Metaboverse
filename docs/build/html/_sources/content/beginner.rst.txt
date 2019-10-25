.. _beginners_link:

################
Beginner's Guide
################

=================================
First Steps
=================================
| If this is your first time doing any programming, congratulations! You are embarking upon a very rewarding path. As with learning any new spoken language, there is a learning curve associated with learning a computer language. While XPRESSpipe is aimed at reducing the majority, if not (hopefully) all of the overhead associated with processing this data, using this software will still require some effort, just as would learning any new language or laboratory technique.

| XPRESSpipe is run through the `command line interface <https://en.wikipedia.org/wiki/Command-line_interface>`_ (or `CLI <https://www.youtube.com/watch?v=kqUR3KtWbTk>`_). This may seem daunting, but luckily, several free online courses are available to quickly catch you up to speed on some of the basics that will be required to use this software. We recommend Codecademy's CLI course, which you can find `here <https://www.codecademy.com/learn/learn-the-command-line>`_ and should take only a couple of hours (Codecademy estimates ~10 hours, but you probably don't need to finish the course to be ready to use XPRESSpipe. The purpose of this is to help you become more comfortable with the command line). We recommend watching the walkthrough videos found on the :doc:`quickstart <./quickstart>` page.

| Once you're ready to jump into the command line, we can get rolling! For the steps below, we're going to assume we are on an Mac operating system and provide examples under this pretext, but this software is compatible with any Linux-like operating system and the syntax is largely the same (sorry Windows users! -- but if you have a newer version of Windows, you may be able to use a Linux-flavored environment).

=================================
Install XPRESSpipe
=================================
| - Please refer to the :doc:`installation documentation <./installation>` or the walkthrough video below:

.. raw:: html

    <embed>
        <script id="asciicast-256347" src="https://asciinema.org/a/256347.js" async></script>
    </embed>

=================================
Generate Reference Files
=================================
| - Before we can process our raw RNA-seq data, we need to create a reference directory (or for a folder, in other terms). In this example, we will be working with human-derived RNA-seq data, so let's perform the following in the command line:

.. code-block:: shell

  $ cd ~/Desktop
  $ mkdir reference_folder
  $ mkdir reference_folder/fasta_files

| 1. The first command helped us navigate to the Desktop. The "~" icon is a shortcut for the User directory, and every directory needs to be separated by a "/"
| 2. The second command created a new folder in the Desktop directory called :data:`reference_folder`
| 3. The third command created a new folder in the reference directory for intermediate reference files

| - Now let's get the reference files. We're going to do this directly in the command line, but if you have trouble with this, I will explain an alternative afterwards. Quick note, because the next lines of code are a bit long, I used the "\" character to indicate I am continuing the command in the next line. You not include these characters when executing the command, they just help make the code more readable. We will first read the retrieval commands into a file which will additionally act as a log file for the version for the genome version we are using.
| - You should modify the the variable calls between the :data:`#` signs. For :data:`GTF_URL`, you should change the URL currently provided to the one appropriate for your organism of interest. Make sure you are downloading the GTF file and NOT the GFF file. For :data:`FASTA_URL`, you should do the same as before with the URL to the chromosome DNA FASTA files, but you should only copy the URL up to "chromosome", but not include the chromosome identifier. For :data:`CHROMOSOMES`, type out the chromosome identifiers you want to download between the " characters with a space between each.

.. code-block:: shell

  $ cd reference_folder/

  ### Change these ###
  $ echo 'GTF_URL=ftp://ftp.ensembl.org/pub/release-97/gtf/homo_sapiens/Homo_sapiens.GRCh38.97.gtf.gz' >> fetch.sh
  $ echo 'FASTA_URL=ftp://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome' >> fetch.sh
  $ echo 'CHROMOSOMES="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y"'
  ####################

  $ echo 'curl -O $GTF_URL' >> fetch.sh
  $ echo 'gzip -d Homo_sapiens.GRCh38.97.gtf.gz' >> fetch.sh
  $ echo 'mv Homo_sapiens.GRCh38.97.gtf transcripts.gtf' >> fetch.sh
  $ echo 'cd fasta_files/' >> fetch.sh
  $ echo 'for X in $CHROMOSOMES; ' >> fetch.sh
  $ echo 'do curl -O ftp://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.${X}.fa.gz; done ' >> fetch.sh
  $ echo 'gzip -d *.gz' >> fetch.sh
  $ echo 'cd ../' >> fetch.sh
  $ bash fetch.sh

| - Let's discuss what we just did:
| 1. We navigated into the reference folder, downloaded a GTF reference file and unzipped it, then navigated to the :data:`fasta_file` directory to download the raw reference data and unzipped it. Finally, we returned to the main reference directory.
| 2. If this didn't work, we can navigate to `Ensembl <https://www.ensembl.org/>`_ to download the relevant data. We need to get the GTF and DNA chromosomal FASTA files for our organism of interest. The link to the chromosome sequence files actually contains more files than we need. We just need the files that start with :data:`Homo_sapiens.GRCh38.dna.chromosome`. You can download them, move them to the appropriate directories within your reference directory, and unzip the files by double-clicking on them.

| - Now we need to curate these references files into something the sequencing alignment software can use. Since we are using ribosome profiling data, we want a reference that will allow us to `avoid mapping to the 5' and 3' ends of genes <https://www.cell.com/cms/10.1016/j.celrep.2016.01.043/attachment/257faf34-ff8f-4071-a642-bfdb531c75b8/mmc1>`_. We also don't want to align to anything but protein coding genes. Finally, we want to quantify to the longest transcript. This last bit just helps the software avoid confusion when a gene has multiple splice variants to choose from. Since this is short read sequencing (let's say we were doing 50 bp single-end sequencing), we also want to factor this into the curation of the reference (see the :data:`--sjdbOverhang` argument below).

.. code-block:: shell

  $ xpresspipe curateReference \
                --output ./ \
                --fasta fasta_files/ \
                --gtf ./transcripts.gtf \
                --protein_coding \
                --truncate \
                --sjdbOverhang 49

  ### or ###

  $ xpresspipe build

  ### And then choose the curate option ###


| - The truncation option is only necessary when using XPRESSpipe to process ribosome profiling samples and their associated RNA-seq samples.
| - If interested in quantifying miRNA, etc, leave out the :data:`--protein_coding` argument.
| - If running sequencing where the read (single-end) or mates not equal to 100 bp, you will want to change the :data:`--sjdbOverhang` argument to be the length of one of the paired-end reads - 1, so if we ran 2x100bp sequencing, we would specify :data:`--sjdbOverhang 99` (although in this case, the default of :data:`--sjdbOverhang 100` is just fine). If you changed this number, remember this for the next steps as you will need to provide it again if changed here.
| - This may take awhile, and as we will discuss later, you may want to run these steps on a supercomputer, but this will serve as a preliminary guide for now.
| - One final consideration -- if we are dealing with an organism with a smaller genome size, we will want to provide the :data:`--genome_size` parameter with the the number of nucleotides in the organism's genome. If you change this parameter in this step, you will need to provide the parameter and value in the :data:`align`, :data:`riboseq`, :data:`seRNAseq`, and :data:`seRNAseq` modules.

=================================
Process Raw Sequencing Files
=================================
| - Now let's get our raw data::
| 1. Make a new folder, something called :data:`raw_data` or whatever you like and place your data there.
| 2. Make sure the files follow proper naming conventions (see naming conventions at :ref:`general_link`)
| 3. Now let's process the data
| 4. Let's also create a folder called something like :data:`output`
| 5. Also, make sure you have the 3' adapter sequence handy used when generating your sequencing library
| 6. We'll feed the program the new GTF file that contains only longest transcript, protein coding, truncated references generating in the reference curation step
| 7. We'll give the experiment a name and also specify what `method of sample normalization <https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/>`_ we want performed on the count data
| 8. We also need to specify the :data:`--sjdbOverhang` amount we fed into the reference curation step, so in this case we will use :data:`--sjdbOverhang 49`

.. code-block:: shell

  $ xpresspipe riboseq --input raw_data/ \
                      --output output/ \
                      --reference reference_folder/ \
                      --gtf reference_folder/transcripts_LCT.gtf
                      --experiment riboseq_test
                      --adapter CTGTAGGCACCATCAAT
                      --method RPKM
                      --sjdbOverhang 49

  ### or ###

  $ xpresspipe build

  ### And then choose the appropriate pipeline to build


| - If you are running a lot of files, especially for human samples, this may take a lot of time. We recommend running this on some kind of server. A situation like yeast with few samples may be feasible to run on a personal computer, but will likely also take some time.

------------------
Sequencing Metrics
------------------
| In your output folder, you will see a file named :data:`riboseq_test_multiqc_report.html`. This file will compile the statistics from each processing step of the pipeline for each sample file you provided as input. Things like read quality, mapping, and quantification statistics can be found here. Just double-click the file or execute the following command to open in your default browser window.

.. code-block:: shell

  $ open riboseq_test_multiqc_report.html

------------------
Library Complexity
------------------
| Within the :data:`complexity` directory in your output folder, you will find summary PDFs for all samples processed analyzing library complexity of each sample.

-------------------
Metagene Analysis
-------------------
| Within the :data:`metagene` directory in your output folder, you will find summary PDFs for all samples processed analyzing the metagene profile of each sample.

--------------------------------
Periodicity (Ribosome Profiling)
--------------------------------
| Within the :data:`periodicity` directory in your output folder, you will find summary PDFs for all samples processed analyzing ribosome periodicity of each of each sample containing reads 28-30nt.

----------------------------------
Count Data and Downstream Analysis
----------------------------------
| Within the :data:`counts` directory in your output folder, you will find individual counts tables for each sample, as well as compiled tables for each sample that was processed.


=======================
Supercomputing
=======================
--------------------
Install
--------------------
| - Much of the same commands will be performed as above, aside from a couple exceptions:
| 1. When installing XPRESSpipe, you need to provide a location for personal storage of the software:

.. code-block:: shell

  $ python setup.py install --prefix ~/.local/bin


| 2. Add this path to your :data:`$PATH`:

.. code-block:: shell

  $ echo 'export PATH="~/.local/bin:$PATH"' >> ~/.bash_profile
  $ echo 'export PYTHONPATH="/uufs/chpc.utah.edu/common/home/u0690617/.local/bin/lib/python3.7/site-packages"' >> ~/.bash_profile

| 3. Let's test this to make sure everything is operating properly:

.. code-block:: shell

  $ cd ~/
  $ xpresspipe test

---------------
Run Data
---------------
| - The commands here are the same as above, but likely the method of execution will be different. A lot of supercomputing clusters manage job submission through a system called `SLURM <https://www.youtube.com/watch?v=RpkAyFI05yY>`_. Likely, the supercomputing cluster you are running your data on will have instructions for how to use this, but briefly, here is an example batch script (should end in the suffix :data:`.sh`):

.. code-block:: shell

  #!/bin/bash
  #SBATCH --time=72:00:00
  #SBATCH --nodes=1
  #SBATCH -o /scratch/general/lustre/$USER/slurmjob-%j
  #SBATCH --partition=this_cluster_has_no_name

  #set up the temporary directory
  SCRDIR=/scratch/general/lustre/$USER/$SLURM_JOBID
  mkdir -p $SCRDIR

  # Provide location of raw data and parent reference directory
  SRA=/scratch/general/lustre/$USER/files/your_favorite_experiment_goes_here
  REF=/scratch/general/lustre/$USER/references/fantastic_creature_reference

  # Send raw data to your Scratch directory
  mkdir $SCRDIR/input
  cp $SRA/*.fastq $SCRDIR/input/.

  # Make an output directory
  mkdir $SCRDIR/output
  cd $SCRDIR/.

  xpresspipe riboseq -i $SCRDIR/input -o $SCRDIR/output/ -r $REF --gtf $REF/transcripts_CT.gtf -e this_is_a_test -a CTGTAGGCACCATCAAT --sjdbOverhang


| - To queue this script into the job pool, you would do the following:

.. code-block:: shell

  $ sbatch my_batch_script.sh

| - To monitor the progress of your job, execute the following:

.. code-block:: shell

  $ watch -n1 squeue -u $USER

| - After the job is finished, you can export the data as shown in the next section.



----------------
Explore the Data
----------------
| Once the data is finished processing, we can start exploring the output. Explanations each quality control analysis can be found in the :ref:`analysis_link` section of the documentation.
| In order to get the data from a HPC to your personal computer, you can use a command like the following:

.. code-block:: shell

  $ cd ~/Desktop # Or any other location where you want to store and analyze the data
  $ scp USERNAME@this_cluster_has_no_name.chpc.university.edu:/full/path/to/files/file_name.suffix ./

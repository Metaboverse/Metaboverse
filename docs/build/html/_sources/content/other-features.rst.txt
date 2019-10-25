############################
Other Features
############################

=================================
Convert Counts Table Gene Names
=================================
| Count tables are often produced with systematic names used to label each gene. The following sub-module will convert the column of systematic gene names to a common name using a reference GTF file

-----------
Arguments
-----------
| The help menu can be accessed by calling the following from the command line:

.. code-block:: shell

  $ xpresspipe convertNames --help

.. list-table::
   :widths: 35 50
   :header-rows: 1

   * - Required Arguments
     - Description
   * - :data:`-i \<path/filename\>, --input \<path/filename\>`
     - Path and file name to sequence dataframe
   * - :data:`-g \</path/transcripts.gtf\>`, :data:`--gtf \</path/transcripts.gtf\>`
     - Path and file name to GTF

.. list-table::
    :widths: 35 50
    :header-rows: 1

    * - Optional Arguments
      - Description
    * - :data:`--orig_name_label \<label\>`
      - Label of original name (usually "gene_id ")
    * - :data:`--orig_name_location \<position\>`
      - Position in last column of GTF where relevant data is found (i.e. 0 would be the first sub-string before the first comma, 3 would be the third sub-string after the second comma before the third comma)
    * - :data:`--new_name_label \<label\>`
      - Label of original name (usually "gene_id ")
    * - :data:`--new_name_location \<position\>`
      - Position in last column of GTF where relevant data is found (i.e. 0 would be the first sub-string before the first comma, 3 would be the third sub-string after the second comma before the third comma)
    * - :data:`--refill \<label\>`
      - In some cases, where common gene names are unavailable, the dataframe will fill the gene name with the improper field of the GTF. In this case, specify this improper string and these values will be replaced with the original name

-----------
Examples
-----------
| **Example 1 -- Convert gene names in count dataframe**

.. ident with TABs
.. code-block:: python

  $ xpresspipe convertNames -i riboprof_out/counts/se_test_counts_table.csv -g se_reference/transcripts.gtf

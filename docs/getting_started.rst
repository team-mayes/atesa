Getting Started with ATESA
==========================

Installation
------------

[Insert finalized installation instructions here]

Usage
-----

ATESA is designed to dynamically handle jobs on a PBS/Torque or Slurm batch system. It is invoked as:

.. code-block:: bash

   atesa <config_file> [<working_directory>]
   
The config_file parameter is required and supports absolute or relative paths. The working_directory parameter is optional, but if it is provided it overrides the value in the configuration file. This option is made available for use on systems that only allocate working space for jobs after they have been initialized.

Although ATESA can be run directly from the command line, because its process continues for the entirety of the job it is usually recommended to submit it as its own batch job. A single core and a modest allocation of memory should be sufficient to run ATESA on most platforms (it is neither memory- nor processor-intensive, although certain features involve a significant amount of I/O).

The Configuration File
----------------------

ATESA comes packaged with a configuration (or "config") file containing the default value for every setting; however, these defaults are set internally to the software and will still be used if the line that sets them is omitted from the configuration file.

Certain settings should be given regardless of the job type, while others are necessary for only a subset of job types. Recommendations for each job type as well as detailed documentation of each setting are available on the :ref:`TheConfigFile` page.
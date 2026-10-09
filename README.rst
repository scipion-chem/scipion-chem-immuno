================================
Immuno scipion plugin
================================

**Documentation under development, sorry for the inconvenience**

Scipion framework plugin for the use of immunoinformatics tools comming from several sources:

- IIITD: for epitope prediction and evaluation
- Vaxign-ML: for protective antigen prediction
- LANL/CATNAP: for cross-referencing epitope candidates against known broadly neutralizing antibodies

===================
Install this plugin
===================

You will need to use `3.0.0 <https://github.com/I2PC/scipion/releases/tag/v3.0>`_ version of Scipion
to run these protocols. To install the plugin, you have two options:

- **Stable version**  

.. code-block:: 

      scipion installp -p scipion-chem-immuno
      
OR

  - through the plugin manager GUI by launching Scipion and following **Configuration** >> **Plugins**
      
- **Developer's version** 

1. **Download repository**:

.. code-block::

            git clone https://github.com/scipion-chem/scipion-chem-immuno.git

2. **Switch to the desired branch** (main or devel):

Scipion-chem-immuno is constantly under development.
If you want a relatively older an more stable version, use main branch (default).
If you want the latest changes and developments, user devel branch.

.. code-block::

            cd scipion-chem-immuno
            git checkout devel

3. **Prerequisites**

- IIITD:
This package uses web browser to access the software servers online.
Therefore, you need to specify which browser to use and its location.
Do so editing the scipion.conf file and add the variables:
    - IIITD_BROWSER = firefox/chrome/chromium  (defines the browser to use, that must already be installed in your computer)
    - IIITD_BROWSER_PATH = <path/to/browser>   (defines the location of the binary for the browser use)

- Vaxign-ML:
This package runs using a docker image. Since docker images need special permission to be downloaded, the user needs to
be included in the "docker" bash group of the machine to be able to use it. Managing these groups need sudo permissions.

To create the "docker" group:
.. code-block::
            sudo groupadd docker

Then add your user to the group:
.. code-block::
            sudo usermod -aG docker $USER

Once your user is in the docker group, the installation can proceed normally.

- LANL/CATNAP reference databases:
These are not installed automatically and are not shipped with the plugin. Their
terms of use do not clearly permit redistribution, and hiv.lanl.gov exposes no
stable download URL for the antibody database, only an interactive search form,
while explicitly discouraging automated traffic. They are therefore downloaded
manually, once.

Required, the antibody database. Open the LANL HIV Molecular Immunology Database
at https://www.hiv.lanl.gov/content/immunology/, go to the antibody search
section, run a search with no filters so that every record is returned, and
export the result as CSV. Save it as 'ab_all.csv'.

Optional, the neutralization panel. Open CATNAP at
https://www.hiv.lanl.gov/components/sequence/HIV/neutralization/ and download the
antibody summary table. It adds the mean panel IC50 and the number of viruses
tested to each match. The protocol runs without it, reporting the potency columns
as empty.

Save both files under the scipion/software/em folder, in a dedicated
lanl-catnap subfolder, which is where this plugin's packages are installed.
Any other readable location also works, as the protocol resolves the files
through the variables below rather than by searching for them. Then add these
variables to the scipion.conf file:
    - LANL_AB_ALL_PATH = <path/to/ab_all.csv>   (required, LANL antibody database export)
    - CATNAP_ABS_PATH = <path/to/catnap_abs>    (optional, CATNAP neutralization panel)

The LANL/CATNAP protocol reports the databases as missing, naming the variable
that is unset or the path that does not exist, instead of failing silently.

4. **Install**:

.. code-block::

            scipion installp -p path_to_scipion-chem-immuno --devel

- **Tests**

To check the installation, simply run the following Scipion test:

===============
Buildbot status
===============

Status devel version: 

.. image:: http://scipion-test.cnb.csic.es:9980/badges/bioinformatics_dev.svg

Status production version: 

.. image:: http://scipion-test.cnb.csic.es:9980/badges/bioinformatics_prod.svg

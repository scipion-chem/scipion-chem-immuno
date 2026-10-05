================================
Immuno scipion plugin
================================

**Documentation under development, sorry for the inconvenience**

Scipion framework plugin for the use of immunoinformatics tools comming from several sources:

- IIITD: for epitope prediction and evaluation
- Vaxign-ML: for protective antigen prediction
- ScanNet and DiscoTope-3.0: for structure-based B-cell epitope prediction
- EpiDope: for sequence-based B-cell epitope prediction
- TMbed and SignalP-6.0: for transmembrane topology and signal peptide annotation
- IApred: for intrinsic antigenicity
- NetCleave: for C-terminal proteasomal cleavage
- StackGlyEmbed: for N-linked glycosylation sequons
- Epitope construct assembly: for building a multi-epitope construct from
  already selected candidates

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

- ScanNet, DiscoTope-3.0, TMbed, EpiDope, IApred, NetCleave and StackGlyEmbed:
These install automatically with their own conda environment, either as part of
the plugin installation or on demand:

.. code-block::

            scipion3 installb ScanNet DiscoTope TMbed epidope IApred NetCleave StackGlyEmbed

- SignalP-6.0:
This one cannot be installed automatically. DTU Health Tech does not allow
redistributing the package, so it has to be downloaded manually from
https://services.healthtech.dtu.dk/services/SignalP-6.0/ (requires an academic
account). Build a dedicated venv for it (Python 3.10, torch>1.7,<2, numpy<2)
and add these variables to the scipion.conf file:
    - SIGNALP_PYTHON_BIN = <path/to/venv/bin/python>  (interpreter that can import the signalp package)
    - SIGNALP_MODEL_DIR = <path/to/models>            (directory containing sequential_models_signalp6/)

The SignalP protocol reports the installation as missing instead of failing
silently if either variable is unset.

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

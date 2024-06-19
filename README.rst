================================
Immuno scipion plugin
================================

**Documentation under development, sorry for the inconvenience**

Scipion framework plugin for the use of immunoinformatics tools comming from several sources:

- IIITD: for epitope prediction and evaluation
- Vaxign-ML: for protective antigen prediction

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

3. **Define browser**

This plugin uses web browser to access the software servers online.
Therefore, you need to specify which browser to use and its location.
Do so editing the scipion.conf file and add the variables:
    - IIITD_BROWSER = firefox/chrome/chromium  (defines the browser to use, that must already be installed in your computer)
    - IIITD_BROWSER_PATH = <path/to/browser>   (defines the location of the binary for the browser use)

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

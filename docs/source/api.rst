Summary
=======

The idea
~~~~~~~~

The electronic properties of molecules can be greatly influenced by their accessible conformation.
As a result, it is essential to have a proper tool to explore the conformation space of small molecules.
Although conformation is usually explored at a classical level using force fields like UFF or MMFF94, these lack descriptions of quantum effects.
Conversely, starting from scratch at a high QM level of theory is not advisable due to two main reasons: the initial structure is crucial to ensure convergence,
and these methods are computationally expensive.

An intermediate level of theory, which is accessible for almost any current computer, is the semi-empirical level of theory.
``aleimi`` attempts to automate the combination of classical, semi-empirical, and QM level theories for the conformational analysis of small molecules.

``aleimi`` uses `RDKit <https://www.rdkit.org/docs/index.html>`_ for stochastic generation of conformer and classical optimization,
`MOPAC <https://github.com/openmopac/mopac>`_ for the semi-empirical part. And `Psi4 <https://psicode.org/>`_, `Gaussian <https://gaussian.com/>`_
or `Orca <https://www.orcasoftware.de/tutorials_orca/index.html#>`_ for the QM simulations. For the last part, ``aleimi`` only generate the files
to launch on the cluster.

General workflow
~~~~~~~~~~~~~~~~

Here we show an example workflow of ``aleimi`` that can be obtained with the following parameter file:

.. code-block:: yaml

    # aleimi.confgen.main
    numConfs: 1000
    rdkit_d_RMSD: 0.2
    UFF: false
    mopac_keywords: PM7 precise ef xyz geo-ok t=3h THREADS = 2 # You can also add `EPS=78.4` to take into account, implicitly, the solvent
    # aleimi.boltzmann.main
    Bd_rmsd: 1.0
    Bd_E: 0.0
    # aleimi.extractor.main
    energy_cut: 2
    engine: psi4

|workflow|


..  |workflow|  image:: https://github.com/ale94mleon/aleimi/blob/main/docs/source/_static/aleimi_workflow.png?raw=true
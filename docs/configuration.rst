.. _configuration:

Configuration Files
===================

SynIM uses **YAML** files to define the parameters of Adaptive Optics (AO) systems.
This standardizes the workflow and ensures seamless integration with the SPECULA framework.

Naming Conventions
------------------
The primary method SynIM uses to recognize various system components (such as wavefront sensors, deformable mirrors, optical sources, and the pupil) is based on the **section name** (the main key) within the YAML file. 

To ensure components are correctly identified and loaded, section names must follow these rules, where ``x`` represents a progressive number from 1 to N:

* **Shack-Hartmann Sensors (WFS):** Must be named ``sh`` (if there is only one sensor) or ``shx`` (e.g., ``sh1``, ``sh2``, or specifying the guide star like ``sh_lgs1``, ``sh_ngs1``, ``sh_ref1``).
* **Deformable Mirrors (DM):** Must be named ``dm`` (if single) or ``dmx`` (e.g., ``dm1``, ``dm2``).
* **Optical Sources:** Define the projection geometry (fundamental for tomographic matrices). They must be named ``sourcex`` or specify the type (e.g., ``source_lgs1``, ``source_ngs1``, ``source_ref1``, ``source_opt1``).
* **Pupil Mask:** Parameters for the pupil mask must be placed in a section named ``pupilstop`` or ``pupil_stop`` (general parameters like pixel sampling go in the ``main`` section instead).

Legacy PASSATA Support (.pro files)
-----------------------------------
Older workflows based on PASSATA used IDL-style configuration files (``.pro``). To maximize performance and keep the codebase clean and simple, **SynIM no longer performs native runtime parsing of .pro files**. 

However, full backward compatibility is maintained via a standalone offline conversion tool provided with the SynIM package.

Converting .pro to .yaml
~~~~~~~~~~~~~~~~~~~~~~~~
You can convert your legacy parameter files using the ``convert_pro_file.py`` file.

.. code-block:: bash

   python tools/convert_pro_file.py params_mcao.pro

This will generate a cleanly formatted ``params_mcao.yaml`` file in the same directory.
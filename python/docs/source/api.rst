.. _api_reference:

API Reference
=============

This is the auto-generated API reference for the ``libtts`` library.

Ground Detection Functions
--------------------------
.. rubric:: Ground Detection

.. automodule:: libtts.ground_detection
   :members:
   :exclude-members: main
   

Tree Detection Functions
------------------------
.. rubric:: Tree Detection

.. automodule:: libtts.tree_detection
   :members:
   :exclude-members: main


Tree Extraction Functions
-------------------------

These functions are implemented in C++ and wrapped with pybind11 for performance.
They are available in the ``_libtts`` submodule and wrapped in the main ``libtts`` library.

**Example Workflow**

.. code-block:: python

   import libtts

   # Example workflow for extracting a single tree
   as_file = libtts.generate_alpha_shape(tls_veg_file, th_alpha_sq)

   segfile = libtts.tls_extract_single_trees(
       as_file,
       treeloc_file,
       th_p2trunk_distance=0.2,
       th_search_radius=0.25
   )
   
   
.. automodule:: _libtts
   :members:
   :exclude-members: als_segment


Points Downsampling Functions
-----------------------------
.. rubric:: Points Downsampling

.. automodule:: libtts.points_downsampling
   :members:
   :exclude-members: main
   

Label Propagation Functions
---------------------------
.. rubric:: Label Propagation

.. automodule:: libtts.label_propagation
   :members:
   :exclude-members: label_pts_from_core_example
   

Dust Database
==============

Calculations are Currently Possible for the following dust particle compositions:

.. literalinclude:: ../dust_mie/optical_dat/dust_dict.yaml

..   :language: yaml

For example, this finds the complex indices of refraction for forsterite at 0.5 microns.

.. code-block:: python

   import dust_mie
   k,n = dust_mie.calc_mie.get_index_refrac(0.5,'Mg2SiO4')
   
You can access the dust dictionary the following way:

.. code-block:: python

   import dust_mie
   print(dust_mie.calc_mie.labelDict)

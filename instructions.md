# Instructions

You are to follow design_doc.md to write and alchemical free energy postprocessing library in rust. The first step should be to parse AMBER output files and integrate TI values (dH/dl) to produce a free energy estimate. After reading the design doc described in the following section, you should outline a plan and write it to outline.md in the current directory. 

# Resources

A design doc is present at design_doc.md. The git repositories for the python packages alchemlyb and pymbar have been cloned in the current directory. These are to be used as a resource for reference and testing. If you need to refer to the documentation for these, they can be accessed at https://alchemlyb.readthedocs.io/en/latest/ https://pymbar.readthedocs.io/en/stable/. Alchemlyb should be preferred in all ways, pymbar is only provided for refernce because alchemlyb uses it for some implementation details. If you need access to alchemlyb to generate ground truth testing values, a conda env has been set up. Access it with `conda activate alchemlyb`. If you want real-world data for testing, some can be found in `/gibbs/helmut/projects/flexible_ligands_GIST/mobley2008_test_set/diverse_representative_20/2-fluorophenol/centroid_TI_GIST/centroid_0/replicate_1/*/centroid_0.prod.out`. You should read this data only, not modify it. 

You have access to cargo, rustdoc, 
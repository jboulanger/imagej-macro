# Adjacency graph and shape

We define the rank between regions as the distance in the region adjacency graph and report a few shape parameters.

The input is a label image saved as a tif file with an active selection corresponding to the reference cell.

![image](https://github.com/jboulanger/imagej-macro/assets/3415561/f015b4c0-ff72-4341-9b35-f2507da6a048)

This macro uses MorpholibJ to compute the region adjacency graph (RAG). MorpholibJ can be installed by enabling the update site Help>upate>IJPB-plugins.

Legland, D., Arganda-Carreras, I., & Andrey, P. (2016). MorphoLibJ: integrated library and plugins for mathematical morphology with ImageJ. Bioinformatics, 32(22), 3532–3534. doi:10.1093/bioinformatics/btw413

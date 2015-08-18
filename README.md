# classify
Multiscale classification using Random Forest in R through GRASS GIS.  Tiles to manage large files.  Incorporates spectral and textural measures.

This script works in Grass6.x (hopefully).

The script uses training data from a vector in Grass to create a random forest
object.  It then creates a multilayer tiff file in Grass, which is tiled and fed
back into R for the prediction, resulting in the output of class tiles, which are
read back into Grass and mosaicked.

Mutliscaling is performed on all predictor layers, whose window sizes are stored
in @scales.  Textural measures are performed at the same scales - the measures
themselves are stored in %textures.  For each scale, aggregates are calculated
according to the functions in @aggregates.

Grass rasters to be used in the spectral analysis should be listed in a text file
called "spectral" in the same directory, and rasters used for the textural 
analysis should be stored in a file called "textural".  Note with the latter, 
the raster should be 8-bit integer (CELL type in Grass-speak).  Also, go easy on
the number of layers - they multiply right up by scales and measures, remember.

You need the R package, with the following packages installed:
randomForest
raster
rgdal
foreign

It's worth the pain - Random Forest is cool.  By the way, the script tiles the
multilayer output for the prediction in $tilesize x $tilesize pixel units, and the RF process
does the prediction in 100,000 pixel batches, both of which measures are designed
to prevent the process falling over no matter how large the raster files...
   ...however...
that's not to say that if your files are too large, the proces won't be in-
credibly slow, so take care.

What ever your region extents are at the time of invoking the script, those are
the maximum extents that will be used (and tiled within).  Set your mask carefully,
it will be maintained.

Cheers

Damien.O'Grady
11 August 2015

# CENITH 0.0.99.1

* new features
add vignette - tutorial part 1 treseg 

# CENITH 0.0.99.0
Prerelease version 0.99

* new features
CENITH - package - add package script with imports and descriptions
add Example Data for Crossvalidation

* bugfixes
some major changes in BestSegVal cat codes.

# CENITH 0.0.0.97
* new features
complete clean up of Development functions (keep stable version of BestSegVal)

TreeSegVal - now iterates over a,b,h,MIN, and chm filters. Advanced ETA estimation added (now uses full iterations of TreeSeg). Further cleand up "cat" codes to more easy see the ETA.

Added all new Example data.

# CENITH 0.0.0.96
* new features
TreeSegValEXPRMTL - advanced version of TreeSegVal which additionally iterates over MIN and supports filtering for chm.

# CENITH 0.0.0.95
* bugfixes
TreeSegVal returns correct values for "Segment quality"

* new features
TreeSegCV - function to perform an n fold cross validation. Uses values (estimated with TreeSegVal recommended) to test the quality on other sides to receive a quality estimation of the inout values for the Segmentation in the AOI.

# CENITH 0.0.0.94
Update for MIN and MAX

* added features
TreeSeg -  Now returns Error, if after clipping no Segments are left. This is used to handle MIN and MAX in TreeSegVal.
TreeSegVal - now uses MIN and MAX in TreeSeg
           - improved results, including hit, oversegmentation and undersegmentation rates and absolut values.

# Cenith 0.0.0.93
Update

* updated features
BestSegValBETA - deleted older check for to high values. To estimate the time to target now up to 3 random iterations are used to sto time. if all 3 fail the function stops. added example with 3/4 iterations to fail.
*TreeSeg -  added MIN MAX from BestSegValBETA to basic Segmentation function. Crops polygons < MIN / >MAX

# CENITH 0.0.0.92
* development version for Validation of Segments

* new feature
BestSegValBETA - development verison contains error catching

# CENITH 0.0.0.91

* changes for development
TreeSeg - now returns Treepos + segments
BestSegVal - now computes more and other validation results

# CENITH 0.0.0.9000
initial version

*new feature

TreeSeg - function for TreeCrown Segmentation
BestSegVal - function to iterate over several a,b,h values to detect best fitting parameters for TreeSeg

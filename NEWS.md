# CENITH 0.1.2
patch version

* bugfixes
With a new rdgal version now the strings for the vector example data contain additional +towgs84=0,0,0,0,0,0,0.
With those differnet strings (even if projcetion is still teh same) some functions and the vignette would not work.
Fixed by setting the crs of the vectors to crs of rasters by hand in example and vignette code.

* update vignette
vignettes now come in hmtl format for direct download via github.
Install without "buildvignette=T" is now recommended.

# CENITH 0.1.1
patch version

* bugfixes general
For all functions using 'TreeSeg' - In TreeSeg now the window diameter is not set to a maximum. This leads to no more errors is the windiameter is wider.
Correct name for CENITH-Package.R in help is used.
Set depending packages in description.

* bugfixes and updates
TreeSeg - default MaxWinDiameter in subfucntion is now set to NULL. 
TreeSeg - corrected example added. Added require() to example.
BestSegVal - now 3rd try in ETA check works. Further now continues even if all 3 trys do not work.
BestSegVal - spelling corrections and added some cat() codes.
BestSegVal - add require() to example.
TreeSegCV - corrections in help
TreeSegCV - spelling corrections

# CENITH 0.1.0
* release version for GitHub (not rdy for CRAN)

# CNEITH 0.0.99.5
* spelling corrections

correction of spelling for all functions.
some added links to other functions.


# CENITH 0.0.99.4
* bugfixes

fixed TreeSegCV - now if a combination leads to an error the error is printen correctly and the function continues.

# CENITH 0.0.99.3
* new features
update vignette - add result='hide' to chunks for increased overview in tutorial.

# CENITH 0.0.99.2

* new features
update vignette - add part 2 (BestSegVal) and 3 (TreeSegCV) to vignette
add example data - chm and rgb for the AOI

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

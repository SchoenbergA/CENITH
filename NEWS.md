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


The files in this directory constitute a C++ implementation of a
region-specific dendritic pruning algorithm that is presented in:

Rishikesh Narayanan, Anusha Narayan, Sumatra Chattarji. A probabilistic
framework for region-specific remodeling of dendrites in three-dimensional
neuronal reconstructions. Neural Comput. 2005 Jan;17(1):75-96.

The specificity of the region can be specified using experimental
data, and the algorithm will account for that using a ratio-of-Gaussian
distribution. For more details, see the above paper. All details regarding
the programs, and the methodology to run them with appropriate inputs are
mentioned below.

Implemented by Rishikesh Narayanan. Contact: rishi.n@gmail.com

/**********************************************************************/

The Makefile and the C++ files in this directory can be used to make two
different executables: "Sholl" and "Prune". 

Program "Sholl"
===============

This program will perform Sholl analysis on an input neuronal morphology,
specified in the SWC format. Usage:

Sholl <Input SWC File> <Output Filename> 

The Output Filename is actually a base filename. Two files will be saved
the same basename, but with two different extensions: .asl and .bsl. The
.asl and .bsl files contain the Sholl analyses for the apical and the basal
dendrites, respectively. For example, when you run the following command,
the console will have the following output.

Sholl Input/CA3b4.swc CA3B4

Centroid of Soma is: -2.70697  13.2238  1.06853
Total dendritic length is: 13669
Total number of BP is: 63

There will be two files CA3B4.asl and CA3B4.bsl saved. These files
will have the dendritic length in each of the Sholl segments (will
be automatically computed based on the annulus length, currently
set at 50 micron in Sholl.h), followed by the number of branch
points in each of these Sholl segments (asl: apical; bsl: basal).

Program "Prune"
===============

This program will prune an input neuronal morphology, specified in
the SWC format, in a region specific manner. The region specificity
can be specified by experimental data as the means and variances
of dendritic length and branch points in control and "treated"
animals. These have to be specified at different Sholl segments,
and the Sholl segment annulus size in the experiments should match
with the Sholl segment annulus size used in the program.  Provided
is an example where mean and variance of dendritic length and branch
points (across 8 apical and 6 basal Sholl segments, with annulus
size at 50 microns) from CA3 neurons are taken from control and
chronic stressed animals (see Vyas et al., below). These data are
then used to generate the probability of pruning dendritic length
and branch points in each Sholl segment, and a pruning algorithm
prunes and saves the data into Output SWC files. For more details
on the algorithm, see the paper referenced above.

Usage: 

Prune <PRN file>

The PRN file contains the links to the files containing the morphology that
is to be pruned, and the files that containing the aforementioned
statistics. The format for the PRN file is given below.

When you run "Prune in.prn", the console will look as follows:

Centroid of Soma is: -2.70697  13.2238  1.06853
Total dendritic length is: 13669
No of Soma pts: 32
No of Stems: 4

Saving file Output/C4_1000.swc after pruning 1015.4 micron of length

Saving file Output/C4_2000.swc after pruning 2004.27 micron of length

Saving file Output/C4_3000.swc after pruning 3006.8 micron of length

Saving file Output/C4_4000.swc after pruning 4009.38 micron of length

Saving file Output/C4_5000.swc after pruning 5005.99 micron of length

Saving file Output/C4_6000.swc after pruning 6006.98 micron of length

Saving file Output/C4_7000.swc after pruning 7001.05 micron of length

Saving file Output/C4_8000.swc after pruning 8000.71 micron of length

Saving file Output/C4_9000.swc after pruning 9002.32 micron of length

Saving file Output/C4_10000.swc after pruning 10005.2 micron of length

Final save in Output/C4_final.swc after pruning 10500.5 micron of length

The output files with the pruned dendritic tree in SWC format may be found
in the filenames that the program mentions (which can be specified in the
.prn file, as mentioned below). These SWC outputs can be converted to
Neuron's HOC code using Neuron. An example output which was generated from
C4_final.swc is appended as C4_final.hoc. The HOC file corresponding to the
unpruned SWC file (Input/CA3B4.swc) is also appended as CA3B4.hoc.


File format for .prn files 
==========================

1. SWC file containing the morphology that should be pruned.

2. Base filename for files containing statistics which should be
used for pruning different Sholl segments. This base filename will
be attached with .cb, .sb, .cd, and .sd; these four files will have
the statistics for pruning branch points and dendritic lengths (see
below for their formats).

3. Number of apical Sholl segments that data in 2 above is available
(this would be the number of apical values that should be provided
in the files mentioned in 2)

4. Number of basal Sholl segments that data in 2 above is available
(this would be the number of basal values that should be provided
in the files mentioned in 2)

5. Dendritic length at which pruning should be stopped.

6. Output basefilename. SWC output files, pruned as per the statistics
specified by files listed in 2 above, will be saved every 1000
microns of dendritic pruning. A number will be attached to this
base filename, and the the SWC file will be saved in that filename.
A final filename named "basefilename"_final.swc will also be saved
after pruning to the number specified in 5 is complete.

An example in.prn file is present in this directory.

/**********************************************************************/

File format for .cb and .sb files 
=================================

These files give the branching point statistics in apical and basal
dendrites in that order. .cb corresponds to branching point statistics
from control animals, and .sb corresponds to branching point
statistics from stress animals. The format is:

Mean value of variation in BP with stress in each apical segment
Variance of variation in BP with stress in each apical segment

Mean value of variation in BP with stress in each basal segment
Variance of variation in BP with stress in each basal segment

Example IL.cb and IL.sb files are appended to the Input subdirectory,
and are used in the in.prn file mentioned above. These data are
from Figure 1 of Vyas et al., J. Neuroscience, 22(15):6810-6818,
2002.  

/**********************************************************************/

File format for .cd and .sd files 
=================================

These files give the dendritic length statistics in apical and basal
dendrites in that order. .cd corresponds to dendritic length
statistics from control animals, and .sd corresponds to dendritic
length statistics from stress animals. The format is:

Mean value of variation in DL with stress in each apical segments
Variance of variation in DL with stress in each apical segments

Mean value of variation in DL with stress in each basal segments
Variance of variation in DL with stress in each basal segments

Example IL.cd and IL.sd files are appended to the Input subdirectory,
and are used in the in.prn file mentioned above. These data are
from Figure 1 of Vyas et al., J. Neuroscience, 22(15):6810-6818,
2002.  

If you don't want to prune the apical side, say, then set all the
apical means and variances in all four statistics files (.cd, .sd,
.cb, .sb) to zero, without altering the basal-side statistics. A 
similar method would hold for not pruning the basal side as well.

/**********************************************************************/

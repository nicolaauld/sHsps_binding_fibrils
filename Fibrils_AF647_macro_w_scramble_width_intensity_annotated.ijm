macro "Alignment macro [h]" {
	
// This macro is suited to imaging fibrils with ASCP dye (emission in 647 channel) and an AF647-labelled protein.
// The experiments were conducted with alternating laser excitation (rather than both at once) so there should be
// 

// No bead image alignment is performed !!! Need a different script with alignment for AF488-labelled protein

// So firstly the cropping, background correction are performed on the fibril image and then the Hsp image. 
// Then the fibril contours/ridges are identified  using the Ridge Detection plugin (which also measures the width) 
// and the intensity of each fibril is measured from the fibril ROIs. 

// Next, a peak-finding step is performed and peaks are saved.

// This is followed by the colocalisation section where the coordinates of the peaks that are found in the fibril and
// Hsp images are compared to identify those that are colocalised with one another. These colocalisations are saved as well

// Then we move to the photobleaching trajectory sections. Here, we open the movies again and the peaks are imported
// into the ROI manager. Then the fluorescence intensity of the peaks over the course of the movie is measured
// and saved into a file. This is done firstly for the colocalised peaks and then the non-colocalised peaks.

// After this I want to get an estimate of the degree of random colocalisation that occurs when two spots on the 
// coverslip are overlapping because of random binding to the same spot on the coverslip, rather than actually being
// in complex with one another. To do this, I use a little script that randomised the coordinates of the sHsp peaks
// and randomises within the maximum and minimum x and y coordinates of the real peaks. Then these randomised
// coordinates are checked for colocalisation with the real coordinates of the fibrils.

// ---------------------------------------IMPORTANT NOTE:------------------------------------------------------

// This script requires the extra files
// "count_colocalized_peaks_Fibril.js"
// "s_m_b.jar" plugin called "Single Molecule Biophysics"
// "Filter_for_colocal.ijm"
// "Filter_for_noncolocal.ijm"
// "results_to_ROI.txt"
// "peak_intensity.txt"
// "scramble_xy.js"

// Some of these files will live in the "plugins" folder" : "count_colocalised_peaks_Andrew.js", "scramble_xy.js"

//-------------------------------------------------------------------------------------------------------------
// Folder setup and manual selection of fibril and Hsp images
//-------------------------------------------------------------------------------------------------------------

var path = File.openDialog ("Please select fibril RAW image");
var dir = File.getParent(path);
var name = File.getName(path);
var upper_dir = File.getParent(dir);

splitter = split(dir, File.separator);
Array.print(splitter);
len = splitter.length ;
trimm = len -1;
ROIarr = Array.slice(splitter,0,trimm);
Array.print(ROIarr);

str = ""
filesep = File.separator;
     for (i=0; i<ROIarr.length-1; i++) 
         str = str + ROIarr[i] + filesep; 
     str = str + ROIarr[ROIarr.length-1] + filesep; 
ROI = filesep + str + "ROI.zip"; //the path for the script to follow to find its' way to the ROIs
print("Raw Image path= " + path)
print("ROI path= " + ROI); 
output = dir + filesep 
print("Output path= " + output)

Hsp_results = output + "Hsp_results.csv";
fibril_results = output + "Fibril_results.csv";
Hsp_image = output + "Hsp.tif";

//SELECT AF647 IMAGE NEXT
//read in the image with the AF647 trajectories within them
var AF647_imagepath = File.openDialog ("Please select Hsp AF647 image");
var dir = File.getParent(AF647_imagepath);
var name = File.getName(AF647_imagepath);

splitter = split(dir, File.separator) 
Array.print(splitter)
len = splitter.length 
trimm = len -1 
ROIarr = Array.slice(splitter,0,trimm) 
Array.print(ROIarr)

str = ""
AF647_filesep = File.separator;
     for (i=0; i<ROIarr.length-1; i++) 
         str = str + ROIarr[i] + filesep; 
     str = str + ROIarr[ROIarr.length-1] + filesep; 
ROI = filesep + str + "ROI.zip"; 
print("af647 Image path= " + AF647_imagepath)
print("ROI path= " + ROI); 
output_647 = dir + filesep 
print("Output path= " + output_647)
run("Set Measurements...", "centroid stack redirect=None decimal=3");

//-------------------------------------------------------------------------------------------------------------
// Defining background images and also defining paths to save the output files
//-------------------------------------------------------------------------------------------------------------

AF647_results = output_647 + "AF647_results.csv";
AF647_image = output_647 + "MAX_Background_corrected_AF647.tif";
bg_AF647= output_647+ "Background_corrected_AF647.tif"

Hsp_background_image = filesep + str + "AF647_background.tif"
fib_background_image= filesep + str + "fibril_background.tif"

fibril = output + "fibril_MAX.tif";

//create output folder
Fibcolocalpath = filesep + str + "/Colocalisation_analysis/";
print(Fibcolocalpath);
File.makeDirectory(Fibcolocalpath);

Fibrilextrapath = upper_dir + File.separator + "Extra_fibril_analysis" + File.separator;
print("Corrected Extra_fibril_analysis Path: " + Fibrilextrapath);
File.makeDirectory(Fibrilextrapath);

split_path = split(path, filesep);  
last_part = split_path[split_path.length - 2];  // Move up one level to get the folder number

// Check if it's a number and store it
if (matches(last_part, "^[0-9]+$")) { 
    image_number = last_part;
} else {
    print("No number found at the end of the path");
}


//-------------------------------------------------------------------------------------------------------------
// Fibril channel image background correction and max projection
//-------------------------------------------------------------------------------------------------------------

//now start actual analysis here
// Fibril channel
open(path); //select raw image
Imagetitle =getTitle();
rename("Raw"); 

roiManager("Open", ROI); //open ROI.zip where you defined channels and saved the file.
roiManager("Select", 0);	//this will be the fibril channel
run("Duplicate...", "title=fibril duplicate");

//Background correction 
open(fib_background_image);
rename("Background_fib");
selectWindow("fibril");
run("32-bit");
run("Beam Profile Correction", "electronic_offset=700 background_image=Background_fib");
setMinAndMax(0, 65536);
run("16-bit");

roiManager("Select", 1); //because these two are the same channel, the crop to fix the edges is called 1, not 2
run("Duplicate...", "title=fibril_background_corrected duplicate");
saveAs("tif", output + "Background_corrected_fib.tif");
selectWindow("Background_fib");
close();

//this will save the max fibril intensity to find fibril ridges
selectWindow("Background_corrected_fib.tif");
run("Z Project...", "stop=20 projection=[Max Intensity]");
run("Properties...", "channels=1 slices=1 frames=1 unit=pixels pixel_width=1 pixel_height=1 voxel_depth=1.0000000");
saveAs("tif", output + "fibril_MAX.tif");
close();
roiManager("reset");
selectWindow("Background_corrected_fib.tif");
close();

//-------------------------------------------------------------------------------------------------------------
// Hsp channel image background correction and max projection
//-------------------------------------------------------------------------------------------------------------



//Hsp channel- don't need to align this because they are in the same channel but separate images
//(i.e. alternate excitation)
open(AF647_imagepath); //select raw image
Imagetitle =getTitle();
Imagetitle =getTitle();
rename("AF647_image"); 

roiManager("Open", ROI);
roiManager("Select", 0);	//this will be the AF647 channel
run("Duplicate...", "title=AF647 duplicate");

//Background correction 
open(Hsp_background_image);
rename("Background_AF647");
selectWindow("AF647");
run("32-bit");
run("Beam Profile Correction", "electronic_offset=550 background_image=Background_AF647");
setMinAndMax(0, 65536);
run("16-bit");
//crop edges again
roiManager("Select", 1);
run("Duplicate...", "title=AF647_background_corrected duplicate");
saveAs("tif", output_647 + "Background_corrected_AF647.tif");
selectWindow("Background_AF647");
close();

//creat max intensity file for peaks to be identified in AF647 channel
list = getFileList(dir); 
selectWindow("Background_corrected_AF647.tif");
run("Z Project...", "stop=20 projection=[Max Intensity] ");
run("Properties...", "channels=1 slices=1 frames=1 unit=pixels pixel_width=1 pixel_height=1 voxel_depth=1.0000000");
saveAs("tiff", output_647+getTitle()); 
selectWindow("MAX_Background_corrected_AF647.tif");
close();
selectWindow("Background_corrected_AF647.tif");
close();

roiManager("reset");

//-------------------------------------------------------------------------------------------------------------
// Fibril ridge detection step
//-------------------------------------------------------------------------------------------------------------
open(fibril);
run("Enhance Contrast", "saturated=0.45");
run("8-bit");
run("Ridge Detection", "line_width=3.50 high_contrast=230 low_contrast=87 estimate_width extend_line displayresults add_to_manager method_for_overlap_resolution=NONE sigma=1.2 lower_threshold=10 upper_threshold=11 minimum_line_length=4 maximum=0");
waitForUser("Delete the junctions from the ROI manager before pressing OK");
selectWindow("fibril_MAX.tif");
	saveAs("tif", output + "Ridge_detection.tif");
	close();

selectWindow("Summary");
//save length data and x y coordinates of peaks
saveAs("Results", output + "Length_data.csv");
run("Close");
selectWindow("Results");
saveAs("Results", output + "Fibril_results.csv");
run("Close");
selectWindow("Junctions");
run("Close");
print("\\Clear");

//-------------------------------------------------------------------------------------------------------------
// Fibril intensity measurements
//-------------------------------------------------------------------------------------------------------------

n = roiManager("Count");  // Count total ROIs in the manager
roiManager("Deselect");   // Deselect any previously selected ROI
run("Set Measurements...", "area mean standard min centroid integrated redirect=None decimal=3");

selectWindow("fibril");
// Loop through and measure each ROI individually
for (i = 0; i < n; i++) {
    roiManager("Select", i);
    roiManager("Measure");
}

// in the results Mean refers to the mean intensity

selectWindow("Results");
saveAs("Results", Fibrilextrapath + "#" + image_number + "_Fibril_results_w_intensity.csv");
run("Close");

roiManager("reset");


//-------------------------------------------------------------------------------------------------------------
// Peak selection in Hsp image
//-------------------------------------------------------------------------------------------------------------

run("Set Measurements...", "centroid stack redirect=None decimal=3");

open(AF647_image);
rename("Hsp");
//Find the peaks within the Hsp image
run("Enhance Contrast", "saturated=0.35");
run("Clear Results");
run("Peak Finder");// need to check the threshold here
waitForUser("Have you selected your peaks? \n \nPress OK to continue");
roiManager("Measure");	//this will allow for a table to be produced
saveAs("Results", output_647 + "Hsp_results.csv"); 
close();
selectWindow("Results"); 
     run("Close"); 
roiManager("reset");


//-------------------------------------------------------------------------------------------------------------
// Colocalisation analysis
//-------------------------------------------------------------------------------------------------------------

open(fibril_results);
IJ.renameResults("Fibril");
open("Hsp_results.csv");
IJ.renameResults("Hsp");
run("count colocalized peaks Fibril", "table_1=Fibril table_2=Hsp maximum_distance=3");
selectWindow("Colocalization");
saveAs("tiff", output + "Colocalization.tiff");
close();
selectWindow("Fibril");
saveAs("Results", output + "Fibril_colocalisation.csv"); 
saveAs("Results", Fibcolocalpath + Imagetitle + ".csv");
IJ.renameResults("Results");
selectWindow("Hsp");
		run("Close");
//now filter these data for only those colocalised not all the spots combined.
runMacro("C:/Nicola/UNIVERSITY/PhD/SM/ImageJ macros/Filter_for_colocal.ijm");
IJ.renameResults("Fibril_colocalisation_only");
selectWindow("Fibril_colocalisation_only");
saveAs("Results", output + "Fibril_colocalisation_only.csv"); 
saveAs("Results", Fibcolocalpath + Imagetitle + "_colocalized.csv");
run("Close");
roiManager("reset");


//-------------------------------------------------------------------------------------------------------------
// Trajectory analysis for colocalised molecules
//-------------------------------------------------------------------------------------------------------------

//create output folder for trajectories that are colocalised (only)
coloc_trajectory_output = filesep + str + "/coloc_trajectories/";
print(coloc_trajectory_output);
File.makeDirectory(coloc_trajectory_output);

//now find the coordinates of the hsps colocalised with fibrils, and get the trajectories at these coordinates
open(bg_AF647);
open(output + "Fibril_colocalisation_only.csv");
Table.rename("Fibril_colocalisation_only.csv", "Results");
runMacro("C:/Nicola/UNIVERSITY/PhD/SM/LRplugins/plugins/Macros/v2_results_to_ROI_X2.txt");
selectWindow("Results");
        run("Close");
runMacro("C:/Nicola/UNIVERSITY/PhD/SM/LRplugins/plugins/Macros/peak_intensity.txt"); // you need to put the location of your script
selectWindow("Results");
        saveAs("Results", coloc_trajectory_output + Imagetitle + "_HSP_colocal_traj.csv");
        run("Close");
roiManager("reset");

selectWindow("Raw");
		run("Close");
selectWindow("fibril");
		run("Close");

selectWindow("AF647");
		run("Close");
selectWindow("AF647_image");
		run("Close");	

//-------------------------------------------------------------------------------------------------------------
// Trajectory analysis for NON-colocalised molecules
//-------------------------------------------------------------------------------------------------------------

//filter Non-colocal tables
noncoloc_output_folder = filesep + str  + "/non-coloc_trajectories/";
print(noncoloc_output_folder);
File.makeDirectory(noncoloc_output_folder);

hsp_noncoloc_trajectories = noncoloc_output_folder + "/647/";
print(hsp_noncoloc_trajectories);
File.makeDirectory(hsp_noncoloc_trajectories);

//Hsp
open(fibril_results);
IJ.renameResults("fibril");
open(output_647 + "Hsp_results.csv");
IJ.renameResults("Hsp");
run("count colocalized peaks Nicola", "table_1=Hsp table_2=fibril maximum_distance=3");
selectWindow("Colocalization");
close();

selectWindow("fibril");
		run("Close");
selectWindow("Hsp");
IJ.renameResults("Results");
runMacro("C:/Nicola/UNIVERSITY/PhD/SM/ImageJ macros/Filter_for_noncolocal.ijm");
IJ.renameResults("Hsp_non-colocalisation");
selectWindow("Hsp_non-colocalisation");
saveAs("Results", output + "647_non-colocalisation.csv"); 
selectWindow("647_non-colocalisation.csv");
     run("Close"); 

open(bg_AF647);	
open(output + "647_non-colocalisation.csv");
IJ.renameResults("Results");
runMacro("C:/Nicola/UNIVERSITY/PhD/SM/LRplugins/plugins/Macros/results_to_ROI.txt"); // you need to put the location of your script
selectWindow("Results");
        run("Close");
runMacro("C:/Nicola/UNIVERSITY/PhD/SM/LRplugins/plugins/Macros/peak_intensity.txt"); // you need to put the location of your script
//selectWindow("Results");
        saveAs("Results", output + "Hsp_non-colocal_traj.csv");
selectWindow("Results");
     	saveAs("Results", hsp_noncoloc_trajectories + Imagetitle + "_Hsp_non-colocal_traj.csv");
        run("Close");
selectWindow("Background_corrected_AF647.tif");
        run("Close");
selectWindow("Background_corrected_AF647-1.tif");
        run("Close");
roiManager("reset");

//-------------------------------------------------------------------------------------------------------------
// Scrambling section where the coordinates of the Hsp foci are randomised then checked for colocalisation
//-------------------------------------------------------------------------------------------------------------

// first make a folder for this stuff to go into
Fibscramblepath = filesep + str + "/Scrambled_fibril_colocalisation_analysis/";
print(Fibscramblepath);
File.makeDirectory(Fibscramblepath);

open(fibril_results);
IJ.renameResults("Fibril");
open("Hsp_results.csv");
IJ.renameResults("Hsp");

// Load the results table from a CSV file
open("Hsp_results.csv");

// Rename the results table
IJ.renameResults("Hsp");

// In this step we are actually doing the scrambling/randomising
runMacro("C:/Nicola/UNIVERSITY/PhD/SM/ImageJ macros/scramble_xy.js");
selectWindow("Hsp");
saveAs("Results", Fibscramblepath + Imagetitle + "_Hsp_scrambled_results.csv");
name_of_scrambled_csv = Imagetitle + "_Hsp_scrambled_results.csv";

//now running the colocalisation macro plugin thing
run("count colocalized peaks Fibril", "table_1=Fibril table_2=Hsp maximum_distance=3");
selectWindow("Colocalization");
close();
selectWindow("Fibril");
//saveAs("Results", output + "Fibril_colocalisation.csv"); 
saveAs("Results", Fibscramblepath + Imagetitle + "_Fibril_scramble_Hsp_results.csv");
IJ.renameResults("Results");
selectWindow("Hsp");
		run("Close");
//now filter these data for only those colocalised not all the spots combined.
runMacro("C:/Nicola/UNIVERSITY/PhD/SM/ImageJ macros/Filter_for_colocal.ijm");
IJ.renameResults("Fibril_colocalisation_only");
selectWindow("Fibril_colocalisation_only");

saveAs("Results", Fibscramblepath + Imagetitle + "_scramble_colocalized.csv");
run("Close");

selectWindow(name_of_scrambled_csv);
run("Close");
roiManager("reset");

//-------------------------------------------------------------------------------------------------------------
// Saving things in a different folder for ease of access
//-------------------------------------------------------------------------------------------------------------

// Lastly just saving the normal Hsp results to the scrambled folder so they are easier to access later

open("Hsp_results.csv");
IJ.renameResults("Hsp");
saveAs("Results", Fibscramblepath + Imagetitle + "_Hsp_results.csv");
run("Close");

//-------------------------------------------------------------------------------------------------------------
// End
//-------------------------------------------------------------------------------------------------------------

waitForUser("All done! :) \n \nPress OK to continue");
}
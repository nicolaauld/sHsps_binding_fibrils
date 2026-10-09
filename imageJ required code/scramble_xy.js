// Import necessary classes from ImageJ
importClass(Packages.ij.IJ);
importClass(Packages.ij.measure.ResultsTable);
importClass(java.util.Random);
importClass(java.lang.Double); // Import Double class explicitly


// Get the Results table by its new name
var rt = ResultsTable.getResultsTable("Hsp");

// Check if the table exists
if (rt == null) {
    IJ.log("No results table found with the name 'Hsp'.");
} else {
    // Get the number of rows in the table
    var nRows = rt.getCounter();

    // Initialize min and max values
    var minX = Double.POSITIVE_INFINITY;
    var maxX = Double.NEGATIVE_INFINITY;
    var minY = Double.POSITIVE_INFINITY;
    var maxY = Double.NEGATIVE_INFINITY;

    // Find the min and max values for X and Y
    for (var i = 0; i < nRows; i++) {
        var x = rt.getValue("X", i);
        var y = rt.getValue("Y", i);

        if (x < minX) minX = x;
        if (x > maxX) maxX = x;
        if (y < minY) minY = y;
        if (y > maxY) maxY = y;
    }

    // Create a Random object for generating random values
    var rand = new Random();

    // Randomize X and Y coordinates within the found min and max values
    for (var i = 0; i < nRows; i++) {
        // Generate random values within the range
        var newX = minX + rand.nextDouble() * (maxX - minX);
        var newY = minY + rand.nextDouble() * (maxY - minY);

        // Update the table with new random values
        rt.setValue("X", i, newX);
        rt.setValue("Y", i, newY);
    }

    // Display the updated Results table
    rt.show("Hsp");

    // Ensure the table is updated in ImageJ
    IJ.setTextPanel(rt);
}

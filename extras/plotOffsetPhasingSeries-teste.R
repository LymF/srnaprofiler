# R code for plotOffsetPhasingSeries.R

# Load command-line arguments
# `trailingOnly = TRUE` ensures that only arguments after `--args` (if Rscript is used) are captured.
args <- commandArgs(trailingOnly = TRUE)

# --- Debug Section (Uncomment to see received arguments) ---
# print("Arguments received by the R script:")
# print(args)
# cat("Number of arguments:", length(args), "\n")
# ----------------------------------------------------------

# Check if the correct number of arguments has been provided
if (length(args) < 2) {
  stop("Usage: Rscript plotOffsetPhasingSeries.R <input_offset_file_path> <output_pdf_prefix>\n",
       "Example: Rscript plotOffsetPhasingSeries.R offset-piRNA.offset offset-piRNA")
}

# Assign arguments to more descriptive variables
input_file <- args[1]
output_prefix <- args[2]

# --- Debug Section ---
# cat("Input file:", input_file, "\n")
# cat("Output prefix:", output_prefix, "\n")
# ---------------------

# Read the offset data file
# Assumes the file is tab-delimited and has a header
if (!file.exists(input_file)) {
  stop(paste("Error: Input file not found:", input_file))
}

off <- tryCatch({
  read.csv(input_file, sep="\t", header=TRUE)
}, error = function(e) {
  stop(paste("Error reading CSV file:", e$message, "\nPlease check file format:", input_file))
})

# --- Debug Section ---
# print("First few rows of the read data:")
# print(head(off))
# print("Structure of the read data:")
# str(off)
# ---------------------

# Open the PDF graphics device
pdf_filename <- paste(output_prefix, ".pdf", sep = "")
pdf(pdf_filename, width=5, height=4)

# Ensure the graphics device is closed even if an error occurs
on.exit(dev.off())

# Check if there is data to plot
if (nrow(off) > 0 && "Count" %in% colnames(off)) {
  
  # Check if the 'Count' column has non-zero values for scaling
  if (all(off$Count == 0)) {
    warning("All values in the 'Count' column are zero. The graph may appear flat or empty.")
    v_off <- rep(0, nrow(off)) # Set all scaled values to 0
  } else {
    v_off <- scale(off$Count) # Scale the data (center and divide by standard deviation)
  }
  
  # --- Debug Section ---
  # cat("Scaled values (v_off):\n")
  # print(head(v_off))
  # ---------------------
  
  # Define Y limits for the plot dynamically
  # Ensure at least a small range is visible
  y_min <- min(v_off, -3, na.rm = TRUE) # Include -3 to ensure visibility of abline
  y_max <- max(v_off, 2, na.rm = TRUE)  # Include 2 to ensure visibility of abline
  
  # Adjust ylim to handle cases where y_min and y_max are too close (e.g., all values are 0)
  if (y_max - y_min < 1) { # If the range is too small, expand it
    y_center <- (y_max + y_min) / 2
    y_min <- y_center - 0.5
    y_max <- y_center + 0.5
  }
  
  # Create the main plot
  plot(off$Overlap, v_off,
       ylim=c(y_min, y_max),         # Dynamic Y limits
       xaxp=c(-30,30,30),            # Define 30 ticks on the X axis from -30 to 30
       type="l",                     # Line type
       lty=1,                        # Line style (solid)
       col="black",                  # Line color
       cex.lab=0.6,                  # Font size of axis labels
       font.lab=1,                   # Font style of axis labels (normal)
       cex.axis=0.6,                 # Font size of tick labels
       xlab="Distance from 5' end (nt)",  # X-axis label
       ylab="Relative freq. of antisense 5' end", # Y-axis label
       main="Phasing Series Offset") # Plot title
  
  # Add vertical grid lines
  abline(v=seq(-29,29,by=2), lty=2, lwd=0.5, col="lightgray")
  abline(v=c(18), lty=2, col="blue", lwd=0.5) # Vertical line at 18
  abline(h=c(2), lty=2, col="red", lwd=0.5)   # Horizontal line at 2
  
} else {
  # If no data or 'Count' column is missing
  warning("The input data frame 'off' is empty or does not contain the 'Count' column. Generating an empty plot.")
  plot(1, type="n", xlim=c(-30, 30), ylim=c(-3, 8),
       xlab="Distance from 5' end (nt)",
       ylab="Relative freq. of antisense 5' end",
       main="Phasing Series Offset (No Data)")
  text(0, 0, "No data to plot or invalid data.", cex=0.8, col="red")
}

# Close the PDF device (on.exit also ensures this)
dev.off()

# End of script
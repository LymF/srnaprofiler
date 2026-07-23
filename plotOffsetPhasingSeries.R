# --- SCRIPT START: DIAGNOSTIC SECTION ---
# This section helps debug argument-passing issues.
cat("--- R SCRIPT STARTING ---\n")
raw_args <- commandArgs(trailingOnly = TRUE)
cat("Raw arguments received by R:", paste(raw_args, collapse=" | "), "\n")

# Filter out any arguments that start with "--" to avoid issues with
# system-injected flags like --no-restore or --no-save.
args <- raw_args[!startsWith(raw_args, "--")]
cat("Arguments after filtering system flags:", paste(args, collapse=" | "), "\n")
# --- END OF DIAGNOSTIC SECTION ---


# Check if the correct number of arguments remains after filtering.
if (length(args) != 4) {
  stop(paste("Error: 4 arguments expected, but", length(args), "received after filtering. Got:", paste(args, collapse=" ")), call. = FALSE)
}

# Define input file names and output prefix from the filtered arguments
phasing_plus_file <- args[1]
phasing_minus_file <- args[2]
offset_file <- args[3]
output_prefix <- args[4]

cat("Input file 1 (plus):", phasing_plus_file, "\n")
cat("Input file 2 (minus):", phasing_minus_file, "\n")
cat("Input file 3 (offset):", offset_file, "\n")
cat("Output prefix:", output_prefix, "\n")


# Define the output PDF file name
pdf_output_file <- paste(output_prefix, ".phasing_offset_line.pdf", sep = "")
cat("PDF will be saved to:", pdf_output_file, "\n")

# --- Start PDF device ---
pdf(pdf_output_file, width = 2, height = 7)

# Set plot layout to 3 rows, 1 column, and adjust margins
par(mfrow = c(3, 1), mar = c(4.5, 4.5, 3, 2)) 

# --- Function to create each plot ---
create_line_plot <- function(file_path, title, color, show_ref_line = FALSE) {
  if (file.exists(file_path) && file.info(file_path)$size > 0) {
    # Read the file. groupBy does not create a header (header=F).
    data <- read.table(file_path, sep = "\t", header = FALSE, col.names = c("Distance", "Count"))
    
    # Filter data to show only distances from 0 to 100
    data <- data[data$Distance >= 0 & data$Distance <= 100, ]
    
    # Check if data remains after filtering
    if (nrow(data) == 0) {
      plot(1, type="n", axes=F, xlab="", ylab="", main = title)
      text(1, 1, "No data in the 0-100 nt range", cex = 1.2)
      return() # Exit the function
    }
    
    # Order data by 'Distance' column for correct line plotting
    data <- data[order(data$Distance), ]
    
    # Calculate Z-score for the counts. Add a small value to avoid sd=0.
    z_score <- scale(data$Count + 1e-9)
    
    # Define Y-axis limits to ensure the reference line is visible
    ylim_range <- range(z_score, na.rm = TRUE)
    if (show_ref_line) {
      ylim_range <- range(c(ylim_range, 2.1), na.rm = TRUE)
    }
    
    # Create the line plot
    plot(data$Distance, z_score,
         type = "l", 
         lwd = 1,
         col = color,
         main = title,
         xlab = "Distance (nt)",
         ylab = "Frequency (Z-score)",
         ylim = ylim_range,
         # Set limits and disable default X-axis
         xlim = c(0, 100),
         xaxt = "n",
         # Add background grid
         panel.first = grid(lty = "dotted", col = "lightgray"))
    
    # Add the custom X-axis
    axis(1, at = seq(1, 81, by = 10))
    
    # Add the red reference line if requested
    if (show_ref_line) {
      abline(h = 2, lty = 2, col = "red", lwd = 1.5)
    }
    
  } else {
    # If the file doesn't exist or is empty, draw a blank plot with a message
    plot(1, type="n", axes=F, xlab="", ylab="", main = title)
    text(1, 1, "Data not found", cex = 1.2)
  }
}

# --- Create the three plots ---

# 1. Phasing Plot - Sense Strand (+)
create_line_plot(
  file_path = phasing_plus_file,
  title = "Frequency sense strand",
  color = "black",
  show_ref_line = TRUE
)

# 2. Phasing Plot - Antisense Strand (-)
create_line_plot(
  file_path = phasing_minus_file,
  title = "Frequency antisense strand",
  color = "black",
  show_ref_line = TRUE
)

# 3. Offset Plot
create_line_plot(
  file_path = offset_file,
  title = "Frequency offset (Z-score)",
  color = "black",
  show_ref_line = TRUE
)

# --- Finalize and close the PDF file ---
dev.off()

cat("--- R SCRIPT FINISHED ---\n")

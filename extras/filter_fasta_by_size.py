#!/usr/bin/env python3
"""
Input: fasta, minimal length, maximal length
Output: fasta
Return sequences whose lengths are within the range.
"""

import sys

def stop_err(msg):
    """Writes an error message to stderr and exits."""
    sys.stderr.write(f"ERROR: {msg}\n")
    sys.exit(1)

def main():
    """Main function to filter FASTA sequences by length."""
    if len(sys.argv) != 5:
        stop_err(f"Usage: {sys.argv[0]} <input.fasta> <min_len> <max_len> <output.fasta>")

    input_filename = sys.argv[1]
    output_filename = sys.argv[4]

    try:
        min_length = int(sys.argv[2])
    except ValueError:
        stop_err("Minimal length of the return sequence requires a numerical value.")
    
    try:
        max_length = int(sys.argv[3])
    except ValueError:
        stop_err("Maximum length of the return sequence requires a numerical value.")

    print(f"Filtering '{input_filename}' for sequences between {min_length} and {max_length} bp...")

    try:
        # Use 'with open' for safe file handling in Python 3
        with open(input_filename, 'r') as input_handle, open(output_filename, 'w') as output_handle:
            sequence_buffer = []
            header = None
            at_least_one_written = False

            for line in input_handle:
                line = line.strip()
                if not line:
                    continue
                
                if line.startswith('>'):
                    # If we have a stored sequence, process it before starting the new one
                    if header and sequence_buffer:
                        seq_len = len("".join(sequence_buffer))
                        if min_length <= seq_len <= max_length:
                            output_handle.write(f"{header}\n")
                            output_handle.write("\n".join(sequence_buffer) + "\n")
                            at_least_one_written = True
                    
                    # Start new sequence
                    header = line
                    sequence_buffer = []
                else:
                    # Add sequence line to buffer
                    sequence_buffer.append(line)
            
            # Process the very last sequence in the file
            if header and sequence_buffer:
                seq_len = len("".join(sequence_buffer))
                if min_length <= seq_len <= max_length:
                    output_handle.write(f"{header}\n")
                    output_handle.write("\n".join(sequence_buffer) + "\n")
                    at_least_one_written = True

            if not at_least_one_written:
                # Use modern print() function for Python 3
                print("Warning: There is no sequence that falls within your range.")

    except FileNotFoundError:
        stop_err(f"Input file not found at '{input_filename}'")
    except Exception as e:
        stop_err(f"An unexpected error occurred: {e}")

    print(f"Filtered sequences saved to '{output_filename}'.")

if __name__ == "__main__":
    main()

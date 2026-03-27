import argparse
from Bio import SeqIO

def fastq_to_fasta(input_file, output_file):
    with open(output_file, "w") as fasta_out:
        for record in SeqIO.parse(input_file, "fastq"):
            SeqIO.write(record, fasta_out, "fasta")

def main():
    parser = argparse.ArgumentParser(description="Converter arquivo FASTQ para FASTA")
    parser.add_argument("-i", "--input", required=True, help="Arquivo FASTQ de entrada")
    parser.add_argument("-o", "--output", required=True, help="Arquivo FASTA de saída")
    
    args = parser.parse_args()
    fastq_to_fasta(args.input, args.output)
    print(f"Conversão concluída: {args.output}")

if __name__ == "__main__":
    main()

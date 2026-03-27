import subprocess
import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import shutil

def checar_ferramentas():
    """Verifica se samtools e bedtools estão no PATH."""
    print("Verificando a presença das ferramentas necessárias...")
    ferramentas = ['samtools', 'bedtools']
    
    for ferramenta in ferramentas:
        if shutil.which(ferramenta) is None:
            print(f"ERRO: A ferramenta '{ferramenta}' não foi encontrada no seu PATH.")
            print("Verifique se ela está instalada e se o ambiente (ex: Conda) está ativado corretamente.")
            exit(1)
    print("Todas as ferramentas foram encontradas com sucesso. 👍")


def obter_cromossomos_do_sam(sam_file):
    """Extrai os nomes dos cromossomos e seus tamanhos do cabeçalho do SAM."""
    print(f"Extraindo nomes de cromossomos do arquivo {sam_file}...")
    cmd = f"samtools view -H {sam_file} | grep '^@SQ' | cut -f 2,3 | sed 's/SN://' | sed 's/LN://'"
    resultado = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
    if not resultado.stdout:
        print(f"ERRO: Não foi possível encontrar informações de cromossomos (@SQ) no cabeçalho de {sam_file}.")
        exit(1)
    genome_data = resultado.stdout.strip().split('\n')
    genome_file_path = 'genome.tmp.genome'
    with open(genome_file_path, 'w') as f:
        for line in genome_data: f.write(line + '\n')
    cromossomos = [line.split('\t')[0] for line in genome_data]
    print(f"Cromossomos encontrados: {', '.join(cromossomos)}")
    return cromossomos, genome_file_path

def processar_cromossomo(sorted_bam, genome_file, cromossomo):
    """Executa bedtools, mescla os resultados e gera o gráfico em PDF."""
    print(f"\nProcessando cromossomo: {cromossomo}...")
    
    # Define os nomes dos arquivos de saída
    saida_senso = f"cobertura_senso_{cromossomo}.tsv"
    saida_antisenso = f"cobertura_antisenso_{cromossomo}.tsv"
    saida_final_tabela = f"cobertura_ambas_fitas_{cromossomo}.tsv"
    # ALTERAÇÃO AQUI: A saída do gráfico agora será .pdf
    saida_final_grafico = f"cobertura_{cromossomo}.pdf"

    try:
        # Calcular coberturas com samtools e bedtools
        print("  Calculando coberturas senso e antissenso...")
        subprocess.run(f"samtools view -bh {sorted_bam} {cromossomo} | bedtools genomecov -ibam stdin -g {genome_file} -strand + -d > {saida_senso}", shell=True, check=True)
        subprocess.run(f"samtools view -bh {sorted_bam} {cromossomo} | bedtools genomecov -ibam stdin -g {genome_file} -strand - -d > {saida_antisenso}", shell=True, check=True)

        # Carregar, mesclar e processar os dados com pandas
        print("  Mesclando tabelas de cobertura...")
        df_senso = pd.read_csv(saida_senso, sep='\t', header=None, names=['cromossomo', 'posicao', 'cobertura_senso'])
        df_antisenso = pd.read_csv(saida_antisenso, sep='\t', header=None, names=['cromossomo', 'posicao', 'cobertura_antisenso'])
        df_merged = pd.merge(df_senso, df_antisenso, on=['cromossomo', 'posicao'])
        df_merged[['cromossomo', 'posicao', 'cobertura_antisenso', 'cobertura_senso']].to_csv(saida_final_tabela, sep='\t', index=False, header=False)
        print(f"  Tabela final salva em: {saida_final_tabela}")

        # Preparar dados para o gráfico
        df_merged['cobertura_antisenso_neg'] = df_merged['cobertura_antisenso'] * -1

        # Gerar o gráfico
        print("  Gerando gráfico em PDF...")
        plt.style.use('classic')
        fig, ax = plt.subplots(figsize=(15, 7))
        
        ax.fill_between(df_merged['posicao'], df_merged['cobertura_senso'], color="#996631", label="Senso (+)")
        ax.fill_between(df_merged['posicao'], df_merged['cobertura_antisenso_neg'], color="#0099CC", label="Antisenso (-)")

        max_abs_cov = max(df_merged['cobertura_senso'].max(), df_merged['cobertura_antisenso'].max()) * 1.1
        ax.set_ylim(-max_abs_cov, max_abs_cov)
        ax.set_xlabel(f"Posição no Cromossomo {cromossomo}", fontsize=14)
        ax.set_ylabel("Densidade de RNA (Cobertura)", fontsize=14)
        ax.set_title(f"Perfil de Cobertura de sRNA - {cromossomo}", fontsize=16, fontweight='bold')
        ax.axhline(0, color='black', linewidth=0.8)
        ax.legend()
        ax.grid(True, linestyle='--', alpha=0.6)
        
        plt.tight_layout()
        # ALTERAÇÃO AQUI: Salva diretamente em PDF
        plt.savefig(saida_final_grafico)
        plt.close()
        print(f"  Gráfico salvo em: {saida_final_grafico}")

    finally:
        # Limpar arquivos intermediários
        if os.path.exists(saida_senso): os.remove(saida_senso)
        if os.path.exists(saida_antisenso): os.remove(saida_antisenso)

def main():
    parser = argparse.ArgumentParser(description="Gerador de gráficos de cobertura de genoma a partir de um arquivo SAM.")
    parser.add_argument("sam_file", help="Caminho para o arquivo SAM de entrada.")
    
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-c", "--cromossomos", nargs='+', help="Lista de nomes de cromossomos a serem analisados.")
    group.add_argument("-f", "--file", help="Caminho para um arquivo de texto com um nome de cromossomo por linha.")
    
    args = parser.parse_args()

    checar_ferramentas()

    sam_file = args.sam_file
    if not os.path.exists(sam_file):
        print(f"ERRO: Arquivo SAM '{sam_file}' não encontrado.")
        exit(1)

    print("Convertendo SAM para BAM ordenado e indexado...")
    sorted_bam_path = sam_file.replace('.sam', '.sorted.bam')
    subprocess.run(f"samtools view -bS {sam_file} | samtools sort -o {sorted_bam_path}", shell=True, check=True)
    subprocess.run(f"samtools index {sorted_bam_path}", shell=True, check=True)
    print(f"BAM ordenado e indexado criado em: {sorted_bam_path}")
    
    todos_cromossomos, genome_file = obter_cromossomos_do_sam(sam_file)
    
    cromossomos_para_analisar = []
    if args.cromossomos:
        cromossomos_para_analisar = args.cromossomos
    elif args.file:
        try:
            with open(args.file, 'r') as f:
                cromossomos_para_analisar = [line.strip() for line in f if line.strip()]
            print(f"Lendo lista de cromossomos do arquivo: {args.file}")
        except FileNotFoundError:
            print(f"ERRO: O arquivo de lista '{args.file}' não foi encontrado.")
            exit(1)
    else:
        print("Nenhuma lista de cromossomos fornecida. Analisando todos os cromossomos encontrados no SAM.")
        cromossomos_para_analisar = todos_cromossomos
    
    for chrom in cromossomos_para_analisar:
        if chrom not in todos_cromossomos:
            print(f"AVISO: Cromossomo '{chrom}' fornecido não foi encontrado no cabeçalho do SAM. Ignorando.")
            continue
        processar_cromossomo(sorted_bam_path, genome_file, chrom)
        
    os.remove(genome_file)
    os.remove(sorted_bam_path)
    os.remove(sorted_bam_path + ".bai")
    print("\nProcesso concluído. Arquivos temporários foram removidos.")

if __name__ == "__main__":
    main()
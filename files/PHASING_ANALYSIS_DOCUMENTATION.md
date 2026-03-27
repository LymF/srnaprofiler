# 📊 Análise de Phasing e Offset - Documentação Completa

## 🎯 Visão Geral

O sistema de análise de **phasing** e **offset** detecta padrões de processamento de piRNAs, especificamente:

1. **Phasing**: Espaçamento periódico entre reads (tipicamente a cada ~10nt em piRNAs)
2. **Offset**: Sobreposição de 10nt na extremidade 5' entre reads sense e antisense (característica do "ping-pong cycle")

## 📚 Três Scripts Principais

### 1. phasing-script-teste-v4-debbug.pl (Pipeline Principal)
**Função**: Orquestra todo o fluxo de análise

### 2. bed_calculate_offsets.pl (Cálculo Core)
**Função**: Calcula offsets e phasing a partir de arquivos BED

### 3. plotOffsetPhasingSeries-teste.R (Visualização)
**Função**: Gera gráficos dos resultados

---

## 🔬 Fluxo Detalhado do Pipeline

### Pipeline Completo (phasing-script-teste-v4-debbug.pl)

```
Input: SAM file
  ↓
Step 1: Split por strand
  ├→ plus.sam  (reads forward)
  └→ minus.sam (reads reverse)
  ↓
Step 2: Filtrar por tamanho
  ├→ plus.21.sam  (apenas reads de 21nt)
  ├→ minus.21.sam
  └→ total.21.sam (ambos strands)
  ↓
Step 3: Converter para BED
  ├→ plus.21.bed
  ├→ minus.21.bed
  └→ total.21.bed
  ↓
Step 4: Calcular Phasing e Offset
  ├→ plus.21.phasing   (phasing do strand +)
  ├→ minus.21.phasing  (phasing do strand -)
  └→ total.21.offset   (overlap sense/antisense)
  ↓
Step 5: Plotar com R
  └→ prefix.21.phasing_offset_line.pdf
```

### Dois Modos de Operação

#### Modo 1: INDIVIDUAL (default)
```bash
perl phasing.pl -s reads.sam -p sample -f results/ \
    -c chr1 -r 1-10000 -si 21 -se 24

# Gera análises separadas:
#   sample.21.* (todos os arquivos para 21nt)
#   sample.22.* (todos os arquivos para 22nt)
#   sample.23.* (todos os arquivos para 23nt)
#   sample.24.* (todos os arquivos para 24nt)
```

#### Modo 2: COMBINED
```bash
perl phasing.pl -s reads.sam -p sample -f results/ \
    -c chr1 -r 1-10000 -si 21 -se 24 --combined

# Gera uma única análise:
#   sample.21-24.* (pooled de todos 21-24nt)
```

---

## 🧮 Algoritmo de Cálculo de Offset (bed_calculate_offsets.pl)

### Conceito Fundamental

O script calcula, para cada read, quantos outros reads começam em posições específicas ao redor dele.

### Input: Arquivo BED

```
chr1	95	125	read1-7	  0	+
chr1	100	150	read2-13  0	+
chr1	100	140	read3-4	  0	+
chr1	80	111	read4-20  0	-
```

Colunas:
1. Cromossomo
2. Start (0-based)
3. End (1-based)
4. Nome (pode ter multiplicidade: "N-M" onde M = número de cópias)
5. Score (ignorado)
6. Strand

### Parâmetros Principais

```perl
--region chr1:start-end  # Região de interesse
--up 80                  # Janela upstream (80nt antes)
--down 80                # Janela downstream (80nt depois)
--mult                   # Usar multiplicidade do nome
--header                 # Print header
```

### Algoritmo Passo a Passo

#### Passo 1: Inicialização

```perl
# Criar arrays vazios para toda a região
@window_positive_sequences[region_size]
@window_negative_sequences[region_size]
@window_positive_reads[region_size]
@window_negative_reads[region_size]
```

#### Passo 2: Calcular Cobertura

Para cada read na região:

```perl
if (strand == "+") {
    # Para forward: posição 5' = start
    window_positive_sequences[start - region_start] += 1
    window_positive_reads[start - region_start] += multiplicity
}
else {
    # Para reverse: posição 5' = end
    window_negative_sequences[end - region_start] += 1
    window_negative_reads[end - region_start] += multiplicity
}
```

**Resultado**: Arrays com contagem de reads por posição

```
Posição:    95  96  97  98  99 100 101 ...
Pos (+):     2   0   0   0   0   2   0  ...
Neg (-):     0   0   0   0   0   0   0  ...
```

#### Passo 3: Scan de Janelas

Para cada read processado:

1. **Determinar posição 5'**:
   - Forward (+): `nuc_position = start`
   - Reverse (-): `nuc_position = end`

2. **Definir janela**:
   ```perl
   if (strand == "+") {
       window_start = nuc_position - upstream    # 80nt antes
       window_end = nuc_position + downstream    # 80nt depois
   }
   else {
       # Inverter para reverse
       window_start = nuc_position - downstream
       window_end = nuc_position + upstream
   }
   ```

3. **Para cada posição na janela**:
   ```perl
   for (window_index in window_start..window_end) {
       # Calcular offset
       if (strand == "+") {
           offset = window_index - nuc_position
       }
       else {
           offset = nuc_position - window_index
       }
       
       # Contar reads na mesma strand
       count_same = window_sequences[strand][window_index]
       
       # Contar reads na strand oposta
       count_opposite = window_sequences[opposite_strand][window_index]
       
       # Output se houver dados
       if (offset == 0) {
           # No offset 0, subtrair o próprio read
           count_same -= 1
       }
       
       print: chrom, genomic_pos, strand, offset, count_same, count_opposite
   }
   ```

### Output: Arquivo de Offsets

```
chrom	genomic_start	strand	window_pos	offset	seq_same	read_same	seq_opposite	read_opposite
chr1	95	+	95	5	2	17	0	0
chr1	100	+	100	-5	1	7	0	0
chr1	100	+	100	0	1	13	0	0
chr1	100	+	100	10	0	0	1	20
```

**Interpretação**:
- `offset = 0`: Posição do próprio read
- `offset = 10`: 10nt downstream (para +) ou upstream (para -)
- `offset = -10`: 10nt upstream (para +) ou downstream (para -)
- `seq_opposite`: Número de reads na strand oposta nessa posição
- `read_opposite`: Total de reads considerando multiplicidade

---

## 📈 Processamento para Plots

### Arquivo de Phasing (plus.phasing / minus.phasing)

Gerado com:
```bash
bed_calculate_offsets.pl ... | \
    grep -v "intervals" | \
    grep -v "chrom" | \
    sed 1d | \
    sort -k5,5 | \
    groupBy -g 5 -c 7 -o sum
```

**O que faz**:
1. Remove linhas "intervals" e "chrom"
2. Remove primeira linha
3. Ordena por offset (coluna 5)
4. Agrupa por offset e soma contagens (coluna 7 = read_count_same_strand)

**Output**:
```
-80	15
-79	12
-78	18
...
0	250
1	35
...
10	189  ← Pico de phasing!
```

### Arquivo de Offset (total.offset)

Gerado com:
```bash
bed_calculate_offsets.pl ... | \
    grep -v "intervals" | \
    grep -v "chrom" | \
    sed 1d | \
    sort -k5,5 | \
    groupBy -g 5 -c 9 -o sum
```

**Diferença**: Usa coluna 9 = `read_count_opposite_strand`

**Output**:
```
-80	5
-79	3
...
10	420  ← Pico de overlap sense/antisense!
```

**Interpretação Biológica**:
- **Pico em offset +10**: "Ping-pong signature"
- Indica que reads antisense começam 10nt downstream dos reads sense
- Característico do ciclo de amplificação de piRNAs

---

## 🎨 Script R de Visualização

### Input

O script R recebe 3 arquivos:

```R
args[1] = plus_phasing   # Phasing do strand +
args[2] = minus_phasing  # Phasing do strand -
args[3] = offset_file    # Overlap sense/antisense
args[4] = output_prefix  # Para nomear PDF
```

**Nota**: O script atual só usa `args[1]` (offset_file) e `args[2]` (output_prefix)

### Processamento

```R
# 1. Ler dados
off <- read.csv(offset_file, sep="\t", header=TRUE)

# Formato esperado:
#   Overlap    Count
#   -80        15
#   -79        12
#   ...

# 2. Escalar valores (normalização Z-score)
v_off <- scale(off$Count)
# v_off = (Count - mean(Count)) / sd(Count)

# 3. Plotar
plot(off$Overlap, v_off,
     type="l",                           # Linha
     xaxp=c(-30,30,30),                 # 30 ticks de -30 a +30
     ylim=c(y_min, y_max),
     xlab="Distance from 5' end (nt)",
     ylab="Relative freq. of antisense 5' end")

# 4. Adicionar linhas de referência
abline(v=seq(-29,29,by=2), lty=2, col="lightgray")  # Grid
abline(v=c(18), lty=2, col="blue")                   # Linha azul em +18
abline(h=c(2), lty=2, col="red")                     # Linha vermelha em y=2
```

### Output

PDF com:
- **Eixo X**: Offset (-30 a +30 nt)
- **Eixo Y**: Frequência relativa escalada
- **Pico esperado**: Em torno de +10 (ping-pong signature)

---

## 🐍 Implementação em Python - Estratégia

### Módulo Proposto: `phasing.py`

```python
class PhasingAnalysis:
    """Análise de phasing e offset para piRNAs"""
    
    def __init__(self, sam_file, region_chrom, region_start, region_end):
        self.sam_file = sam_file
        self.region_chrom = region_chrom
        self.region_start = region_start
        self.region_end = region_end
        self.region_size = region_end - region_start + 1
        
        # Arrays de cobertura
        self.window_positive_sequences = [0] * self.region_size
        self.window_negative_sequences = [0] * self.region_size
        self.window_positive_reads = [0] * self.region_size
        self.window_negative_reads = [0] * self.region_size
        
        self.intervals = []
    
    def calculate_coverage_from_sam(self, size_min=None, size_max=None):
        """
        Calcula cobertura diretamente do SAM.
        Equivalente a: SAM → BED → calculate_offsets
        """
        with open(self.sam_file) as f:
            for line in f:
                if line.startswith('@'):
                    continue
                
                parts = line.strip().split('\t')
                chrom = parts[2]
                start = int(parts[3]) - 1  # SAM é 1-based
                flag = int(parts[1])
                seq = parts[9]
                
                # Filtrar região
                if chrom != self.region_chrom:
                    continue
                
                # Filtrar tamanho
                read_len = len(seq)
                if size_min and read_len < size_min:
                    continue
                if size_max and read_len > size_max:
                    continue
                
                # Calcular end
                end = start + read_len
                
                # Filtrar região
                if start < self.region_start or end > self.region_end:
                    continue
                
                # Determinar strand
                strand = '-' if (flag & 0x10) else '+'
                
                # Multiplicidade do nome (se existir)
                read_name = parts[0]
                multiplicity = self._extract_multiplicity(read_name)
                
                # Posição 5'
                if strand == '+':
                    pos_5prime = start
                else:
                    pos_5prime = end - 1
                
                # Atualizar cobertura
                index = pos_5prime - self.region_start
                if 0 <= index < self.region_size:
                    if strand == '+':
                        self.window_positive_sequences[index] += 1
                        self.window_positive_reads[index] += multiplicity
                    else:
                        self.window_negative_sequences[index] += 1
                        self.window_negative_reads[index] += multiplicity
                
                # Guardar intervalo
                self.intervals.append({
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'multiplicity': multiplicity
                })
    
    def calculate_offsets(self, upstream=80, downstream=80):
        """
        Calcula offsets para cada intervalo.
        Retorna lista de dicionários com offset data.
        """
        results = []
        
        for interval in self.intervals:
            chrom = interval['chrom']
            start = interval['start']
            end = interval['end']
            strand = interval['strand']
            multiplicity = interval['multiplicity']
            
            # Posição 5'
            if strand == '+':
                nuc_position = start - self.region_start
            else:
                nuc_position = end - 1 - self.region_start
            
            # Janela
            if strand == '+':
                window_start = max(0, nuc_position - upstream)
                window_end = min(self.region_size - 1, nuc_position + downstream)
            else:
                window_start = max(0, nuc_position - downstream)
                window_end = min(self.region_size - 1, nuc_position + upstream)
            
            # Scan janela
            for window_index in range(window_start, window_end + 1):
                # Calcular offset
                if strand == '+':
                    offset = window_index - nuc_position
                else:
                    offset = nuc_position - window_index
                
                # Contagens
                if strand == '+':
                    seq_same = self.window_positive_sequences[window_index]
                    read_same = self.window_positive_reads[window_index]
                    seq_opposite = self.window_negative_sequences[window_index]
                    read_opposite = self.window_negative_reads[window_index]
                else:
                    seq_same = self.window_negative_sequences[window_index]
                    read_same = self.window_negative_reads[window_index]
                    seq_opposite = self.window_positive_sequences[window_index]
                    read_opposite = self.window_positive_reads[window_index]
                
                # No offset 0, subtrair próprio read
                if offset == 0:
                    seq_same -= 1
                    read_same -= multiplicity
                
                # Apenas registrar se houver dados
                if seq_same > 0 or read_same > 0 or seq_opposite > 0 or read_opposite > 0:
                    results.append({
                        'chrom': chrom,
                        'genomic_pos': self.region_start + nuc_position,
                        'strand': strand,
                        'offset': offset,
                        'seq_same': seq_same,
                        'read_same': read_same,
                        'seq_opposite': seq_opposite,
                        'read_opposite': read_opposite
                    })
        
        return results
    
    def aggregate_phasing(self, offset_data, by_strand='+'):
        """
        Agrega dados de phasing por offset.
        Equivalente a: groupBy -g 5 -c 7 -o sum
        """
        from collections import defaultdict
        
        aggregated = defaultdict(int)
        
        for row in offset_data:
            if row['strand'] == by_strand:
                offset = row['offset']
                count = row['read_same']
                aggregated[offset] += count
        
        # Converter para lista ordenada
        result = []
        for offset in sorted(aggregated.keys()):
            result.append({'offset': offset, 'count': aggregated[offset]})
        
        return result
    
    def aggregate_offset(self, offset_data):
        """
        Agrega overlap sense/antisense por offset.
        Equivalente a: groupBy -g 5 -c 9 -o sum
        """
        from collections import defaultdict
        
        aggregated = defaultdict(int)
        
        for row in offset_data:
            offset = row['offset']
            count = row['read_opposite']
            aggregated[offset] += count
        
        # Converter para lista ordenada
        result = []
        for offset in sorted(aggregated.keys()):
            result.append({'offset': offset, 'count': aggregated[offset]})
        
        return result
    
    def write_phasing_file(self, phasing_data, output_file):
        """Escreve arquivo de phasing."""
        with open(output_file, 'w') as f:
            for row in phasing_data:
                f.write(f"{row['offset']}\t{row['count']}\n")
    
    def write_offset_file(self, offset_data, output_file):
        """Escreve arquivo de offset."""
        with open(output_file, 'w') as f:
            f.write("Overlap\tCount\n")
            for row in offset_data:
                f.write(f"{row['offset']}\t{row['count']}\n")
    
    def plot_with_r(self, offset_file, output_prefix, r_script_path):
        """Chama script R para gerar plot."""
        import subprocess
        
        cmd = [
            "Rscript",
            r_script_path,
            str(offset_file),
            str(output_prefix)
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        return result.returncode == 0
    
    @staticmethod
    def _extract_multiplicity(read_name):
        """Extrai multiplicidade do nome: 'N-M' → M"""
        if '-' in read_name:
            try:
                parts = read_name.split('-')
                return int(parts[-1])
            except:
                pass
        return 1
```

### Uso do Módulo

```python
from phasing import PhasingAnalysis

# 1. Inicializar
analysis = PhasingAnalysis(
    sam_file="sample.sam",
    region_chrom="chr1",
    region_start=1000,
    region_end=10000
)

# 2. Calcular cobertura (com filtro de tamanho)
analysis.calculate_coverage_from_sam(size_min=24, size_max=30)

# 3. Calcular offsets
offset_data = analysis.calculate_offsets(upstream=80, downstream=80)

# 4. Agregar para plots
phasing_plus = analysis.aggregate_phasing(offset_data, by_strand='+')
phasing_minus = analysis.aggregate_phasing(offset_data, by_strand='-')
offset_agg = analysis.aggregate_offset(offset_data)

# 5. Salvar arquivos
analysis.write_phasing_file(phasing_plus, "output.plus.phasing")
analysis.write_phasing_file(phasing_minus, "output.minus.phasing")
analysis.write_offset_file(offset_agg, "output.offset")

# 6. Gerar plot com R
analysis.plot_with_r(
    offset_file="output.offset",
    output_prefix="output",
    r_script_path="plotOffsetPhasingSeries.R"
)
```

---

## ✅ Checklist de Implementação

- [ ] Criar classe `PhasingAnalysis`
- [ ] Implementar `calculate_coverage_from_sam()`
- [ ] Implementar `calculate_offsets()`
- [ ] Implementar agregação (phasing e offset)
- [ ] Implementar escrita de arquivos
- [ ] Integrar com módulo plots.py
- [ ] Adicionar parâmetros CLI
- [ ] Testes com dados reais
- [ ] Validar contra script Perl original

---

## 📊 Arquivos de Output Esperados

Para análise de 21-24nt:

```
results/
├── sample.21.plus.phasing      # Phasing + para 21nt
├── sample.21.minus.phasing     # Phasing - para 21nt  
├── sample.21.offset            # Overlap para 21nt
├── sample.21.pdf               # Plot para 21nt
├── sample.22.plus.phasing
├── sample.22.minus.phasing
├── sample.22.offset
├── sample.22.pdf
... (e assim por diante para 23 e 24)
```

Ou modo combined:
```
results/
├── sample.21-24.plus.phasing
├── sample.21-24.minus.phasing
├── sample.21-24.offset
└── sample.21-24.pdf
```

---

## 🎯 Próximos Passos

1. **Implementar classe Python** com lógica fiel ao Perl
2. **Integrar no plots.py** como nova função
3. **Adicionar CLI** para phasing standalone
4. **Manter script R original** para visualização
5. **Testes extensivos** comparando com pipeline Perl

**Vantagem**: Todo processamento em Python, apenas visualização em R!

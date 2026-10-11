# Genômica comparativa e ontogenia em Archosauria

**_Gallus gallus_ frente a _Alligator mississippiensis_ e _Crocodylus porosus_**

Pipeline reprodutível que contrasta o genoma da galinha com três genomas de crocodilianos, usando apenas dados públicos do NCBI. Quantifica a aceleração da taxa de substituição na linhagem aviária, a compartimentação GC associada a microcromossomos e os elementos não-codificantes conservados nos quatro clusters Hox.

> *A reproducible comparative-genomics pipeline contrasting the chicken genome with three crocodilian genomes. Report and code comments are in Portuguese; all data come from public NCBI records.*

---

## Resultados principais

| Achado | Valor | Onde |
|---|---|---|
| Aceleração aviária em sítios quase neutros (3ª posição de códon) | **1,90×** | `results/tabela_taxas_lineagem.csv` |
| Aceleração em sítios sob seleção purificadora (proteína / códon 1+2) | 1,28× / 1,26× | idem |
| Teste de taxa relativa de Tajima (*Chrysemys* como grupo externo) | χ² = 460,8; *P* = 3,2 × 10⁻¹⁰² | `results/tabela_tajima_rrt.csv` |
| Microcromossomos (<20 Mb): galinha vs jacaré | **28 vs 0** | `results/tabela2_cariotipo_gc.csv` |
| Correlação tamanho–GC preservada mesmo sem microcromossomos | ρ = −0,96 (*Gallus*), −0,92 (*Alligator*) | idem |
| CNEs galinha × jacaré nos clusters Hox | **251 elementos, 93.031 bp, 79,8% de identidade** | `results/tabela6_cnes_hox.csv` |
| Compactação do genoma aviário | 1,053 Gb vs 2,347 Gb | `results/tabela1_estatisticas_genomicas.csv` |

A leitura de fundo: a aceleração da linhagem aviária concentra-se nos sítios quase neutros e não nos sítios sob restrição funcional, o que é a assinatura esperada de uma redução do tempo de geração — e não de seleção positiva generalizada. No jacaré, a fusão dos microcromossomos apagou-os como entidades citogenéticas sem apagar o compartimento genômico GC-rico que eles carregavam.

![Arquitetura genômica](figures/fig1_arquitetura_genomica.png)

![Taxas moleculares](figures/fig2_taxas_moleculares.png)

![CNEs nos clusters Hox](figures/fig3_hox_cnes.png)

O relatório completo, com as três fases da análise (genômica comparativa, filogenia/tempo profundo e evo-devo) e 57 referências com DOI verificado, está em **[`docs/relatorio_archosauria.md`](docs/relatorio_archosauria.md)**.

---

## Como reproduzir

```bash
conda env create -f environment.yml
conda activate archosaur
export NCBI_CONTACT_EMAIL="voce@exemplo.org"   # o NCBI pede um contato
./run_all.sh
```

Os passos 1, 3 e 7 consultam as APIs do NCBI e levam alguns minutos. Os produtos intermediários já versionados (`data/`, `alignments/`) permitem executar qualquer passo isoladamente sem refazer as buscas.

**Nota de ambiente:** `src/05_phylogeny.sh` exporta `KMP_AFFINITY=disabled`. Sem isso, o runtime OpenMP do IQ-TREE aborta com `pthread_setaffinity_np() failed / Operation not permitted` em contêineres que bloqueiam chamadas de afinidade de CPU.

---

## Estrutura

```
src/          pipeline em 10 passos (ncbi.py é o cliente das APIs)
data/         respostas brutas do NCBI, ortólogos e sequências por gene
alignments/   supermatrizes (proteína, códon, pos. 1+2, pos. 3) e árvores IQ-TREE
results/      tabelas 1-6 em CSV
figures/      figuras 1-3 em PNG 300 dpi
docs/         relatório científico completo
```

| Passo | Script | Produz |
|---|---|---|
| 1 | `01_fetch_assemblies.py` | Tabela 1, relatórios de sequência |
| 2 | `02_karyotype_gc.py` | Tabela 2 |
| 3 | `03_fetch_orthologs.py` | 39 ortólogos 1:1 em 7 táxons |
| 4 | `04_build_supermatrix.py` | supermatrizes (24.375 aa / 73.125 nt) |
| 5 | `05_phylogeny.sh` | árvore ML + ramos em topologia fixa |
| 6 | `06_rates_tajima.py` | Tabelas 3 e 4 |
| 7 | `07_fetch_hox.py` | Tabela 5, regiões Hox |
| 8 | `08_cne_scan.py` | Tabela 6 |
| 9 | `09_figures.py` | Figuras 1-3 |
| 10 | `10_bibliography.py` | resolução de DOIs via CrossRef |

---

## Métodos, em resumo

Ortologia 1:1 obtida do NCBI Datasets para sete amniotas; proteína e CDS pareados por **tradução exata**, de modo que a retrotradução para códons é inequívoca. Alinhamento com MAFFT `--auto`, colunas retidas com ≥6/7 táxons presentes. IQ-TREE estima primeiro a árvore ML sem restrição (que recupera a topologia aceita: Crocodylia UFBoot 100, Longirostres 74) e em seguida os comprimentos de ramo sobre a topologia fixa, em três partições. As distâncias são medidas **a partir do ancestral comum dos arcossauros**, de modo que o tempo decorrido é idêntico para as duas linhagens comparadas. O teste de Tajima (1993) é aplicado à supermatriz de códons com três enraizamentos. Os CNEs vêm de `dc-megablast` entre regiões Hox ortólogas, com os éxons codificantes subtraídos por alinhamento dos CDS Hox contra a própria região.

**Limitações.** As contagens gênicas de *C. porosus* e *G. gangeticus* provêm de montagens de scaffold de 2016 e refletem fragmentação, não conteúdo gênico real. As taxas derivam de 39 genes conservados, o que torna 1,90× um limite **inferior**. Não houve datação molecular. A detecção de CNE por BLAST entre pares é um limite inferior do conjunto real.

---

## Dados e licença

Montagens e anotações: NCBI RefSeq (`GCF_016699485.2`, `GCF_030867095.1`, `GCF_001723895.1`, `GCF_001723915.1`, `GCF_011386835.1`, `GCF_035594765.1`, `GCF_048771995.1`). APIs: NCBI Datasets v2, NCBI E-utilities, UCSC Genome Browser, PubMed/PMC, CrossRef.

Código sob licença MIT; relatório e figuras sob CC BY 4.0. Textos completos de artigos não são redistribuídos aqui.

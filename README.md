# pipelinegenomedemontrationandsimulation
This a  simple repo only for experiments and reviews. You can use to consult and enjoy experiments with GEO DATABases


Passo a Passo da Análise de RNA-seq (Transcriptoma)

Para este experimento, realizaremos uma análise de dados de sequenciamento de RNA (RNA-seq) para identificação de genes diferencialmente expressos e módulos de coexpressão. As seguintes ferramentas foram utilizadas no processo: R, uma linguagem de programação estatística amplamente utilizada em bioinformática, e os pacotes especializados GEOquery, limma e WGCNA.

Para simplificar o processo e assimilar o conteúdo, a base de dados foi limitada a 300 amostras. Isso também serve como um ponto de partida para que outros pesquisadores possam reproduzir e aprimorar o estudo utilizando bases de dados completas disponíveis no GEO (Gene Expression Omnibus).

Para reproduzir este estudo, utilize os seguintes processos:
Passo 1: Instale e Carregue os Pacotes Essenciais

Abra uma nova sessão do RStudio e execute o seguinte código para garantir que todos os pacotes necessários estejam instalados e carregados.
R

# Instale os pacotes (apenas na primeira vez)
install.packages("BiocManager")
BiocManager::install(c("GEOquery", "limma", "WGCNA"))

# Carregue os pacotes na sua sessão de trabalho
library(GEOquery)
library(limma)
library(WGCNA)

Passo 2: Baixe e Prepare os Dados do GEO

Utilize o pacote GEOquery para baixar o conjunto de dados de câncer de mama GSE42568 diretamente do repositório. Em seguida, selecione as 300 primeiras amostras para a análise.
R

# Baixe os dados do GEO
gse <- getGEO("GSE42568", GSEMatrix = TRUE)
eset <- gse[[1]] # Extrai o conjunto de dados

# Crie um subconjunto com as 300 primeiras amostras
eset_subset <- eset[, 1:300]

Passo 3: Análise de Expressão Diferencial com limma

Esta etapa identifica os genes que se expressam de forma significativamente diferente entre os grupos de amostras (tumor vs. controle).
R

# Extraia os dados de expressão e as informações das amostras
exprs_data <- exprs(eset_subset)
pheno_data <- pData(eset_subset)

# Identifique os grupos
groups <- factor(pheno_data$characteristics_ch1.2) # O nome da coluna pode variar

# Crie a matriz de design
design <- model.matrix(~0 + groups)

# Execute a análise limma
fit <- lmFit(exprs_data, design)
contrast.matrix <- makeContrasts(tumor-normal, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Obtenha e salve os resultados
results_limma <- topTable(fit2, number = Inf, adjust.method = "BH")
write.csv(results_limma, "results/tabela_limma_cancer.csv")

Passo 4: Análise de Coexpressão com WGCNA

Nesta fase, módulos (grupos) de genes que trabalham em conjunto são identificados.
R

# Transponha a matriz de dados de expressão para o WGCNA
datExpr <- as.data.frame(t(exprs_data))

# Escolha o soft-thresholding power
sft <- pickSoftThreshold(datExpr, networkType = "unsigned")

# Construa a rede e identifique os módulos
net <- blockwiseModules(datExpr,
                        power = sft$powerEstimate,
                        networkType = "unsigned",
                        minModuleSize = 30,
                        maxBlockSize = 20000,
                        mergeCutHeight = 0.25,
                        verbose = 3)

# Salve os resultados dos módulos
module_colors <- net$colors
table(module_colors)

# Salve os dados de WGCNA
save(net, file = "results/wgcna_cancer_results.RData")



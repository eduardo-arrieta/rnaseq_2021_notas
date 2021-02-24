## Lets build our first SummarizedExperiment object
library("SummarizedExperiment")
## ?SummarizedExperiment

## De los ejemplos en la ayuda oficial

## Creamos los datos para nuestro objeto de tipo SummarizedExperiment
## para 200 genes a lo largo de 6 muestras
nrows <- 200
ncols <- 6
## Números al azar de cuentas
set.seed(20210223)
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
## Información de nuestros genes
rowRanges <- GRanges(
  rep(c("chr1", "chr2"), c(50, 150)),
  IRanges(floor(runif(200, 1e5, 1e6)), width = 100),
  strand = sample(c("+", "-"), 200, TRUE),
  feature_id = sprintf("ID%03d", 1:200)
)
names(rowRanges) <- paste0("gene_", seq_len(length(rowRanges)))
## Información de nuestras muestras
colData <- DataFrame(
  Treatment = rep(c("ChIP", "Input"), 3),
  row.names = LETTERS[1:6]
)
## Juntamos ahora toda la información en un solo objeto de R
rse <- SummarizedExperiment(
  # Las tablas de ASSAY deben ser de las mismas dimensiones
  assays = SimpleList(counts = counts),
  rowRanges = rowRanges,
  colData = colData
)

## Exploremos el objeto resultante
rse

# Cantidad de cromosomas
seqlevels(rse)

## Ejercicio
#---------------------
## Comando 1
# Imprime los primeros dos genes del objeto a lo largo de todas la muestras
rse[1:2, ]
## Comando 2
# Imprimirá todos los genes en las muestras A, D Y F
rse[, c("A", "D", "F")]

#---------------------

# Almacenado en colData Names
rse$Treatment

## Explora el objeto rse de forma interactiva
library("iSEE")
iSEE::iSEE(rse)

# Ejercicio ISEE
## Descarguemos unos datos de spatialLIBD
sce_layer <- spatialLIBD::fetch_data("sce_layer")
iSEE::iSEE(sce_layer)

rowRanges(sce_layer[c('ENSG00000168314', 'ENSG00000183036', 'ENSG00000197971'),])$gene_name
# [1] "MOBP" "PCP4" "MBP"

# Generación de un subset
sce_short_3 <- sce_layer[c('ENSG00000168314', 'ENSG00000183036', 'ENSG00000197971'),]
iSEE::iSEE(sce_short_3)
sce_short_3
iSEE::iSEE(sce_layer[c('ENSG00000168314', 'ENSG00000183036', 'ENSG00000197971'),])

## Load recount3 R package
library("recount3")
## Revisemos todos los proyectos con datos de humano en recount3
human_projects <- available_projects()
## Encuentra tu proyecto de interés. Aquí usaremos
## SRP009615 de ejemplo
proj_info <- subset(
  human_projects,
  project == "SRP009615" & project_type == "data_sources"
)
## Crea un objetio de tipo RangedSummarizedExperiment (RSE)
## con la información a nivel de genes
rse_gene_SRP009615 <- create_rse(proj_info)

## Explora los proyectos disponibles de forma interactiva
proj_info_interactive <- interactiveDisplayBase::display(human_projects)
## Selecciona un solo renglón en la tabla y da click en "send".

## Aquí verificamos que solo seleccionaste un solo renglón.
stopifnot(nrow(proj_info_interactive) == 1)
## Crea el objeto RSE
rse_gene_interactive <- create_rse(proj_info_interactive)

## Convirtamos las cuentas por nucleotido a cuentas por lectura
## usando compute_read_counts().
## Para otras transformaciones como RPKM y TPM, revisa transform_counts().
assay(rse_gene_SRP009615, "counts") <- compute_read_counts(rse_gene_SRP009615)

## Para este estudio en específico, hagamos más fácil de usar la
## información del experimento
rse_gene_SRP009615 <- expand_sra_attributes(rse_gene_SRP009615)
colData(rse_gene_SRP009615)[
  ,
  grepl("^sra_attribute", colnames(colData(rse_gene_SRP009615)))
]
# DataFrame with 12 rows and 4 columns
# sra_attribute.cells sra_attribute.shRNA_expression sra_attribute.source_name sra_attribute.treatment
# <character>                    <character>               <character>             <character>
#   SRR387777                K562                             no                    SL2933               Puromycin
# SRR387778                K562             yes, targeting SRF                    SL2934  Puromycin, doxycycline
# SRR387779                K562                             no                    SL5265               Puromycin
# SRR387780                K562              yes targeting SRF                    SL3141  Puromycin, doxycycline
# SRR389079                K562            no shRNA expression                    SL6485               Puromycin
# ...                       ...                            ...                       ...                     ...
# SRR389082                K562         expressing shRNA tar..                    SL2592  Puromycin, doxycycline
# SRR389083                K562            no shRNA expression                    SL4337               Puromycin
# SRR389084                K562         expressing shRNA tar..                    SL4326  Puromycin, doxycycline
# SRR389077                K562            no shRNA expression                    SL1584               Puromycin
# SRR389078                K562         expressing shRNA tar..                    SL1583  Puromycin, doxycycline

# Ejercicio 3
# -----------------------
# Ploeo de expresión
iSEE::iSEE(rse_gene_SRP009615)
# La imagen se guardó en figuras
#------------------------

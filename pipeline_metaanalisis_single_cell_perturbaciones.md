# Pipeline para un metaanálisis computacional de perturbaciones single-cell y multi-ómicas

## Título propuesto

**Metaanálisis computacional single-cell de perturbaciones genéticas y farmacológicas para inferir lógica regulatoria, respuesta celular y reversión de estados patológicos**

---

## 1. Objetivo general

Construir un marco computacional reproducible que integre datasets públicos de single-cell RNA-seq, single-cell multiome, Perturb-seq, CRISPRi/a, CROP-seq, drug-seq y atlas celulares humanos para responder:

```text
Dado un estado celular basal + contexto biológico + perturbación genética/farmacológica,
¿se puede predecir la respuesta transcriptómica y regulatoria de forma transferible entre estudios?
```

El estudio no requiere generación experimental propia. La contribución principal es una combinación de:

```text
metaanálisis multi-estudio
+ armonización de datasets públicos
+ reconstrucción regulatoria
+ benchmark de modelos predictivos
+ priorización in silico de perturbaciones farmacológicas
```

---

## 2. Hipótesis central

Los datasets públicos de perturbación single-cell y multi-ómica contienen suficiente señal para identificar módulos reguladores y respuestas celulares conservadas entre contextos. Sin embargo, los modelos actuales pueden estar capturando principalmente estructura correlacional y coexpresión, no necesariamente lógica regulatoria causal.

Por tanto, el proyecto debe evaluar explícitamente:

```text
1. conservación de respuestas perturbacionales
2. dependencia del contexto celular
3. capacidad predictiva de modelos actuales
4. recuperación de relaciones regulatorias causales
5. utilidad para priorización farmacológica in silico
```

---

## 3. Diseño general del estudio

El proyecto debe plantearse como:

```text
Systematic computational meta-analysis
+ reproducible benchmark
+ cross-dataset regulatory reanalysis
+ perturbation-response prediction framework
```

No debe plantearse como una revisión narrativa. El valor científico surge del reanálisis cuantitativo y del benchmark reproducible.

---

## 4. Preguntas específicas

1. ¿Qué respuestas transcriptómicas son conservadas entre estudios, tipos celulares y plataformas?
2. ¿Qué respuestas son específicas de tipo celular, enfermedad, tejido, línea celular o batch experimental?
3. ¿Qué modelos predicen mejor perturbaciones no observadas?
4. ¿Qué redes reguladoras explican mejor los cambios transcriptómicos observados?
5. ¿Los modelos fundacionales single-cell capturan lógica regulatoria causal o sólo coexpresión?
6. ¿Qué perturbaciones farmacológicas revierten estados patológicos in silico?
7. ¿Qué nivel de incertidumbre acompaña a cada predicción?

---

## 5. Alcance recomendado

### 5.1 Alcance inicial recomendado

Para un primer proyecto, se recomienda iniciar con un dominio restringido:

```text
Dominio: cáncer / líneas celulares humanas
Modalidades iniciales: scRNA-seq + perturbaciones farmacológicas/genéticas
Extensión posterior: RNA+ATAC / multiome
```

### 5.2 Justificación

El dominio de cáncer y líneas celulares tiene ventajas prácticas:

```text
- abundancia de perturbaciones farmacológicas
- disponibilidad de líneas celulares bien caracterizadas
- existencia de datasets masivos de drug perturbation
- mayor conexión con descubrimiento de fármacos
- posibilidad de usar firmas de reversión patológica
```

### 5.3 Alcance no recomendado para la primera versión

No iniciar con:

```text
todas las células humanas
+ todas las enfermedades
+ todas las modalidades multi-ómicas
+ todas las perturbaciones
```

Ese alcance es demasiado amplio y dificulta producir una contribución controlada.

---

## 6. Fuentes de datos

### 6.1 Perturbaciones genéticas

Fuentes recomendadas:

```text
- Perturb-seq
- CRISPRi
- CRISPRa
- CROP-seq
- Replogle perturbation datasets
- Norman combinatorial perturbations
- Dixit perturbation datasets
- Adamson perturbation datasets
- scPerturb
```

### 6.2 Perturbaciones farmacológicas

Fuentes recomendadas:

```text
- Tahoe-100M
- sci-Plex
- drug-seq
- LINCS / CMap
- datasets compatibles con CPA / chemCPA
```

### 6.3 Multi-ómica regulatoria

Fuentes recomendadas:

```text
- single-cell multiome RNA+ATAC
- Perturb-ATAC
- Spear-ATAC
- CRISPRmap
- scATAC-seq atlases
- datasets multiome pareados RNA+ATAC
```

### 6.4 Atlas de referencia

Fuentes recomendadas:

```text
- Human Cell Atlas
- CZ CELLxGENE
- CELLxGENE Census
- Tabula Sapiens
- atlas específicos de enfermedad
- atlas específicos de tejido
```

---

## 7. Criterios de inclusión

Incluir datasets que tengan:

```text
- conteos crudos disponibles
- controles no perturbados
- identificación clara de la perturbación
- metadatos de célula, tejido, enfermedad o línea celular
- al menos 50–100 células por condición relevante
- anotación de genes compatible con Ensembl/HGNC
- información de batch, estudio, protocolo y plataforma
- descripción experimental interpretable
```

---

## 8. Criterios de exclusión

Excluir datasets con:

```text
- ausencia de controles adecuados
- perturbaciones ambiguas
- baja calidad celular
- anotación incompleta de targets
- sólo datos normalizados sin raw counts
- diseño experimental no interpretable
- fuerte confusión entre batch y perturbación
- número insuficiente de células por condición
- ausencia de metadatos críticos
```

---

## 9. Estructura del repositorio

```text
project/
├── README.md
├── data_raw/
├── data_intermediate/
├── data_processed/
├── metadata/
├── notebooks/
├── scripts/
├── models/
├── results/
├── figures/
├── reports/
├── references/
├── environment/
│   ├── environment.yml
│   ├── requirements.txt
│   └── Dockerfile
└── workflow/
    ├── Snakefile
    └── config.yaml
```

---

## 10. Formato interno de datos

Se recomienda utilizar:

```text
AnnData / MuData / Zarr / TileDB-SOMA
```

### 10.1 Esquema mínimo de `obs`

```text
cell_id
dataset_id
study_id
batch_id
platform
species
tissue
cell_type
cell_state
disease_status
donor_id
cell_line
perturbation_type
perturbation_target
perturbation_id
drug_name
drug_dose
drug_time
control_status
guide_id
MOI
n_counts
n_genes
percent_mito
percent_ribo
```

### 10.2 Esquema mínimo de `var`

```text
gene_id
gene_symbol
ensembl_id
chromosome
start
end
gene_type
highly_variable
```

### 10.3 Modalidades

```text
RNA:       X / layers["counts"]
ATAC:      peaks / gene activity matrix
CITE-seq:  protein abundance
Drug:      SMILES / chemical descriptors / drug embeddings
GRN:       TF-target / enhancer-target graph
```

---

## 11. Control de calidad

### 11.1 QC por dataset

Realizar primero control de calidad por estudio, no global.

Evaluar:

```text
- genes detectados por célula
- UMIs por célula
- porcentaje mitocondrial
- porcentaje ribosomal
- células vacías
- doublets
- multiplets
- células multiperturbadas no deseadas
- eficiencia de knockdown/activación
```

### 11.2 QC por perturbación

Para cada perturbación:

```text
- número de células perturbadas
- número de controles comparables
- distribución de UMIs
- distribución de genes detectados
- cambio de expresión del gen target
- consistencia entre guías
- separación control vs perturbado
```

### 11.3 QC de confusión experimental

Evaluar si existe colinealidad entre:

```text
perturbación ≈ batch
perturbación ≈ día experimental
perturbación ≈ línea celular
perturbación ≈ profundidad de secuenciación
perturbación ≈ donador
```

Si una perturbación sólo aparece en un batch, no se puede distinguir efecto biológico de efecto técnico.

---

## 12. Normalización e integración

### 12.1 Mantener dos representaciones

```text
Matriz A: raw counts
    Uso: expresión diferencial, pseudobulk, modelos probabilísticos

Matriz B: normalizada/integrada
    Uso: embeddings, clustering, visualización, transferencia
```

### 12.2 Métodos sugeridos

```text
- Scanpy
- scVI
- totalVI
- scANVI
- Harmony
- BBKNN
- scArches
- Seurat v5
```

### 12.3 Advertencia metodológica

No se debe sobrecorregir batch si la señal de perturbación está parcialmente alineada con batch. Una integración demasiado agresiva puede borrar el efecto biológico que se quiere estudiar.

---

## 13. Cálculo de efectos perturbacionales

Para cada dataset, tipo celular y perturbación:

```text
control cells vs perturbed cells
```

Calcular:

```text
- logFC por gen
- p-value
- FDR
- effect size
- Δexpression vector
- genes diferencialmente expresados
- pathway enrichment
- TF activity shift
- cell-state shift
```

### 13.1 Representación matemática básica

```text
Δy = expression_perturbed - expression_control
```

### 13.2 Niveles de análisis

```text
- célula individual
- pseudobulk
- perturbación
- tipo celular
- estudio
- enfermedad
- línea celular
```

Se recomienda usar pseudobulk para inferencia estadística robusta y single-cell para estudiar heterogeneidad de respuesta.

---

## 14. Metaanálisis estadístico

### 14.1 Metaanálisis por gen

Para perturbaciones repetidas entre estudios:

```text
gene_i ~ perturbation + covariates
```

Combinar efectos con:

```text
- random-effects meta-analysis
- inverse-variance weighting
- Stouffer Z-score
- Fisher combined p-value
- Bayesian hierarchical model
```

### 14.2 Heterogeneidad

Calcular:

```text
- I²
- τ²
- Cochran Q
- leave-one-study-out sensitivity
```

Interpretación:

```text
bajo I²  → respuesta conservada
alto I²  → respuesta contexto-dependiente
```

### 14.3 Salidas esperadas

```text
- conserved perturbation signatures
- cell-type-specific signatures
- disease-specific signatures
- drug-specific signatures
- context-dependent failures
```

---

## 15. Inferencia regulatoria

### 15.1 RNA-only

Herramientas posibles:

```text
- SCENIC
- GRNBoost2
- GENIE3
- PIDC
```

### 15.2 Multiome RNA+ATAC

Herramientas posibles:

```text
- SCENIC+
- CellOracle
- LINGER
- Cicero
- chromVAR
- ArchR
- Signac
```

### 15.3 Metarred regulatoria

Construir una red consenso:

```text
TF → enhancer/promoter → target gene
```

Asignar pesos según:

```text
- evidencia de motivo
- accesibilidad ATAC
- correlación enhancer-gene
- cambio tras perturbación
- conservación entre datasets
- soporte por bases externas
```

### 15.4 Validación in silico

Para cada TF perturbado:

```text
1. identificar targets predichos
2. evaluar si los targets cambian tras perturbación
3. comparar dirección esperada vs observada
4. calcular precisión, recall y enrichment de targets
```

---

## 16. Benchmark de modelos predictivos

### 16.1 Baselines simples

```text
- mean baseline
- linear regression
- ridge regression
- elastic net
- kNN
- random forest
- XGBoost
- scVI latent linear model
```

### 16.2 Modelos especializados

```text
- CPA
- chemCPA
- GEARS
- scGen
- scGPT fine-tuned
- Geneformer fine-tuned
- STATE-like predictor
- SSM/Mamba experimental model
- diffusion-based perturbation predictor
```

### 16.3 Splits estrictos

No usar sólo train/test aleatorio.

Usar:

```text
Split 1: perturbación vista, células vistas
Split 2: perturbación no vista
Split 3: tipo celular no visto
Split 4: fármaco no visto
Split 5: dataset no visto
Split 6: combinación no vista
Split 7: enfermedad no vista
```

El split más informativo es:

```text
unseen perturbation + unseen cellular context
```

### 16.4 Métricas

```text
- Pearson sobre Δexpression
- Spearman sobre Δexpression
- RMSE logFC
- DE gene overlap
- precision@k para genes diferenciales
- recall@k
- sign accuracy
- pathway overlap
- TF-target recovery
- calibration error
- uncertainty coverage
```

### 16.5 Control obligatorio

Comparar siempre contra baselines simples. Si un modelo complejo no supera modelos simples bajo splits estrictos, su utilidad predictiva debe interpretarse con cautela.

---

## 17. Evaluación de causalidad regulatoria

### 17.1 Pregunta central

```text
¿Los modelos aprenden regulación causal o sólo coexpresión?
```

### 17.2 Pruebas

```text
1. TF perturbado → ¿cambian sus targets esperados?
2. Genes con motivo accesible → ¿responden más que genes sin motivo?
3. Enhancers accesibles → ¿predicen mejor respuesta?
4. Edges GRN → ¿mejoran predicción frente a correlación?
5. Atención/embeddings → ¿recuperan TF-target reales?
```

### 17.3 Comparaciones

```text
- correlation network
- motif-only network
- ATAC-only network
- RNA-only GRN
- RNA+ATAC GRN
- foundation model attention
- SSM/Mamba embeddings
- explicit GRN-constrained model
```

---

## 18. Módulo farmacológico

### 18.1 Pregunta principal

```text
¿Puede una perturbación farmacológica revertir un estado celular patológico?
```

### 18.2 Firmas

```text
disease_signature = disease - healthy
drug_signature = drug_treated - untreated
```

### 18.3 Reversal score

```text
reversal_score = -correlation(disease_signature, drug_signature)
```

### 18.4 Reversal score estratificado

No usar únicamente una firma global.

Calcular:

```text
reversal_score(cell_type, disease, drug, dose, time)
```

### 18.5 Salidas esperadas

```text
- ranking de fármacos
- ranking de genes diana
- ranking de combinaciones
- tipos celulares sensibles
- tipos celulares resistentes
- mecanismo regulatorio propuesto
- incertidumbre de la predicción
```

---

## 19. Metamodelo final

Después del metaanálisis y benchmark, entrenar un modelo propio.

### 19.1 Entrada

```text
basal cell state
+ perturbation token
+ cell type token
+ disease token
+ tissue token
+ optional ATAC/regulatory tokens
```

### 19.2 Backbone

Opciones:

```text
- SSM/Mamba
- Transformer compacto
- híbrido SSM-Transformer
- encoder variacional con prior regulatorio
```

### 19.3 Priors regulatorios

```text
- TF-target graph
- enhancer-gene links
- pathway graph
- motif accessibility
- chromatin accessibility
```

### 19.4 Salidas

```text
- Δexpression
- genes diferencialmente expresados
- TF activity shift
- pathway shift
- drug sensitivity score
- disease reversal score
- uncertainty score
```

### 19.5 Justificación de SSM/Mamba

La inclusión de modelos de estado-espacio debe tratarse como hipótesis empírica, no como afirmación previa. La pregunta es si arquitecturas eficientes para secuencias largas y representaciones regulatorias superan Transformers o modelos clásicos en escenarios de generalización difíciles.

---

## 20. Incertidumbre

Cada predicción debe reportar incertidumbre.

Métodos posibles:

```text
- deep ensembles
- Monte Carlo dropout
- Bayesian last layer
- Gaussian process sobre embeddings
- conformal prediction
- bootstrap por estudio
```

Salida esperada:

```text
predicción
+ intervalo de confianza
+ confianza por gen
+ confianza por pathway
+ confianza por fármaco
```

La incertidumbre es esencial para convertir el modelo en una herramienta de priorización experimental.

---

## 21. Figuras principales del artículo

### Figura 1 — Diseño del metaanálisis

```text
datasets → modalidades → perturbaciones → células → integración
```

### Figura 2 — Atlas de respuestas

```text
latent space / UMAP
cell types
perturbations
disease states
```

### Figura 3 — Conservación vs contexto

```text
heatmap de perturbaciones
I² por gen/pathway
firmas conservadas
firmas específicas
```

### Figura 4 — Redes regulatorias consenso

```text
TF → enhancer → target
comparación RNA-only vs RNA+ATAC
```

### Figura 5 — Benchmark de modelos

```text
CPA vs GEARS vs scGPT vs Geneformer vs SSM/Mamba vs baselines
```

### Figura 6 — Causalidad regulatoria

```text
TF perturbado
targets predichos
targets observados
recuperación de causalidad
```

### Figura 7 — Priorización farmacológica

```text
disease reversal score
drug ranking
cell-type specificity
uncertainty
```

---

## 22. Producto mínimo publicable

### Título en inglés

```text
A cross-study single-cell perturbation meta-analysis reveals context-dependent regulatory logic and benchmarks virtual-cell models
```

### Título en español

```text
Metaanálisis single-cell de perturbaciones revela lógica regulatoria dependiente del contexto y evalúa modelos virtuales de célula
```

### Contribuciones principales

```text
1. Curación y armonización multi-estudio de perturbaciones single-cell.
2. Metaanálisis de respuestas transcriptómicas conservadas y contexto-dependientes.
3. Benchmark estricto de modelos de predicción perturbacional.
4. Evaluación de causalidad regulatoria usando perturbaciones públicas.
5. Reconstrucción de redes reguladoras consenso RNA+ATAC.
6. Priorización in silico de perturbaciones farmacológicas.
```

---

## 23. Cronograma realista

### Meses 0–2: protocolo y datasets

```text
- definir dominio inicial
- registrar protocolo de búsqueda
- descargar datasets iniciales
- construir esquema AnnData/MuData
- definir criterios de inclusión/exclusión
```

### Meses 2–4: armonización y QC

```text
- normalizar metadatos
- mapear genes
- integrar controles
- filtrar células y perturbaciones
- generar reporte QC por dataset
```

### Meses 4–6: firmas perturbacionales

```text
- calcular DE por perturbación
- calcular Δexpression
- calcular pathway enrichment
- calcular response scores
- construir matriz perturbation × gene × context
```

### Meses 6–9: metaanálisis

```text
- random-effects meta-analysis
- análisis de heterogeneidad
- firmas conservadas
- firmas contexto-dependientes
- validación leave-one-study-out
```

### Meses 9–12: redes regulatorias

```text
- inferencia SCENIC+/LINGER/CellOracle
- red consenso TF-target
- validación con CRISPRi/a
- análisis de causalidad regulatoria
```

### Meses 12–15: benchmark predictivo

```text
- entrenar baselines
- evaluar CPA/GEARS/scGPT/Geneformer
- entrenar modelo propio compacto
- comparar splits estrictos
```

### Meses 15–18: fármacos y manuscrito

```text
- disease reversal
- ranking farmacológico
- análisis de incertidumbre
- figuras finales
- preprint/artículo
```

---

## 24. Infraestructura mínima

Para evitar entrenar un modelo fundacional desde cero:

```text
CPU: 32–64 cores
RAM: 256–512 GB
GPU: 1–4 GPUs de 24–80 GB VRAM
Storage: 10–50 TB
Formato: Zarr / TileDB-SOMA / AnnData backed mode
```

Para datasets grandes, trabajar con:

```text
- streaming
- backed mode
- Zarr
- TileDB-SOMA
- subconjuntos estratificados
- embeddings precomputados
```

No cargar datasets masivos completos en memoria RAM.

---

## 25. Software recomendado

### Lenguaje principal

```text
Python
R para metaanálisis estadístico si es necesario
```

### Ecosistema single-cell

```text
scanpy
anndata
muon
scvi-tools
cellxgene-census
squidpy
decoupler
pertpy
```

### Modelado

```text
pytorch
pytorch-lightning
transformers
xgboost
scikit-learn
gpytorch
```

### Redes regulatorias

```text
pySCENIC
SCENIC+
CellOracle
LINGER
ArchR
Signac
chromVAR
```

### Workflows

```text
Snakemake
Nextflow
Docker
Conda/Mamba
DVC
MLflow/W&B
```

---

## 26. Riesgos y mitigaciones

| Riesgo | Impacto | Mitigación |
|---|---:|---|
| Batch effects dominan la señal | Alto | splits por dataset y batch-aware modeling |
| Perturbaciones no comparables | Alto | ontología de perturbaciones y targets |
| Sobreintegración borra biología | Alto | conservar raw counts y analizar DE por estudio |
| Modelos complejos no superan baselines | Alto | baselines fuertes desde el inicio |
| Causalidad débil | Muy alto | validar contra CRISPRi/a y TF-target recovery |
| Multi-ómica incompleta | Medio-alto | modelos con missing modalities |
| Datos farmacológicos heterogéneos | Alto | estratificar por dosis, tiempo y célula |
| Conclusiones demasiado generales | Alto | restringir dominio inicial |
| Costos computacionales altos | Medio-alto | streaming, Zarr, muestreo estratificado |
| Leakage entre train/test | Muy alto | splits por perturbación, dataset, célula y contexto |

---

## 27. Entregables técnicos

```text
1. Dataset armonizado en AnnData/MuData
2. Tabla maestra de perturbaciones
3. Tabla maestra de metadatos celulares
4. Reportes QC por dataset
5. Matriz perturbation × gene × context
6. Firmas transcriptómicas conservadas
7. Red regulatoria consenso
8. Benchmark de modelos predictivos
9. Ranking farmacológico in silico
10. Manuscrito reproducible
11. Repositorio con pipeline Snakemake/Nextflow
12. Figuras listas para publicación
```

---

## 28. Orden de implementación recomendado

### Paso 1

Seleccionar un dominio restringido:

```text
cáncer / líneas celulares
```

### Paso 2

Construir dataset RNA-only con perturbaciones genéticas y farmacológicas.

### Paso 3

Calcular firmas perturbacionales y hacer metaanálisis.

### Paso 4

Evaluar modelos simples y especializados.

### Paso 5

Añadir redes regulatorias RNA-only.

### Paso 6

Añadir multiome RNA+ATAC.

### Paso 7

Construir modelo propio con priors regulatorios.

### Paso 8

Evaluar reversión farmacológica.

### Paso 9

Preparar manuscrito y repositorio reproducible.

---

## 29. Criterios de éxito

El proyecto debe considerarse exitoso si demuestra al menos uno de los siguientes puntos:

```text
1. identificación de respuestas perturbacionales conservadas entre estudios
2. detección de respuestas dependientes de contexto celular
3. benchmark que revele límites claros de modelos actuales
4. mejora de predicción al incorporar priors regulatorios
5. recuperación significativa de TF-targets bajo perturbaciones CRISPRi/a
6. ranking farmacológico plausible para reversión de estados patológicos
7. modelo reproducible que supere baselines simples bajo splits estrictos
```

---

## 30. Conclusión estratégica

Este proyecto puede ser altamente innovador sin generar datos experimentales propios si se ejecuta como:

```text
metaanálisis computacional
+ benchmark reproducible
+ integración regulatoria
+ predicción perturbacional
+ priorización farmacológica
```

La clave es no limitarse a recopilar datasets. El valor científico surge de probar, con controles estrictos, si los modelos actuales realmente aprenden lógica regulatoria transferible o si sólo capturan correlaciones transcriptómicas dependientes del contexto.

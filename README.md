# MBM_G8
## PROYECTO: Ensamblaje de *novo* y validación taxonómica del bacteriófago T4 como alternativa biológica para el control de cepas resistentes de *Escherichia coli* 

## INTEGRANTES:
* Edison Gustavo Agualema Valdéz
* Edisson Xavier Balarezo Cambi
* Christian Andrés Paredes de la Cueva
* Tatiana Estefania Pillco Encalada
  
## OBJETIVO GENERAL:  
Realizar el ensamblaje de *novo* y la validación taxonómica del bacteriófago T4 mediante herramientas bioinformáticas, con el fin de evaluar su potencial aplicación como alternativa biológica para el control de cepas resistentes de *Escherichia coli*.

## CONJUNTO DE DATOS (DATASET):  
### Datos de Secuenciación de Lecturas Crudas (Raw Reads)
Se utilizaron lecturas crudas depositadas en el Sequence Read Archive (SRA) del NCBI, las cuales representan la base experimental para el ensamblaje de *novo* del genoma viral.  
•	**Identificador de Acceso (SRA):** DRR817419.  
•	**Plataforma de Secuenciación:** Illumina NovaSeq 6000.  
•	**Estrategia de Librería:** WGS (Whole Genome Sequencing).  
•	**Configuración de Lecturas:** Paired-end (Lecturas emparejadas).  
•	**Volumen de Datos Crudos:** 1.5 G bases, con un tamaño de archivo comprimido de 457.8 MB.  
•	**Importancia Técnica:** La utilización de la plataforma NovaSeq 6000 garantiza una alta fidelidad en las lecturas, lo que permite un pre-procesamiento riguroso y un ensamblaje de alta calidad.

### Genoma de Referencia (Gold Standard)
Con el fin de evaluar la precisión del ensamblaje generado y realizar la asignación taxonómica definitiva, se empleó la secuencia genómica de referencia oficial del NCBI.

- **Organismo:** *Escherichia virus T4* (Fago T4)
- **Número de Acceso (RefSeq):** NC_000866.4.
- **Base de datos:** NCBI Nucleotide
- **Longitud de Referencia:** 168,903 bp.
- **Topología:** ADN lineal de doble cadena (dsDNA).
- **Enlace directo:** https://www.ncbi.nlm.nih.gov/nuccore/NC_000866.4

**Aplicación:** Esta secuencia actúa como el control positivo para la validación y la identificación taxonómica mediante herramientas de alineamiento local (BLAST).

## FLUJO DE TRABAJO:  
<img width="3000" height="1688" alt="Ensamblaje de novo y validación taxonómica del (1)_page-0001" src="https://github.com/user-attachments/assets/f9951084-5593-4921-829d-ba0ed5e7b7eb" />





## RESULTADOS:  
El análisis bioinformático mixto (Galaxy-Terminal) del dataset DRR817419 permitió el ensamblaje de novo y la validación taxonómica del bacteriófago T4. Tras la curación de lecturas con Trimmomatic (97.25% de bases > Q30), el ensamblaje preliminar en Galaxy generó un set de 88 scaffolds, destacando un bloque mayoritario de 327,290 pb (176.37X).   
El alineamiento en NCBI BLASTn de este contig reveló una identidad del 99.90% con *Escherichia coli*, confirmando una co-secuenciación masiva del hospedero bacteriano. Para eliminar este ruido genómico, se aplicó un pipeline de exclusión en terminal con Bowtie 2, segregando el ADN bacteriano y aislando 28,104 lecturas virales (0.60% del dataset).  
El re-ensamblaje de este set depurado en SPAdes Terminal resolvió con éxito un único scaffold unificado (NODE_1) de 168,129 pb con una cobertura de 24.64X. Finalmente, la validación global en NCBI BLASTn ratificó un Query Cover del 100%, E-value de 0.0 e identidad nucleotídica del 99.98% con la secuencia de referencia de *Escherichia* virus T4 (NC_000866.4). Estos resultados certifican la máxima pureza del genoma viral obtenido y validan la eficiencia del flujo de trabajo para caracterizar fagos con potencial aplicación biotecnológica y terapéutica.


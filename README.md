# MBM_G8
## PROYECTO: Ensamblaje de *novo* y validación taxonómica del bacteriófago T4 como alternativa biológica para el control de cepas resistentes de *Escherichia coli* 

## INTEGRANTES:
* Edison Gustavo Agualema Valdéz
* Edisson Xavier Balarezo Cambi
* Christian Andrés Paredes de la Cueva
* Tatiana Estefania Pillco Encalada
  
## OBJETIVO GENERAL:  
Realizar el ensamblaje de *novo* y la validación taxonómica del bacteriófago T4 mediante herramientas bioinformáticas, con el fin de evaluar su potencial aplicación como alternativa biológica para el control de cepas resistentes de *Escherichia coli*.

## OBJETIVOS ESPECÍFICOS:  
•	Evaluar la calidad de las lecturas del bacteriófago T4 mediante FastQC, identificando parámetros de calidad relevantes para el análisis bioinformático posterior.  
•	Realizar el ensamblaje de *novo* del genoma del bacteriófago T4 utilizando software especializado para reconstruir su secuencia genómica.   
•	Validar taxonómicamente las secuencias ensambladas mediante herramientas de comparación y clasificación molecular.   
•	Analizar la relevancia biológica del bacteriófago T4 como posible alternativa para el control de cepas resistentes de *Escherichia coli*.  

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

## ## FLUJO DE TRABAJO:

```mermaid
graph TD
    %% Estilos Globales
    classDef fase fill:#4b2e83,color:#fff,stroke:#333,stroke-width:2px;
    classDef subfase fill:#fff,stroke:#4b2e83,stroke-width:1px;
    classDef info fill:#f8f9fa,stroke:#ccc,stroke-dasharray: 5 5;

    %% ETAPA 1: OBTENCIÓN
    E1(1. OBTENCIÓN DE DATOS):::fase
    SRA[ID Acceso: DRR317419<br/>Reads Pareados NovaSeq]:::subfase
    REF[Referencia: NC_000866<br/>~169 Kbp]:::info

    E1 --> SRA
    SRA -.-> REF

    %% ETAPA 2: QC Y LIMPIEZA
    E2(2. CONTROL DE CALIDAD Y LIMPIEZA):::fase
    QC[2.1 FastQC:<br/>Diagnóstico de calidad]:::subfase
    TRIM[2.2 Trimmomatic:<br/>Filtrado y remoción de adaptadores]:::subfase
    PARAM[Parámetros:<br/>SLIDINGWINDOW:4:20<br/>MINLEN:50]:::info

    SRA --> E2
    E2 --> QC
    E2 --> TRIM
    TRIM -.-> PARAM

    %% ETAPA 3: ENSAMBLAJE
    E3(3. ENSAMBLAJE DE NOVO):::fase
    SHOV[3.1 Shovill / SPAdes:<br/>Construcción de Contigs]:::subfase
    VER[3.2 Verificación:<br/>Análisis de N50 y Longitud]:::subfase

    TRIM --> E3
    E3 --> SHOV
    SHOV --> VER

    %% ETAPA 4: ENTREGABLES
    E4(4. ENTREGABLES DEL PROYECTO):::fase
    GIT[4.1 Repositorio GitHub:<br/>INFORME.md y Workflow]:::subfase
    DEF[4.2 Defensa Técnica:<br/>Comparativa Galaxy vs CLI]:::subfase

    VER --> E4
    E4 --> GIT
    E4 --> DEF

    %% Layout
    subgraph "Pipeline Bioinformático - Grupo 8"
    E1
    E2
    E3
    E4
    end
```




## RESULTADOS:  



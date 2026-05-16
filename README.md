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
<img width="2753" height="1216" alt="Ensamblaje de novo y validación taxonómica del_page-0001" src="https://github.com/user-attachments/assets/dff40c8d-bcfb-4f67-97fb-cdc31a5823fe" />




## RESULTADOS:  

El análisis bioinformático del dataset DRR817419 permitió obtener un ensamblaje genómico *de novo* de alta calidad mediante la herramienta SPAdes. A partir de las lecturas previamente filtradas y depuradas, se generó un total de **88 scaffolds principales**, evidenciando una reconstrucción genómica adecuada para posteriores análisis taxonómicos y estructurales.

El ensamblaje presentó un **scaffold principal de 327,290 pb**, caracterizado por una elevada cobertura de secuenciación y una alta continuidad ensamblativa, lo que indica una adecuada calidad de las lecturas utilizadas durante el proceso de ensamblaje *de novo*.

Posteriormente, la validación taxonómica realizada mediante BLASTn mostró una **identidad del 99.90% con *Escherichia coli***, confirmando la presencia predominante del hospedero bacteriano dentro del dataset analizado. Este resultado es biológicamente coherente debido a que el bacteriófago T4 infecta específicamente cepas de *E. coli*, por lo que es esperado encontrar material genético bacteriano asociado durante la secuenciación.

Los resultados obtenidos permitieron identificar de manera precisa la relación biológica entre el bacteriófago T4 y su hospedero, proporcionando una base sólida para futuras etapas de filtrado, aislamiento de secuencias virales y validación taxonómica específica del genoma del fago.


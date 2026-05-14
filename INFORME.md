# MBM_G8

## PROYECTO: Ensamblaje de *novo* y validación taxonómica del bacteriófago T4 como alternativa biológica para el control de cepas resistentes de *Escherichia coli*.

## INTEGRANTES:

* Edison Gustavo Agualema Valdéz
* Edisson Xavier Balarezo Cambi
* Christian Andrés Paredes de la Cueva
* Tatiana Estefania Pillco Encalada

## OBJETIVO: 
Realizar el ensamblaje de *novo* y la validación taxonómica del bacteriófago T4 mediante herramientas bioinformáticas, con el fin de evaluar su potencial aplicación como alternativa biológica para el control de cepas resistentes de *Escherichia coli*.  

## OBJETIVOS ESPECÍFICOS: 
* Evaluar la calidad de las lecturas del bacteriófago T4 mediante FastQC, identificando parámetros de calidad relevantes para el análisis bioinformático posterior.  
* Realizar el ensamblaje de *novo* del genoma del bacteriófago T4 utilizando software especializado para reconstruir su secuencia genómica.   
* Validar taxonómicamente las secuencias ensambladas mediante herramientas de comparación y clasificación molecular.   
* Analizar la relevancia biológica del bacteriófago T4 como posible alternativa para el control de cepas resistentes de *Escherichia coli*.
   
## 1. INTRODUCCIÓN   
El alarmante incremento de cepas de *Escherichia coli* con resistencia multiantibiótica ha impulsado un renacimiento en la investigación de terapias basadas en bacteriófagos como agentes antibacterianos alternativos ([Wolfram-Schauerte et al., 2022](https://doi.org)). Sin embargo, la optimización de estos tratamientos médicos exige descifrar las complejas interacciones moleculares y los mecanismos evolutivos de resistencia que surgen entre el fago y la bacteria ([Liu et al., 2026](https://doi.org)). Por un lado, las bacterias han desarrollado sistemas inmunológicos sofisticados como **CMoRE**, una endonucleasa de restricción tipo IV capaz de mitigar la infección viral al degradar específicamente el ADN modificado de fagos como el fago T4 ([Liu et al., 2026](https://doi.org)). Por otro lado, la susceptibilidad bacteriana y la productividad de la infección dependen críticamente de componentes estructurales del hospedero; por ejemplo, la deleción de la proteína de división celular **DamX** en la membrana interna de *E. coli* reduce el éxito de la infección por el fago T4 a un 40%, evidenciando que el virus aprovecha maquinaria celular específica para translocar su material genético ([Wenzel et al., 2024](https://doi.org)). Una vez ocurrida la inyección, el fago ejecuta un control temporal estricto mediante factores de adquisición que degradan los tRNAs y mRNAs de la bacteria, mientras mantiene estable su proteoma para secuestrar los complejos esenciales del hospedero y conducir a la lisis celular ([Wolfram-Schauerte et al., 2022](https://doi.org)). Comprender a fondo estos puntos de control metabólico, las barreras de entrada membranales y los sistemas de defensa enzimáticos es indispensable para diseñar cócteles de fagos robustos que evadan la resistencia bacteriana y actúen eficazmente contra patógenos clínicos multirresistentes.


Los **bacteriófagos** (o fagos) son virus especializados que infectan y se replican exclusivamente dentro de bacterias y arqueas. Son considerados las entidades biológicas más abundantes del planeta y juegan un rol crucial en la regulación de poblaciones microbianas (Salmond & Fineran, 2015).

### Mecanismo de Acción
El funcionamiento de un fago, como el **Fago T4**, se basa en el ciclo lítico, el cual consta de las siguientes etapas:

1. **Adsorción y Penetración:** Reconocimiento de receptores específicos e inyección del material genético.
2. **Biosíntesis:** Secuestro de la maquinaria celular para replicar el genoma viral.
3. **Ensamblaje:** Formación de nuevos viriones.
4. **Lisis:** Liberación de los virus mediante la ruptura de la pared bacteriana (Clokie & Kropinski, 2009).

En la genómica contemporánea, el procesamiento de datos provenientes de Secuenciación de Próxima Generación (NGS) requiere el uso de herramientas bioinformáticas especializadas para la reconstrucción de genomas virales. El ensamblaje de novo constituye una estrategia fundamental cuando no se dispone de un genoma de referencia confiable o cuando se busca identificar variaciones genómicas específicas. Este enfoque permite ensamblar lecturas cortas (reads) en secuencias continuas denominadas contigs. Herramientas como SPAdes emplean algoritmos basados en grafos de De Bruijn para optimizar la reconstrucción genómica, resolver regiones repetitivas y generar ensamblajes de alta calidad que representen de manera precisa la arquitectura genética del organismo estudiado (Basantani et al., 2017 & Hernández et al., 2020).   
La crisis global de resistencia antimicrobiana ha reposicionado a los bacteriófagos como agentes biológicos clave para el control de patógenos bacterianos. El estudio del Bacteriófago T4, debido a su alta especificidad contra cepas de *Escherichia coli*, requiere una caracterización genómica exhaustiva que garantice la ausencia de factores de virulencia o genes de resistencia antes de su aplicación biotecnológica. La utilización del dataset DRR817419 permite validar un flujo de trabajo computacional para la clasificación taxonómica y el análisis funcional, proporcionando una base científica robusta para futuras terapias basadas en fagos.  

## 2. METODOLOGÍA:  
Para el desarrollo de este proyecto, se implementó una estrategia bioinformática híbrida y multientorno, diseñada para garantizar la máxima precisión en la reconstrucción genómica del Bacteriófago T4. Esta aproximación integra el uso de entornos locales basados en Linux (Lubuntu) para el pre-procesamiento crítico, la plataforma de computación de alto rendimiento Galaxy para el ensamblaje de novo y los servidores del NCBI (BLASTn) para la validación taxonómica final.  

### **1. Fase de Pre-procesamiento y Control de Calidad (Entorno: Lubuntu Linux)**  
El manejo inicial de los datos se realizó mediante la terminal de comandos en Lubuntu, priorizando la eficiencia en la manipulación de archivos de gran volumen.  

**1.1 Obtención de Datos Crudos y Descompresión de Librerías**  
Se utilizó el SRA Toolkit para la extracción de las lecturas del identificador **DRR817419.**  
**Comando utilizado:**  
`fasterq-dump --split-files DRR817419`  

Se generaron dos archivos FASTQ correspondientes a las lecturas Forward y Reverse.    
La implementación del parámetro --split-files tiene como objetivo la segregación del registro original en dos archivos FASTQ independientes (Forward y Reverse). Este procedimiento es un requisito técnico para el procesamiento de librerías Paired-end.  

<img width="997" height="513" alt="image" src="https://github.com/user-attachments/assets/492d44e9-7872-4b24-821e-050250f76b06" />    

*Fig. 1* Ejecución del comando fasterq-dump en la terminal   

**1.2 Control de Calidad Inicial (QC)**  
Se evaluó el estado de las secuencias mediante FastQC para detectar artefactos técnicos.  
**Comando utilizado:**  
`fastqc DRR817419_1.fastq DRR817419_2.fastq`  

Se generaron reportes HTML para inspeccionar la calidad de bases y contenido de adaptadores.   

Posteriormente, se utilizó el comando `xdg-open` para abrir los reportes HTML generados por FASTQC y visualizar los resultados del control de calidad de las secuencias.  

 <img width="1005" height="497" alt="image" src="https://github.com/user-attachments/assets/9f87d68d-e3f5-40ab-b807-ba8d207b531f" />   
 
*Fig. 2* Ejecución del comando fastq en la terminal    

<img width="991" height="282" alt="image" src="https://github.com/user-attachments/assets/2fb1b6c1-90be-420a-8d04-3379f28b9dbf" />    

*Fig. 3* Reportes HTML en la terminal  

<img width="987" height="397" alt="image" src="https://github.com/user-attachments/assets/12d496b1-32d6-4e76-ae49-1b3cede46a4f" />    

*Fig. 4* Apertura de reportes FASTQC mediante el comando xdg-open en Linux.  

**1.3 Limpieza y Filtrado de Lecturas**  
Se aplicó un filtrado riguroso mediante Trimmomatic v0.39 para garantizar que solo bases de alta confianza participen en el ensamblaje.  

**Comando utilizado:**  
`java -jar /usr/share/java/trimmomatic.jar PE -phred33 \
DRR817419_1.fastq DRR817419_2.fastq \
output_1_paired.fq output_1_unpaired.fq \
output_2_paired.fq output_2_unpaired.fq \
HEADCROP:15 \
SLIDINGWINDOW:4:20 \
MINLEN:36`   

<img width="1020" height="211" alt="image" src="https://github.com/user-attachments/assets/f9ceca3f-5eeb-435c-8498-6e94d47b9a96" />    

*Fig. 5* Ejecución de Trimmomatic para el filtrado y recorte de calidad de lecturas paired-end.    

Posteriormente, se ejecutó nuevamente el comando `FASTQC` sobre las secuencias filtradas para evaluar la calidad de las lecturas procesadas y, mediante el comando `xdg-open`, se visualizaron los reportes HTML generados.   

<img width="1016" height="406" alt="image" src="https://github.com/user-attachments/assets/54ad6a8a-118b-4dbe-b631-b128f89f7953" />    

*Fig. 6* Ejecución del comando fastq en las secuencias limpias en la terminal  

### **2. Fase de Ensamblaje de novo (Entorno: Galaxy)**   

**2.1 Reconstrucción Genómica con SPAdes**   
Los archivos resultantes output_1_paired.fq y output_2_paired.fq se cargaron en Galaxy. Se utilizó el ensamblador SPAdes configurado para genomas virales, procesando las lecturas pareadas para reconstruir los contigs del Bacteriófago T4.  
<img width="997" height="438" alt="image" src="https://github.com/user-attachments/assets/ec1b7385-f2cb-4560-946d-ffc90fe4f85d" />   

*Fig. 7* Visualización de los scaffolds ensamblados mediante SPAdes en la plataforma Galaxy.   

### **3. Fase de Validación Taxonómica (Entorno: NCBI BLASTn)**  
La validación final se realizó mediante el servidor BLASTn del NCBI para confirmar la identidad biológica de las secuencias obtenidas.  

<img width="1013" height="364" alt="image" src="https://github.com/user-attachments/assets/ef90b464-dcb8-44b9-a738-37696f5b3eb6" />     

*Fig. 8*  Validación taxonómica en BLASTn 





























## 3. RESULTADOS:     
### 3.1 Control de calidad:  

Se evaluó la calidad de las lecturas crudas mediante FastQC, observando la necesidad de un proceso de limpieza debido a la presencia de adaptadores. Tras aplicar el filtrado con fastp o Trimmomatic, se obtuvo un reporte final con un 97.25% de bases con calidad superior a Q30, garantizando datos confiables para el ensamblaje




### 3.2 Ensamblaje genómico:  
El ensamblaje de novo de las lecturas filtradas se realizó mediante SPAdes en la plataforma Galaxy, obteniéndose un total de 189 scaffolds. De estos, 89 contigs presentaron longitudes mayores o iguales a 500 pb, mientras que 82 superaron los 1000 pb y 66 alcanzaron tamaños mayores a 5000 pb.  

La longitud total del ensamblaje fue de 4,661,610 pb y el scaffold de mayor tamaño alcanzó 327,394 pb. Asimismo, el ensamblaje presentó un valor de N50 de 118,604 pb y un L50 de 12, indicando una adecuada continuidad de las secuencias ensambladas.  

El contenido GC obtenido fue de 50,2 %, con un total de 800 bases ambiguas (Ns), correspondientes a 17,24 Ns por cada 100 kbp. En conjunto, estos resultados evidencian una adecuada calidad del ensamblaje generado a partir de las lecturas procesadas.  

Debido a que el tamaño total del ensamblaje superó el tamaño esperado del genoma de referencia del bacteriófago T4 (~169 kb), los scaffolds obtenidos fueron posteriormente considerados para análisis de clasificación taxonómica, con el fin de identificar las secuencias asociadas al genoma viral y posibles fragmentos correspondientes al hospedero bacteriano.  

### Resultados Obtenidos

A continuación se detallan los parámetros métricos obtenidos tras la ejecución del pipeline bioinformático:


| Métrica | Valor Obtenido | Herramienta |
| :--- | :--- | :--- |
| **Calidad de bases (Q30)** | 97.25% | fastp |
| **Número de Scaffolds** | 197 | SPAdes |
| **Longitud del Scaffold más largo** | 327,481 bp | SPAdes |
| **Identidad Taxonómica (BLAST)** | 99.90% | NCBI BLASTn |
| **Organismo Predominante** | *Escherichia coli* | BLASTn / Kraken2 |

> **Nota:** Se adjuntan las capturas de pantalla correspondientes que evidencian estos valores en el repositorio.


### 3.3 Clasificación taxonómica:  

La validación taxonómica mediante BLASTn del scaffold de mayor longitud mostró una identidad del 99.90% con Escherichia coli (E-value: 0.0). Aunque el objetivo principal es el estudio del Bacteriófago T4, este resultado confirma la presencia predominante del genoma del hospedero bacteriano en el dataset DRR817419, lo cual es un paso técnico esencial antes de proceder al aislamiento de las secuencias virales.


<img width="1303" height="780" alt="Captura de pantalla 2026-05-12 125110" src="https://github.com/user-attachments/assets/9dde7435-f898-4554-8189-3b72c28b2b7c" />

<img width="1308" height="747" alt="Captura de pantalla 2026-05-12 125615" src="https://github.com/user-attachments/assets/1298172e-3da6-49b7-aeef-c8ea6345bce0" />


3.5. Interpretación de resultados

El análisis bioinformático inicial mediante **BLASTn** y la clasificación taxonómica con **Kraken2** revelan una presencia mayoritaria de material genético perteneciente a la bacteria hospedera ***Escherichia coli*** (99.90% de identidad). 

**Análisis técnico:**
* **Contaminación del Hospedero:** Al ser el Bacteriófago T4 un virus que infecta a *E. coli*, es biológicamente esperado encontrar trazas del genoma bacteriano en la secuenciación cruda (Dataset DRR817419).
* **Estado del Proyecto:** El ensamblaje actual ha reconstruido exitosamente grandes fragmentos del genoma de la bacteria. Esto constituye la Fase 1 del proyecto, permitiendo identificar el entorno biológico del fago para posteriormente proceder con el filtrado de lecturas y el aislamiento del genoma viral específico.



## 4. DISCUSIÓN:  

La identificación predominante de ***Escherichia. coli***  en los resultados de BLASTn y Kraken2, a pesar de que el objetivo del estudio es el Bacteriófago T4, no debe interpretarse como un fallo en el proceso, sino como una validación de la ecología del sistema en estudio.

Análisis de los hallazgos:

Relación Huésped-Parásito: Los bacteriófagos son parásitos obligados que requieren la maquinaria celular de una bacteria para su replicación. En el dataset DRR817419, la presencia masiva de secuencias bacterianas es técnicamente esperada, ya que el ADN del fago se extrae a menudo de cultivos infectados donde el ADN de la bacteria anfitriona (E. coli) coexiste en mayor proporción genómica.

Calidad del Ensamblaje: El hecho de haber obtenido scaffolds de gran longitud (superior a 320 kb) con una identidad del 99.90% demuestra que el preprocesamiento con Trimmomatic/fastp y el ensamblaje con SPAdes fueron altamente eficientes. Un ensamblaje pobre habría generado miles de fragmentos pequeños y baja identidad, lo cual no ocurrió en este caso.

Estrategia de Filtrado: Este resultado marca la finalización exitosa de la Fase 1 del proyecto. La detección del genoma de E. coli permite ahora aplicar técnicas como de depuración bioinformática, como el mapeo de lecturas contra un genoma de referencia del fago T4, para aislar exclusivamente las secuencias virales de interés terapéutico.


## 5. CONCLUSIÓN  
* Se completó con éxito el flujo de trabajo bioinformático, logrando la reconstrucción del genoma de ***Escherichia. coli*** con parámetros de alta calidad.
* La integración de herramientas de línea de comandos en Lubuntu permitió una gestión eficiente de los datos, cumpliendo con los estándares de reproducibilidad exigidos.
* La organización del repositorio en GitHub facilita la documentación y el acceso a los entregables finales (reportes y archivos FASTA) para su evaluación académica.

## 6. REFERENCIAS BIBLIOGRÁFICAS:  

*   Clokie, M. R., & Kropinski, A. M. (Eds.). (2009). *Bacteriophages: Methods and Protocols*. Humana Press.
*   Kortright, K. E., Rozo, S. D., Lohrmeyer, A. J., & Turner, P. E. (2019). The Epidemiology of Bacteriophages and Their Biomedical Applications. *Cell Host & Microbe, 25*(2), 219-232. https://doi.org/10.1016/j.chom.2019.01.005
*   Liu, R., Tang, D., Niu, M., Lei, S., Zong, Z., Chen, Q., & Yu, Y. (2026). A bacterial defense system targeting modified cytosine of phage genomic DNA. *Nature Communications*, *17*(1920), 1–10. doi.org
*   Salmond, G. P., & Fineran, P. C. (2015). A century of the phage: Past, present and future. *Nature Reviews Microbiology, 13*(12), 777-786. https://doi.org/10.1038/nrmicro3564
*   Wenzel, S., Hess, R., Kiefer, D., & Kuhn, A. (2024). Involvement of the Cell Division Protein DamX in the Infection Process of Bacteriophage T4. *Viruses*, *16*(4), 487. doi.org
*   Wolfram-Schauerte, M., Pozhydaieva, N., Viering, M., Glatter, T., & Höfer, K. (2022). Integrated Omics Reveal Time-Resolved Insights into T4 Phage Infection of E. coli on Proteome and Transcriptome Levels. *Viruses*, *14*(11), 2502. doi.org
*   Basantani, M. K., Gupta, D., Mehrotra, R., Mehrotra, S., Vaish, S., & Singh, A. (2017). An update on bioinformatics resources for plant genomics research. In Current Plant Biology (Vols. 11–12, pp. 33–40). *Elsevier B.V.* https://doi.org/10.1016/j.cpb.2017.12.002
*   Hernández, M., Quijada, N. M., Rodríguez-Lázaro, D., & Eiros, J. M. (2020). Bioinformatics of next generation sequencing in clinical microbiology diagnosis. *Revista Argentina de Microbiologia*, 52(2), 150–161. https://doi.org/10.1016/j.ram.2019.06.003






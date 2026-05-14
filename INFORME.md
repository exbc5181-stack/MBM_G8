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

---

## 2. METODOLOGÍA:  
Para el desarrollo de este proyecto, se implementó una estrategia bioinformática híbrida y multientorno, diseñada para garantizar la máxima precisión en la reconstrucción genómica del Bacteriófago T4. Esta aproximación integra el uso de entornos locales basados en Linux (Lubuntu) para el pre-procesamiento crítico, la plataforma de computación de alto rendimiento Galaxy para el ensamblaje de novo y los servidores del NCBI (BLASTn) para la validación taxonómica final.  

**1. Fase de Pre-procesamiento y Control de Calidad (Entorno: Lubuntu Linux)**  
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

Posteriormenete, se utilizó el comando `xdg-open` para abrir los reportes HTML generados por FASTQC y visualizar los resultados del control de calidad de las secuencias.  

 <img width="1005" height="497" alt="image" src="https://github.com/user-attachments/assets/9f87d68d-e3f5-40ab-b807-ba8d207b531f" />   
 
*Fig. 2* Ejecución del comando fastq en la terminal    

<img width="958" height="282" alt="image" src="https://github.com/user-attachments/assets/0f2a8dd7-6dc4-4fe4-b4c7-7a100c0d6cb9" />   

*Fig. 3* Reportes HTML en la terminal  

<img width="987" height="397" alt="image" src="https://github.com/user-attachments/assets/12d496b1-32d6-4e76-ae49-1b3cede46a4f" />    

*Fig. 4* Apertura de reportes FASTQC mediante el comando xdg-open en Linux.  

**1.3 Limpieza y Filtrado de Lecturas (Trimming)**  
Se aplicó un filtrado riguroso mediante Trimmomatic v0.39 para garantizar que solo bases de alta confianza participen en el ensamblaje.  
**Comando utilizado:**  
`java -jar trimmomatic.jar PE DRR817419_1.fastq DRR817419_2.fastq \
output_1_paired.fq output_1_unpaired.fq \
output_2_paired.fq output_2_unpaired.fq \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 HEADCROP:15 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36`   

<img width="1020" height="408" alt="image" src="https://github.com/user-attachments/assets/5b973129-a329-4c81-8d42-23fc0d579b13" />  

*Fig. 5* Ejecución de Trimmomatic para el filtrado y recorte de calidad de lecturas paired-end.  

















## 3. RESULTADOS:   
## 4. DISCUSIÓN:  
## 5. CONCLUSIÓN  
## 6. REFERENCIAS BIBLIOGRÁFICAS:  

*   Clokie, M. R., & Kropinski, A. M. (Eds.). (2009). *Bacteriophages: Methods and Protocols*. Humana Press.
*   Kortright, K. E., Rozo, S. D., Lohrmeyer, A. J., & Turner, P. E. (2019). The Epidemiology of Bacteriophages and Their Biomedical Applications. *Cell Host & Microbe, 25*(2), 219-232. https://doi.org/10.1016/j.chom.2019.01.005
*   Liu, R., Tang, D., Niu, M., Lei, S., Zong, Z., Chen, Q., & Yu, Y. (2026). A bacterial defense system targeting modified cytosine of phage genomic DNA. *Nature Communications*, *17*(1920), 1–10. doi.org
*   Salmond, G. P., & Fineran, P. C. (2015). A century of the phage: Past, present and future. *Nature Reviews Microbiology, 13*(12), 777-786. https://doi.org/10.1038/nrmicro3564
*   Wenzel, S., Hess, R., Kiefer, D., & Kuhn, A. (2024). Involvement of the Cell Division Protein DamX in the Infection Process of Bacteriophage T4. *Viruses*, *16*(4), 487. doi.org
*   Wolfram-Schauerte, M., Pozhydaieva, N., Viering, M., Glatter, T., & Höfer, K. (2022). Integrated Omics Reveal Time-Resolved Insights into T4 Phage Infection of E. coli on Proteome and Transcriptome Levels. *Viruses*, *14*(11), 2502. doi.org





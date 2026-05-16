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
* Evaluar la calidad de las lecturas del dataset, identificando parámetros y sesgos analíticos para el acondicionamiento de los datos.
* Realizar el ensamblaje de *novo* dirigido del genoma del bacteriófago T4 utilizando SPAdes para reconstruir su secuencia de manera continua.
* Validar taxonómicamente las secuencias ensambladas finales para confirmar la identidad molecular y exactitud del genoma viral obtenido.
   
## 1. INTRODUCCIÓN   

## 1.1 CONTEXTO BIOLÓGICO Y MECANISMO DE ACCIÓN
La crisis global de resistencia antimicrobiana ha reposicionado a los bacteriófagos (o fagos) como agentes biológicos clave para el control de patógenos bacterianos. Los fagos son virus especializados que infectan y se replican exclusivamente dentro de bacterias y arqueas; son considerados las entidades biológicas más abundantes del planeta y juegan un rol crucial en la regulación de poblaciones microbianas *(Salmond & Fineran, 2015)*.

El estudio del Bacteriófago T4, debido a su alta especificidad contra cepas de *Escherichia coli*, requiere una caracterización genómica exhaustiva que garantice la ausencia de factores de virulencia o genes de resistencia antes de su aplicación biotecnológica. El funcionamiento de este virus se basa en el ciclo lítico clásico, el cual consta de cuatro etapas fundamentales *(Clokie & Kropinski, 2009)*:

```mermaid
graph LR
    %% Estilos de Nodos
    classDef etapa fill:#4b2e83,color:#fff,stroke:#333,stroke-width:2px;
    
    E1(1. Adsorción y Penetración):::etapa --> E2(2. Biosíntesis):::etapa
    E2 --> E3(3. Ensamblaje):::etapa
    E3 --> E4(4. Lisis Celular):::etapa
```

1. **Adsorción y Penetración:** Reconocimiento de receptores específicos e inyección del material genético.
2. **Biosíntesis:** Secuestro de la maquinaria celular para replicar el genoma viral.
3. **Ensamblaje:** Formación de nuevos viriones.
4. **Lisis:** Liberación de los virus mediante la ruptura de la pared bacteriana (Clokie & Kropinski, 2009).

## 1.2 INTERACCIONES MOLECULARES Y EVOLUCIÓN DE LA RESISTENCIA

El alarmante incremento de cepas de *Escherichia coli* con resistencia multiantibiótica ha impulsado un renacimiento en la investigación de estas terapias basadas en bacteriófagos como agentes antibacterianos alternativos *(Wolfram-Schauerte et al., 2022)*. Sin embargo, la optimización de estos tratamientos médicos exige descifrar las complejas interacciones moleculares y los mecanismos evolutivos de resistencia que surgen bidireccionalmente entre el fago y la bacteria *(Liu et al., 2026)*.

```mermaid
graph TD
    %% Configuración de Estilos Globales
    classDef caja fill:#4b2e83,color:#fff,stroke:#333,stroke-width:2px;
    classDef proceso fill:#fff,stroke:#4b2e83,stroke-width:1px;
    classDef consecuencia fill:#fff,stroke:#e74c3c,stroke-width:2px;

    %% Bloque Membrana Interna
    subgraph "E. coli: Membrana Interna"
        A[Proteína DamX]:::caja
        B[Si se elimina / deleta]:::proceso
        C["Éxito de infección del Fago T4 disminuye al 40% (Wenzel et al., 2024)"]:::consecuencia
        
        A --> B
        B --> C
    end

    %% Bloque Sistema Inmune
    subgraph "E. coli: Sistema Inmune"
        D[Endonucleasa CMoRE]:::caja
        E[Infección Activa]:::proceso
        F["Degrada específicamente el ADN modificado del fago (Liu et al., 2026)"]:::consecuencia
        
        D --> E
        E --> F
    end

    %% Ajuste visual de texto de alerta
    style C color:#c0392b
    style F color:#c0392b
```
### 1.2.1 Mecanismos de Defensa Bacteriana

Las bacterias han desarrollado sistemas inmunológicos sofisticados para contrarrestar la agresión viral. Destaca entre ellos **CMoRE**, una endonucleasa de restricción tipo IV capaz de mitigar la infección viral al degradar específicamente el ADN modificado de fagos como el fago T4 *(Liu et al., 2026)*.

Asimismo, la susceptibilidad bacteriana y la productividad de la infección dependen críticamente de componentes estructurales del hospedero. Por ejemplo, se ha evidenciado que la deleción de la proteína de división celular **DamX** en la membrana interna de *Escherichia coli* reduce el éxito de la infección por el fago T4 a un **40%**, demostrando que el virus aprovecha maquinaria celular específica para lograr translocar su material genético *(Wenzel et al., 2024)*.

### 1.2.2 Estrategias de Contradefensa del Fago

Una vez ocurrida con éxito la inyección, el fago ejecuta un control temporal estricto mediante factores de adquisición que degradan los tRNAs y mRNAs de la bacteria, mientras mantiene estable su propio proteoma para secuestrar los complejos esenciales del hospedero y conducir inexorablemente a la lisis celular *(Wolfram-Schauerte et al., 2022)*.

**Implicación Clínica:** Comprender a fondo estos puntos de control metabólico, las barreras de entrada membranales y los sistemas de defensa enzimáticos es indispensable para diseñar cócteles de fagos robustos que evadan la resistencia bacteriana y actúen eficazmente contra patógenos clínicos multirresistentes.

---

## 1.3 METODOLOGÍA BIOINFORMÁTICA Y ANÁLISIS GENÓMICO

En la genómica contemporánea, el procesamiento de datos provenientes de Secuenciación de Próxima Generación (NGS) requiere el uso de herramientas bioinformáticas especializadas para la reconstrucción de genomas virales. La utilización del dataset **DRR317419** permite validar un flujo de trabajo computacional para la clasificación taxonómica y el análisis funcional, proporcionando una base científica robusta para futuras terapias basadas en fagos.

El **ensamblaje *de novo*** constituye una estrategia fundamental cuando no se dispone de un genoma de referencia confiable o cuando se busca identificar variaciones genómicas específicas sin sesgar el alineamiento. Este enfoque permite ensamblar lecturas cortas (*reads*) en secuencias continuas denominadas *contigs* *(Basantani et al., 2017; Hernández et al., 2020)*

Herramientas como **SPAdes** emplean algoritmos basados en **grafos de De Bruijn** para optimizar la reconstrucción genómica, resolver regiones repetitivas y generar ensamblajes de alta calidad que representen de manera precisa la arquitectura genética del organismo estudiado *(Basantani et al., 2017; Hernández et al., 2020)*.

## 2. METODOLOGÍA:  
Para el desarrollo de este proyecto, se implementó una estrategia bioinformática híbrida y multientorno, diseñada para garantizar la máxima precisión en la reconstrucción genómica del Bacteriófago T4. Esta aproximación integra el uso de entornos locales basados en Linux (Lubuntu) para el pre-procesamiento crítico, la plataforma de computación de alto rendimiento Galaxy para el ensamblaje de novo y los servidores del NCBI (BLASTn) para la validación taxonómica final.  

### **1. Fase de Pre-procesamiento y Control de Calidad (Entorno: Lubuntu Linux)**  
El manejo inicial de los datos se realizó mediante la terminal de comandos en Lubuntu, priorizando la eficiencia en la manipulación de archivos de gran volumen.  

**1.1 Obtención de Datos Crudos y Descompresión de Librerías**  
Se utilizó el SRA Toolkit para la extracción de las lecturas del identificador **DRR817419.**  
**Comando utilizado:**  
```
fasterq-dump --split-files DRR817419
```  

Se generaron dos archivos FASTQ correspondientes a las lecturas Forward y Reverse.    
La implementación del parámetro --split-files tiene como objetivo la segregación del registro original en dos archivos FASTQ independientes (Forward y Reverse). Este procedimiento es un requisito técnico para el procesamiento de librerías Paired-end.  

<img width="997" height="513" alt="image" src="https://github.com/user-attachments/assets/492d44e9-7872-4b24-821e-050250f76b06" />    

*Fig. 1* Ejecución del comando fasterq-dump en la terminal   

**1.2 Control de Calidad Inicial (QC)**  
Se evaluó el estado de las secuencias mediante FastQC para detectar artefactos técnicos.  
**Comando utilizado:**  
```
fastqc DRR817419_1.fastq DRR817419_2.fastq
```

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
```
java -jar /usr/share/java/trimmomatic.jar PE -phred33 \
DRR817419_1.fastq DRR817419_2.fastq \
output_1_paired.fq output_1_unpaired.fq \
output_2_paired.fq output_2_unpaired.fq \
HEADCROP:15 \
SLIDINGWINDOW:4:20 \
MINLEN:36
```   

<img width="1020" height="211" alt="image" src="https://github.com/user-attachments/assets/f9ceca3f-5eeb-435c-8498-6e94d47b9a96" />    

*Fig. 5* Ejecución de Trimmomatic para el filtrado y recorte de calidad de lecturas paired-end.    

Posteriormente, se ejecutó nuevamente el comando `FASTQC` sobre las secuencias filtradas para evaluar la calidad de las lecturas procesadas y, mediante el comando `xdg-open`, se visualizaron los reportes HTML generados.   

<img width="1016" height="406" alt="image" src="https://github.com/user-attachments/assets/54ad6a8a-118b-4dbe-b631-b128f89f7953" />    

*Fig. 6* Ejecución del comando fastq en las secuencias limpias en la terminal  

**2. Fase de Ensamblaje de *novo* y Depuración Genómica**       

**2.1 Reconstrucción Genómica con SPAdes (Galaxy)**     
El ensamblaje de novo preliminar se ejecutó en Galaxy utilizando el algoritmo SPAdes. Como datos de entrada (input), se emplearon las lecturas paired-end de alta calidad previamente depuradas con Trimmomatic. El software procesó estas secuencias limpias para generar el set inicial de scaffolds estructurales destinados a la evaluación.  

<img width="997" height="438" alt="image" src="https://github.com/user-attachments/assets/ec1b7385-f2cb-4560-946d-ffc90fe4f85d" />   

*Fig. 7* Visualización de los scaffolds ensamblados mediante SPAdes en la plataforma Galaxy.   

**2.2 Depuración Genómica con Bowtie 2 (Terminal)**    
Al identificarse co-secuenciación masiva del hospedero bacteriano en Galaxy, el flujo de trabajo se trasladó a entorno de terminal Linux para ejecutar un filtrado por exclusión. Las lecturas previamente limpias se mapearon mediante la herramienta Bowtie 2 contra el genoma de referencia de *Escherichia coli* para segregar el ruido molecular. Las lecturas remanentes, correspondientes al virus, se sometieron directamente a un segundo proceso de ensamblaje de novo en SPAdes Terminal para generar el scaffold definitivo.  

**Comando utilizado**  

```
PEGA AQUI CRIS TU COMANDO   
```  

### **3. Fase de Validación Taxonómica (Entorno: NCBI BLASTn)**  
La caracterización y validación taxonómica de los scaffolds obtenidos se realizó mediante la herramienta BLASTn contra la base de datos de referencia de nucleótidos estándar (nt/nr) del NCBI. Este alineamiento global se ejecutó en dos etapas independientes: primero, para identificar la naturaleza biológica de los bloques genómicos preliminares y confirmar la presencia del hospedero bacteriano, y segundo, para certificar la pureza, el porcentaje de identidad molecular y el linaje taxonómico oficial del genoma viral aislado tras la depuración.    

<img width="1303" height="780" alt="Captura de pantalla 2026-05-12 125110" src="https://github.com/user-attachments/assets/9dde7435-f898-4554-8189-3b72c28b2b7c" />

*Fig. 8*  Validación taxonómica en BLASTn 


## 3. RESULTADOS:     
### 3.1 Control de calidad:  

Se evaluó la calidad de las lecturas crudas mediante FastQC, observando la necesidad de un proceso de limpieza debido a la presencia de adaptadores. Tras aplicar el filtrado con Trimmomatic, se obtuvo un reporte final con un 97.25% de bases con calidad superior a Q30, garantizando datos confiables para el ensamblaje

<img width="935" height="699" alt="image" src="https://github.com/user-attachments/assets/17636f6d-5e3e-40f4-88ed-cbbb6c01ca38" />

Fig. 9 FastQC: Per Base Sequence Content proporción de cada una de las cuatro bases nitrogenadas (Timina %T, Citosina %C, Adenina %A y Guanina %G) en cada posición a lo largo de las lecturas de secuenciación.

Se realizó una evaluación exhaustiva de la integridad y composición de los datos de secuenciación masiva correspondientes a las lecturas directas (forward) e inversas (reverse) de la muestra DRR817419, comparando su estado crudo inicial con el obtenido tras la curación bioinformática.

La composición nucleotídica global se evaluó a través del módulo Per Sequence Content, analizado tanto en su distribución porcentual como en su recuento absoluto de lecturas. Ambas métricas revelaron un perfil unimodal y simétrico compatible con una distribución normal teórica, posicionando el pico de abundancia máxima en un contenido medio de GC de aproximadamente 50,2%. La alerta preventiva (amarilla) emitida inicialmente por el programa se atribuyó a una ligera asimetría en la cola izquierda de la curva (rango de 30% a 45% de GC), la cual responde a la heterogeneidad natural de las regiones transcritas ricas en bases AT.

El análisis de la composición posicional mediante el módulo Per Base Sequence Content exhibió originalmente una marcada fluctuación en las proporciones de las cuatro bases entre las posiciones 1 y 9 del extremo 5'. Este comportamiento es un sesgo técnico característico derivado del cebado aleatorio durante la preparación de la librería. Con el objetivo de mitigar esta distorsión y evitar penalizaciones o desalineamientos en las herramientas analíticas posteriores, se procedió a realizar un recorte adaptativo estricto de los primeros 15 nucleótidos en el extremo 5' de las lecturas.

Como resultado de este procesamiento, el módulo Per Base Sequence Content logró una validación exitosa (criterio de aprobación verde). Las curvas de abundancia para adenina, timina, citosina y guanina muestran una convergencia absoluta y paralela en torno al 25% cada una, manteniéndose con total estabilidad y linealidad a lo largo de toda la extensión remanente de los fragmentos. Paralelamente, el módulo Sequence Length Distribution reflejó esta modificación metodológica mediante una reducción proporcional en la longitud máxima de lectura, confirmando la remoción homogénea del bloque nucleotídico inicial sesgado.

<img width="975" height="488" alt="image" src="https://github.com/user-attachments/assets/e9300cf0-4f5c-457c-abff-ac226a2b1cfa" />

Fig. 10 Evaluación de calidad de las lecturas del bacteriófago T4 mediante FastQC y MultiQC.

Se evaluó la composición y el estado de los datos de secuenciación masiva correspondientes a las lecturas directas (forward) e inversas (reverse) de la muestra DRR817419, tanto en su estado crudo inicial como posterior al proceso de curación bioinformática.

El perfil de recuento de secuencias obtenido mediante la herramienta FastQC (integrado en MultiQC) reveló un volumen inicial aproximado de 4.8 millones de lecturas por cada archivo pareado (DRR817419_forward y DRR817419_reverse).

Tras la aplicación del software Trimmomatic para la eliminación de adaptadores y el filtrado de bases de baja calidad, se observó una reducción marginal en el número total de lecturas, estabilizándose en aproximadamente 4.7 millones de lecturas retenidas por archivo. Esta pérdida controlada valida la especificidad del proceso de limpieza, garantizando que no se descartó información biológica masiva de forma errónea, sino únicamente secuencias artefactuales o de calidad insuficiente.

Este nivel de duplicación es consistente y el rendimiento cuantitativo y la retención de datos tras el trimado confirman que las muestras procesadas poseen la integridad y el volumen necesarios para continuar con las etapas posteriores de ensamblaje o alineamiento contra referencia.

### 3.2 Ensamblaje genómico:  
#### 3.2.1 Evaluación Estructural del Ensamblaje Preliminar (QUAST)
El ensamblaje de *novo* a partir de las lecturas filtradas se ejecutó mediante el algoritmo SPAdes dentro de la plataforma Galaxy. Para evaluar la continuidad, fragmentación y éxito general de la reconstrucción molecular, se analizaron las métricas estadísticas estructurales obtenidas a través de la herramienta QUAST:

<img width="622" height="637" alt="image" src="https://github.com/user-attachments/assets/70f3f5ad-c09d-4dbe-8a0b-993a21149705" />

Fig. 11 Resultados del ensamblaje de *novo* del bacteriófago T4 obtenidos mediante SPAdes. Se observa un contig principal con elevada cobertura y longitud.  

Respecto al rendimiento de la reconstrucción, el proceso generó un total de 88 scaffolds totales, de los cuales únicamente 81 scaffolds presentaron una longitud útil mayor o igual a 1,000 pb ($\ge$ 1,000 pb). Esta relación indica que el algoritmo operó con una limpieza técnica sobresaliente, puesto que la inmensa mayoría de los scaffolds construidos superaron el umbral de las mil bases, reduciendo casi por completo la presencia de fragmentos menores dispersos o ruido bioinformático.

En su conjunto, la longitud acumulada de este ensamblaje alcanzó una extensión total de 4,638,873 pb (aproximadamente 4.64 Mb). Esta escala macroscópica de nucleótidos representa una sobredimensión crítica respecto al tamaño biológico esperado para el genoma de referencia del Bacteriófago T4 (~169 kb), evidenciando que el software reconstruyó una masa cromosómica casi 30 veces mayor a la del virus debido a una masiva co-secuenciación de material genético celular.

Esta hipótesis de contaminación por el hospedero se corrobora matemáticamente al analizar el contenido de Guanina y Citosina, el cual registró un promedio del 50.2% para el total de las secuencias moleculares obtenidas. Dado que el %GC teórico del Bacteriófago T4 es característicamente bajo (~34%) y el de *Escherichia coli* ronda el ~50%, este sesgo composicional actúa como una firma molecular irrefutable de que el ensamblaje preliminar está constituido primordialmente por el genoma de la bacteria hospedera, enmascarando las secuencias del virus.

Por otra parte, al evaluar la continuidad del ensamblaje mediante las métricas estadísticas N50 (118,604 pb) y L50 (12), se evidencia una alta estabilidad técnica en el proceso. En el contexto bioinformático, el valor L50 de 12 actúa como un indicador cuantitativo de orden que certifica que la mitad de la masa total de este gigantesco genoma (más de 2.3 Mb) se encuentra concentrada de forma eficiente en apenas 12 scaffolds de gran tamaño. Asimismo, el valor N50 complementa este criterio de calidad al establecer un umbral de longitud, certificando que el menor de los fragmentos dentro de este bloque principal mide 118,604 pb y asegurando que el algoritmo SPAdes no fragmentó en exceso la secuencia consenso. Dentro de esta distribución, destacó de forma individual la resolución del contig de máxima extensión, el cual alcanzó los 327,290 pb. En conclusión, este primer escrutinio estructural demuestra que el pipeline procesó y acopló con éxito los datos crudos, pero expone la necesidad estricta de ejecutar un paso posterior de discriminación molecular y filtrado taxonómico para aislar las lecturas virales de los bloques bacterianos predominantes.   

#### 3.3.1 Depuración Genómica mediante Mapeo de Lecturas (Bowtie 2)  
Debido a la masiva co-secuenciación del hospedero *Escherichia coli* evidenciada en el análisis de QUAST, se procedió a ejecutar una etapa de filtrado por exclusión en entorno de terminal para aislar las lecturas correspondientes al bacteriófago T4. Para optimizar el pipeline bioinformático, se tomaron las lecturas de alta calidad previamente procesadas por Trimmomatic/fastp y se mapearon directamente con la herramienta Bowtie 2 contra el genoma de referencia de *Escherichia coli*.    

La aplicación de este mapeo de alta sensibilidad arrojó una tasa de alineamiento específica del 0.60%, logrando capturar e identificar de forma exacta un total de 28,104 lecturas verdaderamente virales. Esta estrategia permitió discriminar con éxito el contenido del virus frente al predominante ruido molecular bacteriano sin necesidad de re-evaluar la calidad general de los datos.

Este set optimizado de 28,104 lecturas remanentes, correspondientes exclusivamente al virus, fue sometido directamente a un segundo proceso de ensamblaje de *novo* en SPAdes Terminal. Este filtrado resolvió con éxito un único scaffold unificado (NODE_1) de 168,129 pb, libre de contaminación biológica y listo para su correspondiente caracterización taxonómica.   

### 3.3 Caracterización y Validación Taxonómica del Genoma Viral Aislado:   
Una vez obtenido el scaffold definitivo NODE_1 de 168,129 pb mediante SPAdes Terminal, se procedió a realizar su caracterización biológica y validación taxonómica. Para comprobar la veracidad estructural del genoma viral reconstrucido y descartar cualquier residuo del hospedero, se ejecutó un alineamiento nucleotídico local mediante la herramienta BLASTn contra la secuencia de referencia oficial del bacteriófago T4 (NC_000866.4) depositada en la base de datos del NCBI.  

Los parámetros métricos oficiales obtenidos en este análisis global se presentan consolidados en la siguiente tabla:  

| Parámetro Métrico | Valor Obtenido | Herramienta / Plataforma |
| :--- | :--- | :--- |
| **Identidad Taxonómica (Hit Principal)** | **99.98%** (*Escherichia virus T4*) | NCBI BLASTn |
| **Cobertura de Consulta (Query Cover)** | **100%** | NCBI BLASTn |
| **Longitud del Scaffold Viral (NODE_1)** | **168,129 pb** | SPAdes Terminal |
| **Profundidad de Cobertura Genómica** | **24.64x** | SPAdes Terminal |
| **E-value Estadístico** | **0.0** | NCBI BLASTn |
| **Bases de Calidad Obtenidas (Q30)** | **97.25%** | fastp / Trimmomatic |  

*Tabla. 1* Parámetros oficiales obtenidos.    

El análisis arrojó un valor numérico esperado ($E\text{-value}$) de 0.0 y una cobertura de consulta (Query Cover) del 100%, confirmando una coincidencia molecular exacta. Asimismo, se registró un porcentaje de identidad del 99.98% con la secuencia completa de *Escherichia* virus T4. La diferencia marginal de apenas ~77 pb respecto al genoma de referencia internacional (168,903 pb) evidencia la alta fidelidad del pipeline bioinformático y la robustez del algoritmo de ensamblaje por grafos de De Bruijn a partir de lecturas cortas pareadas (paired-end).  

Es decir, el re-ensamblaje enfocado únicamente en estas lecturas purificadas resolvió por completo el ruido del hospedero, generando un scaffold definitivo excepcional.  

<img width="1303" height="780" alt="Captura de pantalla 2026-05-12 125110" src="https://github.com/user-attachments/assets/9dde7435-f898-4554-8189-3b72c28b2b7c" />

<img width="1308" height="747" alt="Captura de pantalla 2026-05-12 125615" src="https://github.com/user-attachments/assets/1298172e-3da6-49b7-aeef-c8ea6345bce0" />  

*Fig. 12* Resultado del alineamiento local BLASTn del scaffold preliminar obtenido en Galaxy, evidenciando una homología nucleotídica del 99.88% con el genoma cromosómico de *Escherichia coli strain C*.

<img width="1302" height="792" alt="Captura de pantalla 2026-05-14 134801" src="https://github.com/user-attachments/assets/a7f878cb-e289-4f51-b286-f84bd087f532" />

<img width="1321" height="804" alt="Captura de pantalla 2026-05-14 134945" src="https://github.com/user-attachments/assets/996a315b-1516-466b-ae35-d04efb00a3d0" />   

*Fig. 13* Resultado del alineamiento BLASTn del scaffold definitivo NODE_1 (168,129 pb) tras la depuración con Bowtie 2, certificando una identidad molecular del 99.99% con el genoma de referencia de *Tequatrovirus* T4.  


### 3.5 Clasificación Taxonómica Oficial:   

Esta categorización sistemática concluye la fase analítica del proyecto, certificando que el fago ensamblado corresponde con absoluta pureza al organismo objetivo del estudio y proporcionando una secuencia consenso de alta resolución molecular para posteriores análisis funcionales. El perfil de clasificación molecular adscribe el scaffold definitivo **168,129 pb** (*NODE_1*) de forma inequívoca dentro de la jerarquía taxonómica oficial aprobada por el Comité Internacional de Taxonomía de Virus (ICTV):

| Nivel Taxonómico | Clasificación Científica |
| :--- | :--- |
| **Dominio Viral** | *Viruses* |
| **Clase** | *Caudoviricetes* |
| **Género** | *Tequatrovirus* |
| **Especie** | *Escherichia virus T4* (antes *Escherichia phage* T4) |

La detección robusta de linajes específicos como *Tequatrovirus* T4 y *Escherichia virus T4* con un **99.98% de identidad** valida con éxito el flujo de trabajo implementado. 

Desde una perspectiva biotecnológica y microbiológica, la confirmación de esta identidad es un pilar fundamental. El fago T4 es un sistema modelo ampliamente estudiado en la biología molecular y la genómica viral debido a su estricto ciclo lítico. Los datos genómicos limpios obtenidos en este proyecto respaldan su viabilidad y seguridad como un candidato biológico óptimo para el desarrollo de terapias fágicas avanzadas y el control epidemiológico de cepas multirresistentes de *Escherichia coli*.  

## 4. DISCUSIÓN:  

Los resultados obtenidos en el presente estudio demuestran que el flujo de trabajo bioinformático implementado permitió realizar exitosamente el ensamblaje de novo y la validación taxonómica del bacteriófago T4. La adecuada calidad de las lecturas observada mediante FastQC y MultiQC evidenció que el dataset poseía características óptimas para análisis genómicos posteriores, lo cual coincide con lo reportado por (Hernández et al.,2020), quienes destacan que la calidad inicial de las secuencias es un factor determinante para garantizar resultados confiables en estudios de secuenciación de nueva generación.

Análisis de los hallazgos:

Relación Huésped-Parásito: Los bacteriófagos son parásitos obligados que requieren la maquinaria celular de una bacteria para su replicación. En el dataset DRR817419, la presencia masiva de secuencias bacterianas es técnicamente esperada, ya que el ADN del fago se extrae a menudo de cultivos infectados donde el ADN de la bacteria anfitriona (*E. coli*) coexiste en mayor proporción genómica.

Calidad del Ensamblaje: El hecho de haber obtenido scaffolds de gran longitud (superior a 320 kb) con una identidad del 99.90% demuestra que el preprocesamiento con Trimmomatic/fastp y el ensamblaje con SPAdes fueron altamente eficientes. Un ensamblaje pobre habría generado miles de fragmentos pequeños y baja identidad, lo cual no ocurrió en este caso.

Estrategia de Filtrado: Este resultado marca la finalización exitosa de la Fase 1 del proyecto. La detección del genoma de *E. coli* permite ahora aplicar técnicas como de depuración bioinformática, como el mapeo de lecturas contra un genoma de referencia del fago T4, para aislar exclusivamente las secuencias virales de interés terapéutico.

El ensamblaje de novo permitió reconstruir contigs de alta cobertura, siendo el principal contig identificado de 327,290 pb con una cobertura aproximada de 176.37X. Estos resultados reflejan una elevada representación del material genético viral dentro de la muestra analizada y respaldan la estabilidad del ensamblaje obtenido. De acuerdo con (Basantani et al., 2017), las herramientas bioinformáticas modernas permiten reconstrucciones genómicas altamente eficientes, facilitando el análisis funcional y taxonómico de organismos complejos, incluidos bacteriófagos.

La validación taxonómica realizada mediante Kraken2 confirmó la presencia de *Escherichia coli* T4 y *Tequatrovirus* T4, validando taxonómicamente el ensamblaje generado. Además, el alto porcentaje de lecturas clasificadas (98.07 %) indica una adecuada correspondencia entre las secuencias ensambladas y las bases de datos de referencia utilizadas. (Hernández et al., 2020) señalan que las herramientas de clasificación taxonómica basadas en secuenciación masiva permiten identificar organismos con elevada precisión, especialmente cuando se dispone de datasets de buena calidad.

La abundante presencia de secuencias asociadas a *Escherichia coli* observada durante la clasificación taxonómica es consistente con la biología del bacteriófago T4, ya que este utiliza cepas de *E. coli* como hospedero natural. Estudios recientes han demostrado que el proceso de infección del bacteriófago T4 involucra complejas interacciones moleculares con proteínas bacterianas específicas, incluyendo proteínas relacionadas con la división celular bacteriana (Wenzel et al., 2024).
Asimismo, (Wolfram-Schauerte et al., 2022) describen que la infección de *E. coli* por el bacteriófago T4 genera importantes cambios transcriptómicos y proteómicos en la bacteria hospedera, lo que evidencia el elevado grado de adaptación biológica entre ambos organismos. Esto explica la fuerte asociación taxonómica detectada entre el bacteriófago y su hospedero durante el análisis bioinformático.

Los bacteriófagos han adquirido creciente importancia como alternativas terapéuticas frente a bacterias resistentes a antibióticos. (Kortright et al., 2019) indican que los fagos representan herramientas promisorias en aplicaciones biomédicas debido a su especificidad y capacidad de eliminar bacterias patógenas resistentes. En este contexto, la correcta identificación del bacteriófago T4 obtenida en el presente estudio respalda su posible aplicación biotecnológica para el control de cepas resistentes de *Escherichia coli*.

De igual manera, (Salmond y Fineran., 2015) destacan que los bacteriófagos constituyen uno de los sistemas biológicos más relevantes en microbiología moderna, tanto por su valor en investigación molecular como por su potencial en terapias antimicrobianas. La información genómica obtenida mediante ensamblaje y validación taxonómica constituye una base fundamental para futuras investigaciones funcionales y aplicaciones terapéuticas.

Por otra parte, (Liu et al., 2026) reportan que algunas bacterias poseen sistemas de defensa dirigidos específicamente contra modificaciones del ADN genómico de fagos, demostrando la complejidad evolutiva de la interacción bacteria-fago. Estos mecanismos resaltan la importancia de caracterizar adecuadamente los genomas virales mediante herramientas bioinformáticas antes de considerar aplicaciones terapéuticas.

Finalmente, (Clokie y Kropinski., 2009) señalan que la caracterización molecular y genómica de bacteriófagos constituye un paso esencial para comprender su biología, diversidad y potencial aplicación en microbiología clínica y ambiental. En este sentido, los resultados obtenidos en el presente trabajo demuestran que el ensamblaje de novo y la validación taxonómica son estrategias efectivas para la identificación y análisis de bacteriófagos con potencial biotecnológico.



## 5. CONCLUSIÓN    

* El análisis estadístico estructural con QUAST demostró que el ensamblaje preliminar en Galaxy arrojó una marcada sobredimensión genómica (4.64 Mb) y un sesgo composicional de %GC alto (50.2%), actuando como una firma molecular inequívoca de una co-secuenciación masiva del hospedero bacteriano *Escherichia coli*.

* La implementación del pipeline de depuración en terminal mediante Bowtie 2 demostró ser una estrategia de alta sensibilidad y especificidad, permitiendo segregar con éxito el ruido molecular de la bacteria para aislar un set optimizado de 28,104 lecturas verdaderamente virales (tasa de alineamiento del 0.60%).

* El re-ensamblaje con SPAdes Terminal resolvió de forma exitosa un único scaffold unificado (NODE_1) de 168,129 pb, cuya validación global en NCBI BLASTn ratificó una identidad nucleotídica del 99.98% con el genoma de referencia de *Escherichia* virus T4, certificando la máxima pureza biológica de la secuencia consenso obtenida.

* Los resultados obtenidos reafirman el potencial del bacteriófago T4 como una alternativa biológica viable para el control de cepas resistentes de *Escherichia coli*, consolidando un flujo de trabajo bioinformático mixto (Galaxy-Terminal) altamente eficiente para la caracterización genómica de especímenes virales.


## 6. REFERENCIAS BIBLIOGRÁFICAS:  

*   Clokie, M. R., & Kropinski, A. M. (Eds.). (2009). *Bacteriophages: Methods and Protocols*. Humana Press.
*   Kortright, K. E., Rozo, S. D., Lohrmeyer, A. J., & Turner, P. E. (2019). The Epidemiology of Bacteriophages and Their Biomedical Applications. *Cell Host & Microbe, 25*(2), 219-232. https://doi.org/10.1016/j.chom.2019.01.005
*   Liu, R., Tang, D., Niu, M., Lei, S., Zong, Z., Chen, Q., & Yu, Y. (2026). A bacterial defense system targeting modified cytosine of phage genomic DNA. *Nature Communications*, *17*(1920), 1–10. doi.org
*   Salmond, G. P., & Fineran, P. C. (2015). A century of the phage: Past, present and future. *Nature Reviews Microbiology, 13*(12), 777-786. https://doi.org/10.1038/nrmicro3564
*   Wenzel, S., Hess, R., Kiefer, D., & Kuhn, A. (2024). Involvement of the Cell Division Protein DamX in the Infection Process of Bacteriophage T4. *Viruses*, *16*(4), 487. doi.org
*   Wolfram-Schauerte, M., Pozhydaieva, N., Viering, M., Glatter, T., & Höfer, K. (2022). Integrated Omics Reveal Time-Resolved Insights into T4 Phage Infection of E. coli on Proteome and Transcriptome Levels. *Viruses*, *14*(11), 2502. doi.org
*   Basantani, M. K., Gupta, D., Mehrotra, R., Mehrotra, S., Vaish, S., & Singh, A. (2017). An update on bioinformatics resources for plant genomics research. In Current Plant Biology (Vols. 11–12, pp. 33–40). *Elsevier B.V.* https://doi.org/10.1016/j.cpb.2017.12.002
*   Hernández, M., Quijada, N. M., Rodríguez-Lázaro, D., & Eiros, J. M. (2020). Bioinformatics of next generation sequencing in clinical microbiology diagnosis. *Revista Argentina de Microbiologia*, 52(2), 150–161. https://doi.org/10.1016/j.ram.2019.06.003






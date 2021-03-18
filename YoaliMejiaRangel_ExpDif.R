##########REPOSITORIO GITHUB PARA EXPRESIÓN DIFERENCIAL##########
##YOALI MEJIA RANGEL

library("sleuth")
library("shiny")
library("ensembldb")
#Esta función permite mapear, a partir de la base de datos de homo sapiens 
tx2gene <- function(){
  
  mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  mart
  t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                       "external_gene_name"), mart = mart)
  t2g <- dplyr::rename(t2g, target_id = "ensembl_transcript_id",
                       ens_gene = "ensembl_gene_id", ext_gene = "external_gene_name")
  return(t2g)
}

t2g <- tx2gene() #esta función será utilizada para los análisis posteriores. Se le asigna la función al vector t2g


base_dir<-"C:/Users/HP/Documents/6to semestre/Genómica/samples"
base_dir #se debe primero localizar el lugar donde están mis archivos, en este caso están dentro d ela carpeta "samples"

samples <- paste0("sample", c("1","2","3","10","11","12")) #se coloca el nombre de la carpeta y los archivos a tomar en cuenta 
samples #en este caso, la carpeta es sample y se tomaron los primeros y últimos 3 archivos (1:3 y 10:12) 

kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
kal_dirs #asignar y ver las carpetas que se utilizaran "sample1, sample2, sample3... con el lugar de origen de base_dir

s2c <- data.frame(path=kal_dirs, sample=samples, muestras=c("control","control","control", "muted",
                                                            "muted", "muted"), stringsAsFactors=FALSE)
s2c #se asigna el control y mutantes para que sean comparados. Los controles son los primeros 3 archivos y las mutantes los últimos 3
so <- sleuth_prep(s2c, ~muestras, target_mapping = t2g,extra_bootstrap_summary = TRUE) #este vector ayuda a normalizar los datos
#se lleva a cabo el bootstrap , normalización de los datos . Datos son leídos en kallisto 
so <- sleuth_fit(so)
so <- sleuth_wt(so, which_beta="muestrasmuted") #condición de referencia, comparadas con las mutantes                  
sleuth_live(so)
#library (shiny) #se necesita cargar a la librería shiny para realizar el siguiente análisis
#se abre una nueva ventana arrojando los resultados, heatmaps, volcano plots, tablas, etc. 

setwd("6to semestre/Genómica/samples/") #nuevamente se cambia de directorio para nuestra carpeta con los archivos 

#una vez que se abre la pestaña con los resultaods (para esto la librería shiny), se descarga la tabla de test table
#y se carga en R para realizar los análisis siguientes : 
resultados<-read.table("test_table.csv",sep=",",
                       header=TRUE) #aquí se carga el archivo descargado 
significativos<-which(resultados$qval<0.1) #separa a los resultados con base en un q valor<0.1
significativos<-resultados[significativos,] #muestra los genes que son significativos
upregulated<-which(significativos$b>0) #muestra los genes más regulados 
upregulated<-significativos[upregulated,]
downregulated<-which(significativos$b<0) #muestra los genes menos regulados con base en los análisis 
downregulated<-significativos[downregulated,]

write.table(upregulated,file="~/6to semestre/Genómica/samples/Upregulated_N2vsCyg-25.txt",sep="\t")
write.table(downregulated,file="~/6to semestre/Genómica/samples/Downregulated_N2vsCyg-25.txt",sep="\t")

upregulated$ext_gene

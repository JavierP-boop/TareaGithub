


# Se cargan las librerías que se vayan a utilizar

library(sleuth)
library (devtools)
library(ensembldb)
library(edgeR)
library(rhdf5)


# Con esta funcion se mapea a traves de la base de datos

tx2gene <- function(){
  
  
    #     Dataset you want to use. To see the different datasets available within a biomaRt yo$
    #     host
    
    mart <- biomaRt::useMart(biomart = "ensembl", dataset = "celegans_gene_ensembl")
  
  t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                       "external_gene_name"), mart = mart)
  t2g <- dplyr::rename(t2g, target_id = "ensembl_transcript_id",
                       ens_gene = "ensembl_gene_id", ext_gene = "external_gene_name")
  return(t2g)
}



# Con esta funcion podemos acceder a las bases de datos que tiene Ensembl y 
# ya con esto podemos usar sus transcritos

t2g <- tx2gene()


# Se asigna un objeto a la carpeta donde se extrajeron los samples 

base_dir <- "TareaGitK/"


# Se seleccionan las muestras 
# Asiganmos a un objeto un numero de muestras (no muchas porque se puede trabar)

samples <- paste0("sample", c("1","2","3","4","5"))



# Se genera una ruta de acceso a las muestras

kal_dirs <- sapply(samples, function(id) file.path (base_dir, id))



#Selected samples

s2c <- data.frame(path=kal_dirs, sample=samples, muestras = c("sample6","sample6",
                                                              "sample7", "sample7"), stringsAsFactors=FALSE)
# Se definen condiciones experimentales
# el codigo da error si solo se maneja una condicion, por ende se repiten (sample6 y sample7)
# Con esto generamos la data frame

so <- sleuth_prep(s2c, ~ muestras, target mapping = t2g, extra_bootstrap_summary = TRUE)
# no estoy seguro que hace esto


so <- sleuth_fit(so)
# Se ajustan modelos de error, varianza biologíca 

so <- sleuth_wt(so, which_beta="muestrassample6")                   
# Se realiza el test de Wald

sleuth_live(so)

# Nos muestra varias cosas, entre ellas un volcano plot, histogramas, graficas de componentes principales
# Muestra varias distribuciones y una tabla de kalisto



resultados<-read.table("TareaGitK/test_table.csv",sep=",", header=TRUE)

# Se asigna a un objeto los resultados del analisis



significativos<-which(resultados$qval<0.1)    # Nos da valores significativos
significativos<-resultados[significativos,]  # Ordena los valores significativos 
upregulated<-which(significativos$b>0)       # Nos da los resultados sobreregulados (mayor a 0)
upregulated<-significativos[upregulated,]    # Ordena los valores sobreregulados
downregulated<-which(significativos$b<0)    # Nos da los resultados que estan subregulados  (menor a 0)
downregulated<-significativos[downregulated,]  # Ordena los resultados que estan subregulados





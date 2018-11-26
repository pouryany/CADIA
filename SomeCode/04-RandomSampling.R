## Random testing evaluation


library(CADIA)
library(dplyr)
library(ggplot2)


# Loading the list of the ENTREZ ID from a HG-U133
tT.all.names <- readRDS("data-raw/allGeneNames")


err.mat.ora   <- matrix(NA,143,500)
err.mat.cdist <- matrix(NA,143,500)
err.mat.ssc   <- matrix(NA,143,500)

for(i in 1:50){
    
    for(j in 1:10){
        
        k <- ((i-1)*10) + j
        
        set.seed(k)
        tT.de.names   <- sample(tT.all.names,i*100)
        tT.pathways   <- causalDisturbance(tT.de.names,tT.all.names,iter = 2000,
                                           alpha = 0.1,statEval = 1)
        
        
        err.mat.ora[,k]   <- tT.pathways$P_ORA
        err.mat.ssc[,k]   <- tT.pathways$P_SSC
        err.mat.cdist[,k] <- tT.pathways$cadia
        
        print(paste0(k,"\n"))
    }
    
}



sum((!(as.vector(err.mat.ssc) != 1) == (as.vector(err.mat.ora) != 1)))


ggplot( data = data.frame(x = as.vector(err.mat.ora)), aes(x=x)) + 
        geom_histogram(color="black", fill="white",binwidth=0.02) +
        geom_hline(yintercept=length(as.vector(err.mat.ora)) /50,
                   size = 1, linetype="dashed", color = "red")+
        theme_bw()+
        labs(y = "", 
             x = "ORA P-values") +
        theme(legend.position = c(0.9, 0.5),
              axis.title = element_text(size = 20),
              legend.text = element_text(size = 15),
              axis.text=element_text(size= 20),
              axis.text.y= element_text(size = 20),
              legend.title=element_blank()) 


ggplot( data = data.frame(x = as.vector(err.mat.ssc)), aes(x=x))  + 
        geom_histogram(color="black", fill="white",binwidth=0.02) +
        geom_hline(yintercept=length(as.vector(err.mat.ora)) /50,
                   size = 1, linetype="dashed", color = "red") +
        theme_bw()+
        labs(y = "", 
             x = "SSC P-values") +
        theme(legend.position = c(0.9, 0.5),
              axis.title = element_text(size = 20),
              legend.text = element_text(size = 15),
              axis.text=element_text(size= 20),
              axis.text.y= element_text(size = 20),
              legend.title=element_blank()) 


# Correlation analysis, for the pathways for which 
# SSC and ORA p-values are non-1 (calculated)

not.one     <- (as.vector(err.mat.ora) != 1 & as.vector(err.mat.ssc) != 1)
pvals_table <- data.frame(ORA = as.vector(err.mat.ora)[not.one],
               SSC =as.vector(err.mat.ssc)[not.one])


ggplot(data = pvals_table, aes(x= ORA, y=SSC)) + geom_point(size=0.1) +
       theme_bw()+
       labs(y = "SSC P-values", 
             x = "ORA P-values") +
       theme(legend.position = c(0.9, 0.5),
             axis.title = element_text(size = 20),
             legend.text = element_text(size = 15),
             axis.text=element_text(size= 20),
             axis.text.y= element_text(size = 20),
             legend.title=element_blank()) 


cor.test(as.vector(err.mat.ora)[not.one],as.vector(err.mat.ssc)[not.one])






err.fdr.ora <- apply(err.mat.ora, 2, p.adjust)
err.fdr.ora <- colSums(err.fdr.ora < 0.05) /143

err.fdr.ssc <- apply(err.mat.ssc, 2, p.adjust)
err.fdr.ssc <- colSums(err.fdr.ssc < 0.05) /143


err.fdr.cdist <- colSums(err.mat.cdist < 0.05) / 143


zz1<- vector()
zz2<- vector()
zz3<- vector()
for (j in 1:50){
    zz1[j] <- mean(err.fdr.ora[(10*(j-1)+1):(10*j)])
    zz2[j] <- mean(err.fdr.ssc[(10*(j-1)+1):(10*j)])
    zz3[j] <- mean(err.fdr.cdist[(10*(j-1)+1):(10*j)])
}


plot_table <- data.frame(sample_size = (1:50), CADIA = zz3 , ORA = zz1,
                         SSC = zz2)

plot_table <- melt(plot_table,id.vars =  c(1))


ggplot(data = plot_table, aes(x = sample_size, y = value, fill = variable)) + 
    geom_bar(stat = 'identity', position = 'dodge') +  
    geom_hline(yintercept=0.05, size = 1, linetype="dashed", color = "red")+
    labs(y = "Average False Positive Rate", 
         x = "Random Perturbation Sample Size (X100)") +
    theme_bw() +
    theme(legend.position = c(0.9, 0.5),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        axis.text=element_text(size= 20),
        axis.text.y= element_text(size = 20),
        legend.title=element_blank()) 




saveRDS(err.mat.cdist, file = "data-raw/err_cdist")
saveRDS(err.mat.ora, file = "data-raw/err_ora")
saveRDS(err.mat.ssc, file = "data-raw/err_ssc")



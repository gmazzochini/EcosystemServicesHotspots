#rst = raster of the ecosystem service or biodiversity
#hsquant = quantile to be measured inside protected areas
#n.ite = number of interactions
#uc.shp = shapefile of protected areas
#tab.name = name of the organized table to be used on the analyses
#table.ready = Logical. If table has been already created table.ready = T

hotspotsUC <- function(rst, hsquant, n.ite, uc.shp, tab.name, table.ready = F){
  ### MONTAR TABELA DE  DADOS #####
  ID_UC <- uc.shp@data$ID_UC
  
  if(table.ready==FALSE){
  ## Pegar coords dos pixels
  tab <- data.frame(round(xyFromCell(rst,which(!is.na(rst[]))),4),
                  hot_spot=as.numeric(na.omit(getValues(rst>quantile(rst,hsquant)))),
                  uc=rep(0,length(na.omit(rst[])))) ## criar tab de dados
  
    for (i in unique(ID_UC)){
      uc <- uc.shp[uc.shp$ID_UC==i,]
      uc_rst  <-  rasterize(uc, rst,1) 
      uc.xy <- data.frame(round(xyFromCell(uc_rst,which(uc_rst[]>0)),4))
      tab[interaction(tab[,1:2])%in%interaction(uc.xy),"uc"] <- i
    }
  write.table(tab,file=paste(tab.name),row.names=F)
  cat("Table is ready and saved.\n")
  }
  
  if(table.ready==TRUE){
    tab <- read.table(tab.name,h=T)
  }
  
  Rmat <- numeric(n.ite)
  for (j in 1:n.ite){
    uc.aleat.tab <- data.frame(tab[,1:3],uc.aleat=rep(0,nrow(tab)))
    for (i in unique(ID_UC)){
      xy.uc <- tab[tab[,4]==i,1:2]
      OK=TRUE
      x.dif <- xy.uc[,1]-xy.uc[1,1]
      y.dif <- xy.uc[,2]-xy.uc[1,2]
      while(OK==TRUE){ ### verificar se a nova posicao da UC esta fora da caatinga
        xy.aleat <- tab[sample(1:nrow(tab),1),1:2]
        uc.aleat <- round(data.frame(x=xy.aleat[,1]+x.dif,y=xy.aleat[,2]+y.dif),4)
        OK=any(is.na(rst[cellFromXY(rst,uc.aleat[,1:2])])) ### se OK= T, ent?o,continua while
      }
      uc.aleat.tab[interaction(uc.aleat.tab[,1:2])%in%interaction(uc.aleat[,1:2]),"uc.aleat"] <- 1
    }
    Rmat[j] <- (sum(c(uc.aleat.tab$hot_spot+as.numeric(uc.aleat.tab$uc>0)==2))/sum(uc.aleat.tab$hot_spot))*100
    
    if(j==n.ite){cat(paste(j,". Finished!\n"),sep="")}
    else{cat(paste(j," ",sep=""))}
  }

  obs <- (sum(c(tab$hot_spot+as.numeric(tab$uc>0)==2))/sum(tab$hot_spot))*100
  p <- sum(Rmat>obs)/n.ite

  results <- data.frame(Obs=round(obs,2),P=p)
  windows(200,100)
  par(mfrow=c(1,2))
  if(results$P<0.05)a <- hist(round(Rmat,2),xlab="Obs",xlim=c(0,obs+obs/2))
  if(results$P>=0.05)a <- hist(round(Rmat,2),xlab="Obs") 
  arrows(x0=obs,y1=0,y0=max(a$counts)/2,col="red",lwd=3)
  plot(rst>quantile(rst,hsquant))
  plot(uc.shp,add=T,lwd=2)
  cat("\n")
  print.table(results)
  return(list(Obs=results[,1],P=results[,2],random.values=Rmat))
}

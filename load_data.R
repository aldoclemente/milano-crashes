
load_data <-function(){
  
  incidenti = read.csv(file = "IncidentiStradaliLombardia2018/Database Incidenti stradali con lesioni alle persone - Lombardia 2018.csv",
                       sep = ";", dec=",")
  
  COD_provincia = 15  # (MI)
  COD_comune    = 146 # Milano
  
  idxs = which(incidenti$COD_provincia == COD_provincia & 
                 incidenti$COD_comune    == COD_comune )
  
  
  incidenti_MI = incidenti[idxs,]
  
  idxs = ( is.na(incidenti_MI$latitudine) | is.na(incidenti_MI$longitudine))
  incidenti_MI = incidenti_MI[!idxs,]
  
  # elimino outliers per semplicità
  long_boxplot = boxplot(incidenti_MI$longitudine, plot=FALSE)
  long_idxs = match(long_boxplot$out, incidenti_MI$longitudine)
  
  lat_boxplot = boxplot(incidenti_MI$latitudine, plot=FALSE)
  lat_idxs = match(lat_boxplot$out, incidenti_MI$latitudine)
  
  idxs = unique( c(long_idxs, lat_idxs) )
  
  incidenti_MI = incidenti_MI[-idxs,]
  
  return(incidenti_MI)
}



library(oce)
library(sea)
library(xtable)


rootfold <- "~/data/SEA/jp/jpdata_edit"
cruiseIDs <- c("C187B","C193A","C199A","C205G","C211A","C218A","C223A","C230A","C235A","C241A","C248B")


ran <- c(39.5,41,-71.3,-70.6)


datefirst <- datelast <- nstat <- NULL
for (ID in cruiseIDs) {
  
  fold <- file.path(rootfold,ID,"cnv")
  foldcsv <- file.path(rootfold,ID,"csv")
  plotfold2 <- file.path(rootfold,ID)
  files <- list.files(path=fold, pattern="*.cnv")
  
    if(ID ==  cruiseIDs[8]){
      X<-readLatLon(filein=file.path(fold,files[2]))
    } else {
      X<-readLatLon(filein=file.path(fold,files[1]))
    }
    datel <- X$r[grep("date",X$r,ignore.case=T)[1]]
    datefirst <- append(datefirst,paste(tail(strsplit(datel,"[, ]+|[ ]+|Date|[-]|[*]")[[1]],3),collapse = " "))
    
    X<-readLatLon(filein=file.path(fold,files[length(files)]))
    datel <- X$r[grep("date",X$r,ignore.case=T)[1]]
    datelast <- append(datelast,paste(tail(strsplit(datel,"[, ]+|[ ]+|Date|[-]|[*]")[[1]],3),collapse = " "))
    
    nstat <- append(nstat,length(files))
}

df1 = strptime(datefirst,format="%d %B %Y",tz="EST")
df2 = strptime(datefirst,format="%d %B %y",tz="EST")
yrf = as.double(format(df1,"%Y"))
df1[yrf<1000] <- df2[yrf<1000]

dl1 = strptime(datelast,format="%d %B %Y",tz="EST")
dl2 = strptime(datelast,format="%d %B %y",tz="EST")
yrl = as.double(format(dl1,"%Y"))
dl1[yrl<1000] <- dl2[yrl<1000]

year = format(df1,"%Y")
date1 = format(df1,"%d-%B")
date2 = format(dl1,"%d-%B")

nw = as.integer(c(7,9,7,7,8,8,7,0,5,9,8))
ne = as.integer(c(8,10,7,5,6,7,7,7,7,9,7))

comments = rep("",length(cruiseIDs))
comments[1] = "Two casts in western section not used"
comments[8] = "NB: Cruise in late-July"

dfo <- data.frame(cruiseID = cruiseIDs,
                 year = year,
                 date1 = date1,
                 date2 = date2,
                 nstat = nstat,
                 nw = nw,
                 ne = ne,
                 comments = comments,
                 stringsAsFactors = F)

# caption <- "Hello"

# colnames(dfo) <- c('Cruise ID','Year','Start Date','End Date','','Number or Stations','','Comments')
firstRow <- c('Cruise ID','Year','Start Date','End Date','\\multicolumn{3}{c}{Number of Stations}','Comments')
secondRow <-c('','','','','Total','West','East','')

fileout <- "~/Documents/SEA/jp/overviewPaper/tables/cruiseOverview.tex"
# printTable(dfo,caption=caption,SigF=NULL,secondRow=secondRow,fileout=fileout)

r <- NULL
# r <- append(r,"\\begin{table}[!ht]")
# r <- append(r,paste0("\\caption{",caption,"}"))
# r <- append(r,"\\centering")
# r <- append(r,paste0("\\begin{tabular}{",paste(rep('c',ncol(dfo)),collapse = " "),"}"))
# r <- append(r,"\\hline")
r <- append(r,paste0(paste(firstRow,collapse=" & "),"\\\\"))
r <- append(r,paste0(paste(secondRow,collapse=" & "),"\\\\"))
r <- append(r,"\\hline")

for (i in 1:nrow(dfo)) {
  r <- append(r,paste0(paste(dfo[i,],collapse= " & "),"\\\\"))
}

# r <- append(r,"\\hline")
# r <- append(r,"\\end{tabular}")
# r <- append(r,"\\end{table}")
#   
fileConn<-file(fileout)
writeLines(r, fileConn)
close(fileConn)
           

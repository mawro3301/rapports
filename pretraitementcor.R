#library(rgeos)
#library(sp)
#library(rgdal)
library(readr)
library(base)
library(stats)
library(dplyr)
library(Rcpp)

setwd("~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/DATA/")

DatacarpoEV06_correction <- read.csv("DatacarpoEV06_correction.CSV",header="F")
View(DatacarpoEV06_correction)
donneesGeno= DatacarpoEV06_correction
dim(donneesGeno)
DatacarpoEV07 <- read.csv("~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/
                          DATA/DatacarpoEV07.csv")
View(DatacarpoEV07)
donneesGeno07=DatacarpoEV07



# donnees genetiques et parcelle d'origine
dP06Brut<-donneesGeno[,c(7,19:ncol(donneesGeno))]
dP07Brut<-donneesGeno07[,c(4,11:ncol(donneesGeno07))]
dP06Brut
dP07Brut
View(dP06Brut)
View(dP07Brut)

#fonction pour retirer les erreurs de génotypages
#cad les alles ac des zeros

zerosGoHome<-function(donnees) {
  listerLigne=0
  for (ligne in 1:nrow(donnees)) {
    for (colonne in 2:ncol(donnees)) {
      if (donnees[ligne,colonne]==0) {
        listerLigne<-cbind(listerLigne,ligne)
        break
      }
    }
  }
  print (listerLigne)
  k=0
  for (i in 2:ncol(listerLigne)) {
    donnees<-donnees[-c(listerLigne[i]-k),]
    k=k+1
  }
  return (donnees)
}

#on l'applique sur nos donnees
dP07T=zerosGoHome(dP07Brut)
dP06T=zerosGoHome(dP06Brut)
dP07T
dP06T
dP07Td<-data.frame(dP07T)
dP06Td<-data.frame(dP06T)

dPTotalT<-bind_rows(dP06T,dP07T)
View(dP07T)
View(dP06T)
View(dPTotalT)
dPTotalMoyen=dPTotalT
dPTotalMoyen[,1]=1000


#pour les freq alléliques initial, on prends les données 2006,
#augmenté d'un verger moyen (notre verger virtuel)
dPTotal2006WithMedian=bind_rows(dP06T,dPTotalMoyen)

#ListeVerger=unique(as.vector(as.matrix(dPTotalT[,1])))
ListeVerger=unique(as.vector(as.matrix(dPTotal2006WithMedian[,1])))
ListeCpNa=unique(as.vector(as.matrix(dPTotal2006WithMedian[,c(1+1,2+1)])))
ListeCpNa=sort(ListeCpNa)
ListeCp3.180=unique(as.vector(as.matrix(dPTotal2006WithMedian[,c(3+1,4+1)])))
ListeCp3.180=sort(ListeCp3.180)
ListeCp3.169=unique(as.vector(as.matrix(dPTotal2006WithMedian[,c(5+1,6+1)])))
ListeCp3.169=sort(ListeCp3.169)
ListeCp2.129=unique(as.vector(as.matrix(dPTotal2006WithMedian[,c(7+1,8+1)])))
ListeCp2.129=sort(ListeCp2.129)
ListeCp1.60=unique(as.vector(as.matrix(dPTotal2006WithMedian[,c(9+1,10+1)])))
ListeCp1.60=sort(ListeCp1.60)
ListeCp1.62=unique(as.vector(as.matrix(dPTotal2006WithMedian[,c(11+1,12+1)])))
ListeCp1.62=sort(ListeCp1.62)
ListeCp5.24=unique(as.vector(as.matrix(dPTotal2006WithMedian[,c(13+1,14+1)])))
ListeCp5.24=sort(ListeCp5.24)
ListeCp4.S=unique(as.vector(as.matrix(dPTotal2006WithMedian[,c(15+1,16+1)])))
ListeCp4.S=sort(ListeCp4.S)
ListeCp6.46=unique(as.vector(as.matrix(dPTotal2006WithMedian[,c(17+1,18+1)])))
ListeCp6.46=sort(ListeCp6.46)
ListeCp2.131=unique(as.vector(as.matrix(dPTotal2006WithMedian[,c(19+1,20+1)])))
ListeCp2.131=sort(ListeCp2.131)
ListeCp6.32=unique(as.vector(as.matrix(dPTotal2006WithMedian[,c(21+1,22+1)])))
ListeCp6.32=sort(ListeCp6.32)
ListeCp5.M=unique(as.vector(as.matrix(dPTotal2006WithMedian[,c(23+1,24+1)])))
ListeCp5.M=sort(ListeCp5.M)
ListeCp2.39=unique(as.vector(as.matrix(dPTotal2006WithMedian[,c(25+1,26+1)])))
ListeCp2.39=sort(ListeCp2.39)
ListeCp4.129=unique(as.vector(as.matrix(dPTotal2006WithMedian[,c(27+1,28+1)])))
ListeCp4.129=sort(ListeCp4.129)
ListeCpCons39=unique(as.vector(as.matrix(dPTotal2006WithMedian[,c(29+1,30+1)])))
ListeCpCons39=sort(ListeCpCons39)
ListeCPFS9VN=unique(as.vector(as.matrix(dPTotal2006WithMedian[,c(31+1,32+1)])))
ListeCPFS9VN=sort(ListeCPFS9VN)
ListeCPG3PA3=unique(as.vector(as.matrix(dPTotal2006WithMedian[,c(33+1,34+1)])))
ListeCPG3PA3=sort(ListeCPG3PA3)
ListeCPG5WUS=unique(as.vector(as.matrix(dPTotal2006WithMedian[,c(35+1,36+1)])))
ListeCPG5WUS=sort(ListeCPG5WUS)
ListeCPGHYJE=unique(as.vector(as.matrix(dPTotal2006WithMedian[,c(37+1,38+1)])))
ListeCPGHYJE=sort(ListeCPGHYJE)
ListeCPGVCLU=unique(as.vector(as.matrix(dPTotal2006WithMedian[,c(39+1,40+1)])))
ListeCPGVCLU=sort(ListeCPGVCLU)
ListeCPGZNC4=unique(as.vector(as.matrix(dPTotal2006WithMedian[,c(41+1,42+1)])))
ListeCPGZNC4=sort(ListeCPGZNC4)
ListeVerger06=unique(donneesGeno$Parcelle)
ListeVerger06=sort(ListeVerger06)
ListeVerger07=unique(donneesGeno07$Parcelle)
ListeVerger07=sort(ListeVerger07)
ListeVerger=unique(dPTotalT$Parcelle)
ListeVerger=sort(ListeVerger)

ListeCpNa
ListeCp3.180
ListeCp3.169
ListeCp2.129
ListeCp1.60
ListeCp1.62
ListeCp5.24
ListeCp4.S
ListeCp6.46
ListeCp2.131
ListeCp6.32
ListeCp5.M
ListeCp2.39
ListeCp4.129
ListeCpCons39
ListeCPFS9VN
ListeCPG3PA3
ListeCPG5WUS
ListeCPGHYJE
ListeCPGVCLU
ListeCPGZNC4

#
# A appeler pour chaque locus
# les entrees sont:
# - listeVerger: les index de vergers
# - liste Alleles: les alléles assosiés au locus
# - valeurs: les individus génotypés, avec les valeurs de alléles et la parcelle d origine
# return le tableau (un par locus au final), verger X alléle , des fréquences alléliques
# 

tableauVergersAllelesFreq<-function(listeVerger, listeAlleles, valeurs) {
  now=length(listeVerger)
  nol=length(listeAlleles)
  tableau<-matrix(0,nrow=now,ncol=nol)
  now2=nrow(valeurs)
  nwt=nrow(tableau)
  nct=ncol(tableau)
  for (ligne in 1:now2) { #pour chaque indiv
    for (collone in 2:3) { # pour chaque valeur d alléle
      for (l in 1:nwt) { #l: index de verger candidat
        for (c in 1:nct) { #c: index de valeur d alléle candidate
          #est-ce le bon couple (verger,allèles)?
          if (listeVerger[l]==valeurs[ligne,1] && listeAlleles[c]==valeurs[ligne,collone]) {
            if (listeAlleles[c]!=0) {
              tableau[l,c]=tableau[l,c]+1
            }
            else  {
              print ("bug in tableauVergersAllelesFreq, still null allèle")
              tableau[,c]<-1000
            }
          }
        }
      }
    }
  }
  for (lin in 1:nwt) {
    TT=which(tableau[lin,]!=1000)
    S=sum(tableau[lin,1:length(TT)])
    if (S==0) {S=1}
    for (co in 1:nct) {
      if (tableau[lin,co]==1000) {tableau[lin,co]<-1}
      else {tableau[lin,co]=tableau[lin,co]/S}
    }
  }
  return (tableau)
}



#application


TableauCpNaFreq<-tableauVergersAllelesFreq(ListeVerger, ListeCpNa, dPTotal2006WithMedian[,c(1:3)])
TableauCpNaFreq

TableauCp1.60Freq<-tableauVergersAllelesFreq(ListeVerger, ListeCp1.60, dPTotal2006WithMedian[,c(1,10,11)])
TableauCp1.60Freq
TableauCp1.62Freq<-tableauVergersAllelesFreq(ListeVerger, ListeCp1.62, dPTotal2006WithMedian[,c(1,12,13)])

TableauCp2.129Freq<-tableauVergersAllelesFreq(ListeVerger, ListeCp2.129, dPTotal2006WithMedian[,c(1,8,9)])

TableauCp2.131Freq<-tableauVergersAllelesFreq(ListeVerger, ListeCp2.131, dPTotal2006WithMedian[,c(1,20,21)])

TableauCp2.39Freq<-tableauVergersAllelesFreq(ListeVerger, ListeCp2.39, dPTotal2006WithMedian[,c(1,26,27)])

TableauCp3.169Freq<-tableauVergersAllelesFreq(ListeVerger, ListeCp3.169, dPTotal2006WithMedian[,c(1,6,7)])

TableauCp3.180Freq<-tableauVergersAllelesFreq(ListeVerger, ListeCp3.180, dPTotal2006WithMedian[,c(1,4,5)])

TableauCp4.129Freq<-tableauVergersAllelesFreq(ListeVerger, ListeCp4.129, dPTotal2006WithMedian[,c(1,28,29)])

TableauCp4.SFreq<-tableauVergersAllelesFreq(ListeVerger, ListeCp4.S, dPTotal2006WithMedian[,c(1,16,17)])

TableauCp5.24Freq<-tableauVergersAllelesFreq(ListeVerger, ListeCp5.24, dPTotal2006WithMedian[,c(1,14,15)])

TableauCp5.MFreq<-tableauVergersAllelesFreq(ListeVerger, ListeCp5.M, dPTotal2006WithMedian[,c(1,24,25)])

TableauCp6.32Freq<-tableauVergersAllelesFreq(ListeVerger, ListeCp6.32, dPTotal2006WithMedian[,c(1,22,23)])

TableauCp6.46Freq<-tableauVergersAllelesFreq(ListeVerger, ListeCp6.46, dPTotal2006WithMedian[,c(1,18,19)])

TableauCpCons39Freq<-tableauVergersAllelesFreq(ListeVerger, ListeCpCons39, dPTotal2006WithMedian[,c(1,30,31)])

TableauCPG3PA3Freq<-tableauVergersAllelesFreq(ListeVerger, ListeCPG3PA3, dPTotal2006WithMedian[,c(1,34,35)])

TableauCPFS9VNFreq<-tableauVergersAllelesFreq(ListeVerger, ListeCPFS9VN, dPTotal2006WithMedian[,c(1,32,33)])

TableauCPG5WUSFreq<-tableauVergersAllelesFreq(ListeVerger, ListeCPG5WUS, dPTotal2006WithMedian[,c(1,36,37)])

TableauCPGHYJEFreq<-tableauVergersAllelesFreq(ListeVerger, ListeCPGHYJE, dPTotal2006WithMedian[,c(1,38,39)])

TableauCPGVCLUFreq<-tableauVergersAllelesFreq(ListeVerger, ListeCPGVCLU, dPTotal2006WithMedian[,c(1,40,41)])

TableauCPGZNC4Freq<-tableauVergersAllelesFreq(ListeVerger, ListeCPGZNC4, dPTotal2006WithMedian[,c(1,42,43)])



#write.table(TableauCpNaFreq, "~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL/FREQUENCES/2006/1.csv", row.names=FALSE, sep=" ",dec=".", na="0")
#write.table(TableauCp3.180Freq, "~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL/FREQUENCES/2006/2.csv", row.names=FALSE, sep=" ",dec=".", na="0")
#write.table(TableauCp3.169Freq, "~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL/FREQUENCES/2006/3.csv", row.names=FALSE, sep=" ",dec=".", na="0")
#write.table(TableauCp2.129Freq, "~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL/FREQUENCES/2006/4.csv", row.names=FALSE, sep=" ",dec=".", na="0")
#write.table(TableauCp1.60Freq, "~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL/FREQUENCES/2006/5.csv", row.names=FALSE, sep=" ",dec=".", na="0")
#write.table(TableauCp1.62Freq, "~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL/FREQUENCES/2006/6.csv", row.names=FALSE, sep=" ",dec=".", na="0")
#write.table(TableauCp5.24Freq, "~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL/FREQUENCES/2006/7.csv", row.names=FALSE, sep=" ",dec=".", na="0")
#write.table(TableauCp4.SFreq, "~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL/FREQUENCES/2006/8.csv", row.names=FALSE, sep=" ",dec=".", na="0")
#write.table(TableauCp6.46Freq, "~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL/FREQUENCES/2006/9.csv", row.names=FALSE, sep=" ",dec=".", na="0")
#write.table(TableauCp2.131Freq, "~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL/FREQUENCES/2006/10.csv", row.names=FALSE, sep=" ",dec=".", na="0")
#write.table(TableauCp6.32Freq, "~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL/FREQUENCES/2006/11.csv", row.names=FALSE, sep=" ",dec=".", na="0")
#write.table(TableauCp5.MFreq, "~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL/FREQUENCES/2006/12.csv", row.names=FALSE, sep=" ",dec=".", na="0")
#write.table(TableauCp2.39Freq, "~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL/FREQUENCES/2006/13.csv", row.names=FALSE, sep=" ",dec=".", na="0")
#write.table(TableauCp4.129Freq, "~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL/FREQUENCES/2006/14.csv", row.names=FALSE, sep=" ",dec=".", na="0")
#write.table(TableauCpCons39Freq, "~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL/FREQUENCES/2006/15.csv", row.names=FALSE, sep=" ",dec=".", na="0")
#write.table(TableauCPFS9VNFreq, "~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL/FREQUENCES/2006/16.csv", row.names=FALSE, sep=" ",dec=".", na="0")
#write.table(TableauCPG3PA3Freq, "~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL/FREQUENCES/2006/17.csv", row.names=FALSE, sep=" ",dec=".", na="0")
#write.table(TableauCPG5WUSFreq, "~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL/FREQUENCES/2006/18.csv", row.names=FALSE, sep=" ",dec=".", na="0")
#write.table(TableauCPGHYJEFreq, "~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL/FREQUENCES/2006/19.csv", row.names=FALSE, sep=" ",dec=".", na="0")
#write.table(TableauCPGVCLUFreq, "~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL/FREQUENCES/2006/20.csv", row.names=FALSE, sep=" ",dec=".", na="0")
#write.table(TableauCPGZNC4Freq, "~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL/FREQUENCES/2006/21.csv", row.names=FALSE, sep=" ",dec=".", na="0")

numerot<-function(liste) {
  numbers=rep(0,length(liste)) 
  for (i in 1:length(liste)) {
    numbers[i]=i
  }
  numbers=rbind(numbers,liste)
  return (numbers)
}

ListevergerTransf=numerot(ListeVerger)
ListeCpNaTransf=numerot(ListeCpNa)
ListeCp3.180Transf=numerot(ListeCp3.180)
ListeCp3.169Transf=numerot(ListeCp3.169)
ListeCp2.129Transf=numerot(ListeCp2.129)
ListeCp1.60Transf=numerot(ListeCp1.60)
ListeCp1.62Transf=numerot(ListeCp1.62)
ListeCp5.24Transf=numerot(ListeCp5.24)
ListeCp4.STransf=numerot(ListeCp4.S)
ListeCp6.46Transf=numerot(ListeCp6.46)
ListeCp2.131Transf=numerot(ListeCp2.131)
ListeCp6.32Transf=numerot(ListeCp6.32)
ListeCp5.MTransf=numerot(ListeCp5.M)
ListeCp2.39Transf=numerot(ListeCp2.39)
ListeCp4.129Transf=numerot(ListeCp4.129)
ListeCpCons39Transf=numerot(ListeCpCons39)
ListeCPFS9VNTransf=numerot(ListeCPFS9VN)
ListeCPG3PA3Transf=numerot(ListeCPG3PA3)
ListeCPG5WUSTransf=numerot(ListeCPG5WUS)
ListeCPGHYJETransf=numerot(ListeCPGHYJE)
ListeCPGVCLUTransf=numerot(ListeCPGVCLU)
ListeCPGZNC4Transf=numerot(ListeCPGZNC4)



listerlesfich<-function(n) {
  f=c()
  for (i in 1:n){
    f[i]<-paste(i,".csv",sep="")
  }
  return(f)
}
listerlesfich(21)

preprecal<-function(donnees, listeVergers, number){
  n<-nrow(donnees)
  m<-length(listeVergers)
  tableProba<-matrix(rep(0,n*m),nrow=n,ncol=m)
  return(tableProba)
}
####le premier précalcul de boris
preCalculProba<-function(donnees, listeVergers, number) {
  T1<-Sys.time() 
  tableProba=preprecal(donnees, listeVergers, number)
  for (f in 1:number) {
    assign(paste("tableau",f,sep=""),read.csv(paste("~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL//FREQUENCES/2006/",f,".csv",sep=""),header=F,sep=" ",dec=".",quote="\""))
    assign(paste("transf",f,sep=""),read.csv(paste("~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL//CORRESPONDANCES/",f,".csv",sep=""),header=F,sep=" ",dec=".", quote="\""))
  }
  prob1=1
  prob2=1
  #pour chaque génotype 'genotype'
  for (genotype in 1:nrow(donnees)) {
    result=1
    #pour chaque sources 'source'
    for (source in 1:length(listeVergers)) {
      result=1
      file=1
      #on boucle sur les locus 'locus' 
      for (locus in seq(1,ncol(donnees)-1, by=2)) {
        
        tableauA<-get(paste("tableau",file,sep=""))
        transfA<-get(paste("transf",file,sep=""))
        file=file+1
        testfindAllele=0
        prob1=0
        prob2=0
        #on cherche la frequence des alléles via la table de correspondance
        for (al in 1:ncol(transfA)) {
          #premiere alléle
          if (donnees[genotype,locus]==transfA[2,al]) {
            prob1<-tableauA[source,al]
            testfindAllele=testfindAllele+1
            
          }
          #seconde alléle
          if (donnees[genotype,locus+1]==transfA[2,al]) {
            prob2<-tableauA[source,al]
            testfindAllele=testfindAllele+1
          }
        }
        result<-result*prob1*prob2
        if (result==0){
          break
        }
      }
      tableProba[genotype,source]<-result
    }
  }
  T2<-Sys.time()
  print (T2-T1)
  write.table(tableProba, "~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL//preCalculProba.csv", row.names=FALSE, sep=" ",dec=".", na=" ")
  return (tableProba)
}
dP07=dP07T[,2:ncol(dP07T)]
preCalculTableau=preCalculProba(dP07, ListeVerger, 21)
View(preCalculTableau)
## a rajouter
########
listerlesrep<-function(number2, ech){
  l=c()
  for (i in 1:number2) {
    l[i]=paste(ech,i,sep="")
    
    }
  return (l)
       
} 
L=listerlesrep(100,"bb")
L
listerlesfichrep<-function(number,number2,chemin,rep1){
  l=listerlesrep(number2,rep1)
  m=matrix(rep(0,number2*number),nrow=number2,ncol=number)
  for (i in 1:number2){
    for (j in 1:number){
      m[i,j]=paste(chemin,rep1,l[i],"/loc",j,"_",l[i],".txt",sep="")
    }
    
  }
    return(m)
}


#write.table(noms.alleles[1], "~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL/CORRESPONDANCES2/1.csv", row.names=FALSE, sep=" ",dec=",", na=" ")
#write.table(noms.alleles[2], "~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL/CORRESPONDANCES2/2.csv", row.names=FALSE, sep=" ",dec=",", na=" ")
#write.table(noms.alleles[3], "~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL/CORRESPONDANCES2/3.csv", row.names=FALSE, sep=" ",dec=",", na=" ")
#write.table(noms.alleles[4], "~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL/CORRESPONDANCES2/4.csv", row.names=FALSE, sep=" ",dec=",", na=" ")
#write.table(noms.alleles[5], "~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL/CORRESPONDANCES2/5.csv", row.names=FALSE, sep=" ",dec=",", na=" ")
#write.table(noms.alleles[6], "~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL/CORRESPONDANCES2/6.csv", row.names=FALSE, sep=" ",dec=",", na=" ")
#write.table(noms.alleles[7], "~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL/CORRESPONDANCES2/7.csv", row.names=FALSE, sep=" ",dec=",", na=" ")
#write.table(noms.alleles[8], "~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL/CORRESPONDANCES2/8.csv", row.names=FALSE, sep=" ",dec=",", na=" ")
#write.table(noms.alleles[9], "~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL/CORRESPONDANCES2/9.csv", row.names=FALSE, sep=" ",dec=",", na=" ")
#write.table(noms.alleles[10], "~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL/CORRESPONDANCES2/10.csv", row.names=FALSE, sep=" ",dec=",", na=" ")
#write.table(noms.alleles[11], "~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL/CORRESPONDANCES2/11.csv", row.names=FALSE, sep=" ",dec=",", na=" ")
#write.table(noms.alleles[12], "~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL/CORRESPONDANCES2/12.csv", row.names=FALSE, sep=" ",dec=",", na=" ")
#write.table(noms.alleles[13], "~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL/CORRESPONDANCES2/13.csv", row.names=FALSE, sep=" ",dec=",", na=" ")
#write.table(noms.alleles[14], "~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL/CORRESPONDANCES2/14.csv", row.names=FALSE, sep=" ",dec=",", na=" ")
#write.table(noms.alleles[15], "~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL/CORRESPONDANCES2/15.csv", row.names=FALSE, sep=" ",dec=",", na=" ")
#write.table(noms.alleles[16], "~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL/CORRESPONDANCES2/16.csv", row.names=FALSE, sep=" ",dec=",", na=" ")
#write.table(noms.alleles[17], "~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL/CORRESPONDANCES2/17.csv", row.names=FALSE, sep=" ",dec=",", na=" ")
#write.table(noms.alleles[18], "~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL/CORRESPONDANCES2/18.csv", row.names=FALSE, sep=" ",dec=",", na=" ")
#write.table(noms.alleles[19], "~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL/CORRESPONDANCES2/19.csv", row.names=FALSE, sep=" ",dec=",", na=" ")
#write.table(noms.alleles[20], "~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL/CORRESPONDANCES2/20.csv", row.names=FALSE, sep=" ",dec=",", na=" ")



# soit on le fait en 6000 ligne soit on l'automatise, 
listerlesrep2<-function(n,chemin,rep) {
  f=c()
  for (i in 1:n){
    f[i]<-paste(chemin,rep,i,sep="")
  }
  return(f)
}


listerlesfichrep<-function(number,number2,chemin,rep1){
  l=listerlesrep2(number2,rep1)
  m=matrix(rep(0,number2*number),nrow=number2,ncol=number)
  for (i in 1:number2){
    for (j in 1:number){
      m[i,j]=paste(chemin,rep1,l[i],"/loc",j,"_",l[i],".txt",sep="")
    }
    
  }
  return(m)
}



listfichcomp1<-function(n,chemin,m){
  f=c()
  for (k in 1:n){
    a<-paste(chemin,"FAcons_loc",k,"_ech",m,".txt",sep="")
    f[k]<-a
  }
  return(f)
}



listerlesfichfusion<-function(number,number2, chemin,rep1){
  l=listerlesrep2(number2, chemin)
  chemin2=paste(chemin,rep1,sep="")
  print("marche")
  m=listerlesfichrep(number,number2,chemin,rep)
  print("marche2")
  m2=matrix(rep(0,number2*number),nrow=number2,ncol=number)
  print("marche3")
  for (j in 1:number2){
    l2=listfichcomp(number,chemin,j)
    print("marche4")
    for (i in 1:number){
      print("pr")
      b=read.table(m[i,j],header=TRUE, quote="\"")
      print("marche5")
      a=read.table(l2[i], header=TRUE, quote="\"")
      c=rbind(b,a)
    }
  }
  return(c)
}




listerlesfich1<-function(n) {
  f=c()
  for (i in 1:n){
    f[i]<-paste(i,".csv",sep="")
  }
  return(f)
}

listerlesfich(100)
listerlesrep2<-function(number, rep){
  l=c()
  for (i in 1:number) {
    l[i]=paste(rep,i,sep="")
    
  }
  return (l)
  
} 



listerlesfichrep<-function(number,number2,chemin,rep1){
  l=listerlesrep(number2,rep1)
  m=matrix(rep(0,number2*number),nrow=number2,ncol=number)
  for (i in 1:number2){
    for (j in 1:number){
      m[i,j]=paste(chemin,l[i],"/loc",j,"_",l[i],".txt",sep="")
    }
    
  }
  return(m)
}



listfichcomp2<-function(n,chemin,m){
  f=c()
  for (k in 1:n){
    a<-paste(chemin,"/FAcons_loc",k,"_ech",m,".txt",sep="")
    f[k]<-a
  }
  return(f)
}



#création des dossier
#for (j in 1:100){
# e=paste("~/Bureau/Mohamed/echfinal/ech", j,sep="")
#dir.create(e)
#}

listerlesfichfusion<-function(number,number2, chemin,rep1){
  l=listerlesrep(number2, chemin)
  chemin2=paste(chemin,rep1,sep="")
  m=listerlesfichrep(number,number2,chemin,rep1)
  m2=matrix(rep(0,number2*number),nrow=number2,ncol=number)
  for (i in 1:number){
    d=paste(chemin2,i,sep="")
    l2=listfichcomp(number,d,i)
    g=paste(i,".csv",sep="")
    for (j in 1:number2){
      b=read.table(m[j,i],header=TRUE, quote="\"")
      a=read.table(l2[i], header=TRUE, quote="\"")
      c=rbind(b,a[,1])
      e=paste("~/Bureau/Mohamed/echfinal/ech", j,sep="")
      f=paste(e,g,sep="/")
      write.table(c,f , row.names=FALSE, sep=" ",dec=".", na="0")
      
      
    }
  }
  return(c)
}




d=paste("~/Bureau/Mohamed/pour_Olivier_et_Mohammed/ech",2,sep="")
e=paste("~/Bureau/Mohamed/echfinal/ech", 1,sep="")
g=paste(2,".csv",sep="")
f=paste(e,g,sep="/")
f
p<-listerlesfichfusion(20,100, "~/Bureau/Mohamed/pour_Olivier_et_Mohammed/","ech")
p
rbind( loc1_ech1,  FAcons_loc1_ech1[,1])



listerlesfich<-function(n) {
  f=c()
  for (i in 1:n){
    f[i]<-paste(i,".csv",sep="")
  }
  return(f)
}

listerlesfich(20)
listerlesrep<-function(number, rep){
  l=c()
  for (i in 1:number) {
    l[i]=paste(rep,i,sep="")
    
  }
  return (l)
  
} 


listerlesfichrep3<-function(number,number2,chemin){
  l=listerlesrep(number2,"ech")
  l2=c()
  m=matrix(rep(0,number2*number),nrow=number2,ncol=number)
  for (i in 1:number2){
    l2[i]<-paste(l[i],"/",sep="")
    for (j in 1:number){
      m[i,j]=paste(chemin,l2[i],j,".csv",sep="")
    }
    
  }
  return(m)
}


preCalculProbafinal<-function(donnees, listeVergers, number, number2) {
  T1<-Sys.time() 
  m<-listerlesfichrep(number,number2)
  tableProba=preprecal(donnees, listeVergers, number)
  m  
  for (f in 1:number){
    assign(paste("transf",f,sep=""),read.csv(paste("~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL//CORRESPONDANCES2/",f,".csv",sep=""),header=TRUE, quote="\""))
    for (i in 1:number2){
      assign(paste("tableau",f,sep=""),read.table(paste(m[i,f]),header=TRUE, quote="\""))
    }}
  #fct correctement
  prob1=1
  prob2=1 
  #pour chaque génotype 'genotype'
  for (genotype in 1:nrow(donnees)) {
    result=1
    #pour chaque sources 'source'
    for (source in 1:length(listeVergers)) {
      result=1
      file=1
      #on boucle sur les locus 'locus' 
      for (locus in seq(1,ncol(donnees)-1, by=2)) {
        
        tableauA<-get(paste("tableau",file,sep=""))
        transfA<-get(paste("transf",file,sep=""))
        file=file+1
        testfindAllele=0
        prob1=0
        prob2=0
        #on cherche la frequence des alléles via la table de correspondance
        for (al in 1:nrow(transfA)) {
          #premiere alléle
          if (donnees[genotype,locus]==transfA[al,1]) {
            prob1<-tableauA[source,al]
            testfindAllele=testfindAllele+1
            
          }
          #seconde alléle
          if (donnees[genotype,locus+1]==transfA[al,1]) {
            prob2<-tableauA[source,al]
            testfindAllele=testfindAllele+1
          }
        }
        
        #culcul du produit.
        result<-result*prob1*prob2
        
      }
      tableProba[genotype,source]<-result
    }
  }
  T2<-Sys.time()
  print (T2-T1)
  write.table(tableProba, "~/Bureau/Mohamed/ETUDE_DONNEES_DURANCE/PRECALCUL//preCalculProba3.csv", row.names=FALSE, sep=" ",dec=".", na=" ")
  return (tableProba)
}

dP07=dP07T[,4:ncol(dP07T)]
View(dP07)
preCalculTableau=preCalculProbafinal(dP07, ListeVerger, 20,100)
View(preCalculTableau)





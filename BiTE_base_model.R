# Immunological Synapse Dynamics Induced by Bispecific T-cell 
# Engagers Predict Clinical Pharmacodynamics and Tumor Evolution 
# across Anatomical Sites

# Base model create by Can Liu
# August 15, 2022

rm(list=ls())
library(dplyr)  
library(data.table)
library(deSolve) 
set.seed(1234)

# cell density
Nr<- 1E+06  # Raji cell number /mL
Nj<- 1E+06  # Jurkat cell number /mL

# receptor expression
CD19 <-  144866   # mean CD19 per cell; CD19M 
CD3 <-  66299     # mean CD3 per cell; CD3L 
GM19 <-  130670   # Geo mean CD19; CD19M 
GM3 <-  60053     # Geo mean CD3; CD3L  
SD19 <- 65092     # SD CD19M population; CD19M
SD3 <- 26763      # SD CD3L population; CD3L 

# drug information
Bconc <- 20       # Blincyto conc ng/mL
Na<- 6.02E+23     # Avogadro constant
MW <- 55000       # Blincyto Mol weight
KDA<-  2.6*10^-7  # M;KD CD3  affinity
KDB<-  1.49*10^-9 # M;KD CD19 affinity

# cell information
Rr<- 6     # um  radius of Raji
Rj <- 5    # um  radius of Jurkat
Sr= 4*pi*(Rr^2)*1.8 # um^2, surface area of Raji
Sj= 4*pi*(Rj^2)*1.8 # um^2, surface area of Jurkat 
Sc1 <- 5            # um^2, transicent contact area; 
Sc2 <- 65           # um^2, full synapse area

# Encounter probability
D<- 0.83  # um^2/s;aka 50 um^2/min, cell diffusion coefficient (50-100 um^2/min);  
Rs<- 6200  # um; spherical diameter of reaction system (1mL= (4/3)*pi*R^3)
dis= Rr+Rj    # um; distance between cell centers when encounter 
alpha=1/((Rs^3)/(3*D*dis)-0.6*(Rs^2/D)) 
fTBT <- 0.75  # encounter space coefficient for ET+E (in this file, T means E (effector cell), B means T (target cell))
fBTB <- 0.66  # encounter space coefficient for ET+T
fTBTB<- 0.75  # encounter space coefficient for ETE+T
fTBTT<- 0.5   # encounter space coefficient for ETE+E
fBTBB<- 0.33  # encounter space coefficient for ETT+T
fBTBT<- 0.67  # encounter space coefficient for ETT+E

# other
kint<- 0.002      # CD19 internalization rate constant
beta <- 0.033     # sensitive coeffcient of binding probability
NTB <- 0          # set inital value = 0, number of formed synapse (ET) 
n <- 1            # for filling Synapse dataframe
u <- 1            # for filling bond dataframe


##############################################################################################################

## drug-antigen binary binding equilibrium 
Ytot= Bconc*1000*1E-09/MW  # mol/L;drug
Atot = Nj*CD3*1000/Na  # mol/L; CD3 
Btot= Nr*CD19*1000/Na    # mol/L; CD19

aA= Atot/Ytot
aB= Btot/Ytot
ztest= aA+aB
kA= KDA/Ytot
kB= KDB/Ytot
b=  -(2+kA+kB+aA+aB) # z^3+b*z^2+c*z+d=0
c=  1+2*aA+2*aB+kA+kB+kB*aA+kA*aB+kA*kB
d= -(kB*aA+kA*aB+aA+aB)
Q=  (3*c-b^2)/9
R= (9*b*c-27*d-2*(b^3))/54
Th= acos(R/sqrt(-Q^3))
z1= -b/3+2*sqrt(-Q)*cos(Th/3) 
z2= -b/3+2*sqrt(-Q)*cos((Th+2*pi)/3)
z3= -b/3+2*sqrt(-Q)*cos((Th+4*pi)/3)

if (z1<ztest & z1>0 & z1<1){z = z1} # 0<z<ztest and z<1
if (z2<ztest & z2>0 & z2<1){z = z2}
if (z3<ztest & z3>0 & z3<1){z = z3}
print(z)

# A+Y->(KDA)->AY; B+Y->(KDB)->YB
A=Atot*KDA/(KDA+Ytot*(1-z))   # conc of free CD3   mol/L
B=Btot*KDB/(KDB+Ytot*(1-z))   # conc of free CD19
AY=(Ytot*A/KDA)/(1+(A/KDA)+(B/KDB))  # conc of AY
YB=(Ytot*B/KDB)/(1+(A/KDA)+(B/KDB))  # conc of YB
print (c(AY,YB,A,B))

# mol/L --> molecular/cell
A2 =  (A*Na)/(1000*Nj) 
B2 = (B*Na)/(1000*Nr)
AY2 = (AY*Na)/(1000*Nj) 
YB2 = (YB*Na)/(1000*Nr)
print(c(AY2,YB2,A2,B2))

## log-normal distribution of antigen expression
miux=log(GM3)  
miuy=log(GM19)   # median = GM =exp(miu)
sigx= sqrt((log(CD3)-miux)*2)     # mean = EXP(u+(sig^2)/2) 
sigy= sqrt((log(CD19)-miuy)*2)    # variance=SD^2= (EXP(sig^2)-1)*EXP(2u+sig^2)
print(c(miux,sigx,miuy,sigy))

# log-normal sampling
x <- rlnorm(n=Nj, miux, sigx) # CD3 Jurkat
y <- rlnorm(n=Nr, miuy, sigy) # CD19 Raji 
round(c(n=length(x), mean=mean(x), sd=sd(x), median=median(x)), 2) 
round(c(n=length(y), mean=mean(y), sd=sd(y), median=median(y)), 2)

# assign to each cell (molecular/um^2)
xA <- x*A2/(CD3*Sj)
xAY<- x*AY2/(CD3*Sj)
yB<- y*B2/(CD19*Sr)
yYB<- y*YB2/(CD19*Sr)

## cell-cell adhesion
# 3D Kd
kons<- 0.0001     # 3D kon (fixed); um^3/molecular*s^-1 fixed equals to 6E+04 M^-1*s^-1
koffAs = kons*KDA*Na/1E+15   # 3D koff CD3 s^-1 --> kons: um^3/molecular*s^-1  KoffAs CD3
koffBs = kons*KDB*Na/1E+15   # 3D koff CD19 s^-1 --> kons: um^3/molecular*s^-1  koffBs CD19
print (c(koffAs,koffBs))

# 3D Diffusion coefficient
Ds19 = Ds3 = DsY = 50  # 3D diffusion coefficient um^2/s;  5E-07 cm^2/s for Ab
Dm19<- 0.006      # um^2/s; 2D membrane diffusion coefficient of CD19
Dm3<- 0.01-Dm19   # Dm19+Dm3=0.01
Rab <- 0.005      # um; encounter distance

# 2D diffusion rate constant
dsf = 4*pi*(Ds19+DsY)*Rab  # um^3/s forward rate constant in solution
dsr = 3*(Ds19+DsY)/(Rab^2)  # s^-1  reverse rate constant in solution
dmf= 2*pi*(Dm19+Dm3)    # 2D on membrane, forward 
dmr= 2*(Dm19+Dm3)/(Rab^2)  # 2D on membrane, reverse 
print(c(dsf,dsr,dmf,dmr))

# rotation
esr=(3*dsr)/(4*pi^2)    # rotation/conformation  
esf=0.04*esr            # Es= esf/esr=0.04
emr=(3*dmr)/(4*pi^2)
emf=0.04*emr            # Em= emf/emr=0.04

# calculate 2D kon koff
rsfA= (kons*dsr*esr)/(dsf*esf-kons*(dsr+esf))
rsfB= (kons*dsr*esr)/(dsf*esf-kons*(dsr+esf))
rsrA= koffAs*(dsr*esr+rsfA*(dsr+esf))/(dsr*esr)
rsrB= koffBs*(dsr*esr+rsfB*(dsr+esf))/(dsr*esr)
print(c(rsfA,rsrA,rsfB,rsrB))

konA = dmf*rsfA*emf/(dmr*emr+rsfA*(dmr+emf))   # CD3
koffA = dmr*rsrA*emr/(dmr*emr+rsfA*(dmr+emf))   # CD3
konB = dmf*rsfB*emf/(dmr*emr+rsfB*(dmr+emf))    # CD19
koffB = dmr*rsrB*emr/(dmr*emr+rsfB*(dmr+emf))   # CD19

params  <- c(konA,konB,koffA,koffB)
print(params)

# ODE function for ternary complex (CD3-BiTE-CD19)
ode_binding <- function (time,x,params){
  AYB <- x[1]   
  A <- x[2]     # previously A B AY and YB have been used to indicate free and binded antigen conc after binary binding in solution
  AY <- x[3]    # should be okay to re-define those terms here.
  B <- x[4]
  YB <- x[5]
  dAYBdt = konA*A*YB+konB*B*AY-koffA*AYB-koffB*AYB
  dAdt   = koffA*AYB-konA*A*YB
  dAYdt  = koffB*AYB-konB*B*AY
  dBdt   = koffB*AYB-konB*B*AY
  dYBdt  = koffA*AYB-konA*A*YB
  
list(c(dAYBdt,dAdt,dAYdt,dBdt,dYBdt))}

timeseries <- seq(from = 0, to= 60, by= 0.1)

###############################################################################################################

# dataframes for Raji, Jurkat, synapse and output
Nmax <- pmax(Nj,Nr) 
Jurkat <- data.frame(1:Nj, rep("single",Nj),xA,xAY) 
Raji <- data.frame(1:Nr, rep("single",Nr),yB,yYB) 
synapse <- as.data.frame(matrix(NA,  nrow = Nmax, ncol = 15))
bond <- as.data.frame(matrix(NA,  nrow = Nmax, ncol = 6))     # to check bond (ternary complex) generation
output <- as.data.frame(matrix(NA,  nrow = 60, ncol = 30))    # nrow number of rounds (60 rounds/hour)

names(Jurkat)<- c("ID","status","xA","xAY")
names(Raji)<- c("ID","status","yB","yYB")
names(synapse) <-c("IDT1","IDB1","status","xA","xAY","yB","yYB","IDT2","IDB2","xA2","xAY2","yB2","yYB2","IDT3","IDB3")
names(bond)<- c("Pb","AYB1","AYB2","NAYB1","NAYB2","bondtime")
names(output)<- c("NTB","NTBT","NBTB","Nmulti","Nsynapse","Tperc","Bperc","TBperc","Tmulti",
                  "Bmulti","TBmulti","multiofTB","TBratio","TBmultiratio","Pse","PmeTBT","PmeBTB","Pb","bond5s","bondmean","bondsum","bondTime",
                  "MeanxA","MeanxAY","MeanyB","MeanyYB","NTBTT","NTBTB","NBTBT","NBTBB")

################################################################################################################

for (t in 1:60) {   # 60 rounds for 1 hour simulation
  
  print(t)   
  print (output)   
  
  # CD3 downmodulation renew
  Jurkat$xA <- xA/(1+(xAY^0.9)*(t^0.7)*0.1)
  Jurkat$xAY <- xAY/(1+(xAY^0.9)*(t^0.7)*0.1)
  synapseTB <- synapse[which(synapse$status=='TB'),] 
  synapseTBT <- synapse[which(synapse$status=='TBT'),] 
  
  if(nrow(synapseTB)>0)  { 
    names(Jurkat)[names(Jurkat)=="ID"] <- "IDT1"
    Jurkat$IDT1 <- as.character(Jurkat$IDT1)
    Jurkat$xA<- as.character(Jurkat$xA)
    Jurkat$xAY<- as.character(Jurkat$xAY)
    synapse <- left_join(synapse, Jurkat, by = "IDT1")
    synapse <- select(synapse, IDT1, IDB1, status = status.x,xA=xA.y,xAY=xAY.y,yB,yYB,IDT2,IDB2,xA2,xAY2,yB2,yYB2,IDT3,IDB3)
    names(Jurkat)[names(Jurkat)=="IDT1"] <- "ID" 
    
     if(nrow(synapseTBT)>0) {
       names(Jurkat)[names(Jurkat)=="ID"] <- "IDT2"
       Jurkat$IDT2 <- as.character(Jurkat$IDT2)
       synapse <- left_join(synapse, Jurkat, by = "IDT2")
       synapse <- select(synapse, IDT1, IDB1, status = status.x,xA=xA.x,xAY=xAY.x,yB,yYB,IDT2,IDB2,xA2=xA.y,xAY2=xAY.y,yB2,yYB2,IDT3,IDB3)
       names(Jurkat)[names(Jurkat)=="IDT2"] <- "ID"}
   
    Jurkat$xA<-as.numeric(Jurkat$xA)
    Jurkat$xAY<-as.numeric(Jurkat$xAY)} 
  
  # CD19 internalization renew
  Raji$yYB <- Raji$yYB*exp(-kint)  # yYB is not original yYB, was renewed by last loop, so not exp(-kint*t)
  synapse$yYB <- as.numeric(as.character(synapse$yYB)) 
  synapse$yYB2 <- as.numeric(as.character(synapse$yYB2))
  synapse$yYB<- synapse$yYB*exp(-kint)
  synapse$yYB2<- synapse$yYB2*exp(-kint)
  
  # renew Rajisingle and Jurkatsingle if any change in previous rounds to calculate Probability
  Rajisingle <- Raji[which(Raji$status=='single'),]
  Jurkatsingle <- Jurkat[which(Jurkat$status=='single'),] 
  synapseTB <- synapse[which(synapse$status=='TB'),] 
  synapseTBT <- synapse[which(synapse$status=='TBT'),]
  synapseBTB <- synapse[which(synapse$status=='BTB'),]
  
  # renew cell-cell encounter probability for each round 
  Pse  <-  1-exp(-alpha*nrow(Rajisingle)*60)                 # encounter probability for ET 
  PmeBTB <-  (1-exp(-alpha*nrow(Rajisingle)*60))*fBTB*2      # encounter probability for ET+T
  PmeTBTB <-  (1-exp(-alpha*nrow(Rajisingle)*60))*fTBTB*3    # encounter probability for ETT+E
  PmeBTBB <-  (1-exp(-alpha*nrow(Rajisingle)*60))*fBTBB*3    # encounter probability for ETT+T
  PmeTBT <-  (1-exp(-alpha*nrow(Jurkatsingle)*60))*fTBT*2    # encounter probability for ET+E
  PmeTBTT <-  (1-exp(-alpha*nrow(Jurkatsingle)*60))*fTBTT*3  # encounter probability for ETE+E
  PmeBTBT <-  (1-exp(-alpha*nrow(Jurkatsingle)*60))*fBTBT*3  # encounter probability for ETE+T
  
  if(NTB>0){  # initiate only TB synapse number > 0 to form three-cell cluster, skip when NTB=0 at 1st round
    
    for (i in 1:nrow(synapseTB)){   # for ET + additional T or E
      
      test1 <- runif(1) 
      
      if(test1 <= PmeBTB){   # encounter ET+T
        
        Rajisingle <- Raji[which(Raji$status=='single'),]
        if(nrow(Rajisingle)==0){break}  
        
        Rrow <- Rajisingle[sample(nrow(Rajisingle),1,replace=TRUE), ] # randomly sample a free single Raji cell
        
        xstart <- c(0,as.numeric(synapseTB[i,4]),as.numeric(synapseTB[i,5]),Rrow[,3],Rrow[,4])
        out<- as.data.frame(ode(times = timeseries,
                                y=xstart,
                                func = ode_binding, 
                                parms = params))
        
        bondtime <-round(runif(1,min=0.1,max=5),digits=1) # random encounter duration among 0.1-5s
        outrow <- bondtime*10+1   # convert to row # in out[,] 
        AYB1 <- out[outrow,2]     # AYB bond conc at bondtime
        AYB2 <- out[601,2]        # AYB conc at 60s
        NAYB1  = AYB1*Sc1         # total AYB bond at bondtime
        NAYB2  = AYB2*Sc2         # total final AYB bond at 60s # just check, not used to calculate adhesion prob
        Pb1 <- 1-exp(-beta*NAYB1) # binding probability 
        
        test2 <- runif(1)
        
        if (test2 < Pb1){      # cell adhesion  ET+T
          synapse$status <- as.character(synapse$status)
          synapse[as.numeric(rownames(synapseTB[i,])),]$status<-"BTB"  # change status from "TB" to "BTB"
          synapse[as.numeric(rownames(synapseTB[i,])),9]<- Rrow$ID     # fill binded Raji cell ID to synapse dataframe 
          synapse[as.numeric(rownames(synapseTB[i,])),12]<- Rrow[,3]   # fill "yB2"
          synapse[as.numeric(rownames(synapseTB[i,])),13]<- Rrow[,4]   # fill "yYB2"
          Raji$status <- as.character(Raji$status)  
          Raji[as.numeric(rownames(Rrow)),]$status <- "BTB"  # change Raji cell status from "single" to "BTB" in Raji dataframe
          bond[u,]<- c(Pb1,AYB1,AYB2,NAYB1,NAYB2,bondtime)   # fill bond dataframe 
          u<- u+1
          
        }    # if Pb1    
      }     # if PmeBTB
      
      if (test1 > PmeBTB & test1 <= PmeBTB+PmeTBT){  # encounter ET+E
        
        Jurkatsingle <- Jurkat[which(Jurkat$status=='single'),] 
        if(nrow(Jurkatsingle)==0){break}  
        
        Jrow <- Jurkatsingle[sample(nrow(Jurkatsingle),1,replace=TRUE), ] 
        xstart <- c(0,Jrow[,3],Jrow[,4],as.numeric(synapseTB[i,6]),as.numeric(synapseTB[i,7]))
        out<- as.data.frame(ode(times = timeseries,
                                y=xstart,
                                func = ode_binding, 
                                parms = params))
        
        bondtime <-round(runif(1,min=0.1,max=5),digits=1) 
        outrow <- bondtime*10+1 
        AYB1 <- out[outrow,2]  
        AYB2 <- out[601,2] 
        NAYB1  = AYB1*Sc1  
        NAYB2  = AYB2*Sc2 
        Pb2 <- 1-exp(-beta*NAYB1) 
        
        test3 <- runif(1)
        
        if (test3 < Pb2){  # cell adhesion ET+E
          synapse$status <- as.character(synapse$status)
          synapse[as.numeric(rownames(synapseTB[i,])),]$status<-"TBT"
          synapse[as.numeric(rownames(synapseTB[i,])),8]<- Jrow$ID
          synapse[as.numeric(rownames(synapseTB[i,])),10]<- Jrow[,3]   # fill "xA2"
          synapse[as.numeric(rownames(synapseTB[i,])),11]<- Jrow[,4]    # fill "xAY2
          Jurkat$status <- as.character(Jurkat$status) 
          Jurkat[as.numeric(rownames(Jrow)),]$status <- "TBT"  
          bond[u,]<- c(Pb2,AYB1,AYB2,NAYB1,NAYB2,bondtime) 
          u<- u+1
          
        }   # if Pb2  
      }   # if PmeTBT
    }     #  for i
  }    # NTB>0
  
  # renew Jurkatsingle & Rajisingle if any change 
  Jurkatsingle <- Jurkat[which(Jurkat$status=='single'),]  
  Rajisingle <- Raji[which(Raji$status=='single'),]
  synapseTBT <- synapse[which(synapse$status=='TBT'),]
  synapseBTB <- synapse[which(synapse$status=='BTB'),]
  synapseTB <- synapse[which(synapse$status=='TB'),] 
  
  # form four-cell cluster
  if(nrow(synapseTBT)>0){  # add the additional E or T to ETE
    
    for (p in 1:nrow(synapseTBT)){   
      
    Rajisingle <- Raji[which(Raji$status=='single'),]
    Jurkatsingle <- Jurkat[which(Jurkat$status=='single'),] 
    if(nrow(Rajisingle)==0){break}
    if(nrow(Jurkatsingle)==0){break} 
    
    test4 <- runif(1)
    
    if(test4 <= PmeTBTB){  # cell encounter ETE + T
      
      Rrow <- Rajisingle[sample(nrow(Rajisingle),1,replace=TRUE), ] 
      
      test5 <- runif(1)
      
        if (test5<= 0.5) {  # decide which E binds to additional T
          
      xstart <- c(0,as.numeric(synapseTBT[p,10]),as.numeric(synapseTBT[p,11]),Rrow[,3],Rrow[,4])
      out<- as.data.frame(ode(times = timeseries,
                              y=xstart,
                              func = ode_binding, 
                              parms = params))
      
      bondtime <-round(runif(1,min=0.1,max=5),digits=1) 
      outrow <- bondtime*10+1 
      AYB1 <- out[outrow,2]   
      AYB2 <- out[601,2]    
      NAYB1  = AYB1*Sc1       
      NAYB2  = AYB2*Sc2         
      Pb3 <- 1-exp(-beta*NAYB1) 
                     
                          } # test5
      
      if (test5>0.5) { 
      xstart <- c(0,as.numeric(synapseTBT[p,4]),as.numeric(synapseTBT[p,5]),Rrow[,3],Rrow[,4])
        out<- as.data.frame(ode(times = timeseries,
                                y=xstart,
                                func = ode_binding, 
                                parms = params))
        
        bondtime <-round(runif(1,min=0.1,max=5),digits=1) 
        outrow <- bondtime*10+1    
        AYB1 <- out[outrow,2]     
        AYB2 <- out[601,2]       
        NAYB1  = AYB1*Sc1        
        NAYB2  = AYB2*Sc2         
        Pb3 <- 1-exp(-beta*NAYB1) 
            
                         }  # else
      
          test6 <- runif(1)
          
          if (test6 < Pb3){   # cell adhesion ETE+T
            synapse$status <- as.character(synapse$status)
            synapse[as.numeric(rownames(synapseTBT[p,])),]$status<-"TBTB"
            synapse[as.numeric(rownames(synapseTBT[p,])),15]<- Rrow$ID
            Raji$status <- as.character(Raji$status)  
            Raji[as.numeric(rownames(Rrow)),]$status <- "TBTB"  
            bond[u,]<- c(Pb3,AYB1,AYB2,NAYB1,NAYB2,bondtime) 
            u<- u+1
                          }  #  if Pb3
                          }  #  if PmeTBTB
  
          if (test4> PmeTBTB & test4 <= PmeTBTB+PmeTBTT){     # cell encounter ETE+E
            
            Jrow <- Jurkatsingle[sample(nrow(Jurkatsingle),1,replace=TRUE), ]
          
            xstart <- c(0,Jrow[,3],Jrow[,4],as.numeric(synapseTBT[p,6]),as.numeric(synapseTBT[p,7]))
            out<- as.data.frame(ode(times = timeseries,
                                    y=xstart,
                                    func = ode_binding, 
                                    parms = params))
            
            bondtime <-round(runif(1,min=0.1,max=5),digits=1) 
            outrow <- bondtime*10+1 
            AYB1 <- out[outrow,2]  
            AYB2 <- out[601,2] 
            NAYB1  = AYB1*Sc1  
            NAYB2  = AYB2*Sc2 
            Pb4 <- 1-exp(-beta*NAYB1) 
            
            test7<- runif(1)
            
            if (test7 < Pb4){   # cell adhesion ETE+E
              synapse$status <- as.character(synapse$status)
              synapse[as.numeric(rownames(synapseTBT[p,])),]$status<-"TBTT"
              synapse[as.numeric(rownames(synapseTBT[p,])),14]<- Jrow$ID
              Jurkat$status <- as.character(Jurkat$status) 
              Jurkat[as.numeric(rownames(Jrow)),]$status <- "TBTT"  
              bond[u,]<- c(Pb4,AYB1,AYB2,NAYB1,NAYB2,bondtime) 
              u<- u+1
              
                             } # if Pb4
                             } # PmeTBTB+PmeTBTT
    }  # p
    }  # if nrow(synapseTBT)>0
  
  # renew Jurkatsingle & Rajisingle if any change
  Jurkatsingle <- Jurkat[which(Jurkat$status=='single'),]  
  Rajisingle <- Raji[which(Raji$status=='single'),]
  synapseTBT <- synapse[which(synapse$status=='TBT'),]
  synapseBTB <- synapse[which(synapse$status=='BTB'),]
  synapseTB <- synapse[which(synapse$status=='TB'),] 
  

  if(nrow(synapseBTB)>0){   # add an additional E or T to TET
    
    for (q in 1:nrow(synapseBTB)){   
      
      Rajisingle <- Raji[which(Raji$status=='single'),]
      Jurkatsingle <- Jurkat[which(Jurkat$status=='single'),] 
      if(nrow(Rajisingle)==0){break}
      if(nrow(Jurkatsingle)==0){break} 
      
      test8 <- runif(1)
      
      if(test8 <= PmeBTBT){  # cell encounter TET+E
        
        Jrow <- Jurkatsingle[sample(nrow(Jurkatsingle),1,replace=TRUE), ]
        
        test9 <- runif(1)
        
        if (test9<= 0.5) {  # decide which T binds to additional E
          
          xstart <- c(0,Jrow[,3],Jrow[,4],as.numeric(synapseBTB[q,6]),as.numeric(synapseBTB[q,7]))
          out<- as.data.frame(ode(times = timeseries,
                                  y=xstart,
                                  func = ode_binding, 
                                  parms = params))
          
          bondtime <-round(runif(1,min=0.1,max=5),digits=1) 
          outrow <- bondtime*10+1 
          AYB1 <- out[outrow,2]  
          AYB2 <- out[601,2] 
          NAYB1  = AYB1*Sc1  
          NAYB2  = AYB2*Sc2 
          Pb5 <- 1-exp(-beta*NAYB1) 
                       } # test9
        
        if (test9 > 0.5) { 
          
          xstart <- c(0,Jrow[,3],Jrow[,4],as.numeric(synapseBTB[q,12]),as.numeric(synapseBTB[q,13]))
          out<- as.data.frame(ode(times = timeseries,
                                  y=xstart,
                                  func = ode_binding, 
                                  parms = params))
          
          bondtime <-round(runif(1,min=0.1,max=5),digits=1) 
          outrow <- bondtime*10+1 
          AYB1 <- out[outrow,2]  
          AYB2 <- out[601,2] 
          NAYB1  = AYB1*Sc1  
          NAYB2  = AYB2*Sc2 
          Pb5 <- 1-exp(-beta*NAYB1) 
                         } # else
        
          test10 <- runif(1)
          
          if (test10 < Pb5){   # cell adhension TET+E
          
          synapse$status <- as.character(synapse$status)
          synapse[as.numeric(rownames(synapseBTB[q,])),]$status<-"BTBT"
          synapse[as.numeric(rownames(synapseBTB[q,])),14]<- Jrow$ID
          Jurkat$status <- as.character(Jurkat$status) 
          Jurkat[as.numeric(rownames(Jrow)),]$status <- "BTBT"  
          bond[u,]<- c(Pb5,AYB1,AYB2,NAYB1,NAYB2,bondtime) 
          u<- u+1
                          } # Pb5
                          } # PmeBTBT
      
      if (test8> PmeBTBT & test8 <= PmeBTBT+PmeBTBB){     # cell encounter  TET+T
        
        Rrow <- Rajisingle[sample(nrow(Rajisingle),1,replace=TRUE), ] 
        
        xstart <- c(0,as.numeric(synapseBTB[q,4]),as.numeric(synapseBTB[q,5]),Rrow[,3],Rrow[,4])
        out<- as.data.frame(ode(times = timeseries,
                                y=xstart,
                                func = ode_binding, 
                                parms = params))
        bondtime <-round(runif(1,min=0.1,max=5),digits=1) 
        outrow <- bondtime*10+1    
        AYB1 <- out[outrow,2]     
        AYB2 <- out[601,2]        
        NAYB1  = AYB1*Sc1        
        NAYB2  = AYB2*Sc2         
        Pb6 <- 1-exp(-beta*NAYB1) 
        
        test11 <- runif(1)
        
        if (test11 < Pb6){  # cell adhesion TET+T
          synapse$status <- as.character(synapse$status)
          synapse[as.numeric(rownames(synapseBTB[q,])),]$status<-"BTBB"
          synapse[as.numeric(rownames(synapseBTB[q,])),15]<- Rrow$ID
          Raji$status <- as.character(Raji$status)  
          Raji[as.numeric(rownames(Rrow)),]$status <- "BTBB"  
          bond[u,]<- c(Pb6,AYB1,AYB2,NAYB1,NAYB2,bondtime) 
          u<- u+1
          
                         } #Pb6
                         } # PmeBTBT+PmeBTBB
    }  # q
    }  # if nrow(synapseBTB)>0
  
  # renew Jurkatsingle & Rajisingle if any change 
  Jurkatsingle <- Jurkat[which(Jurkat$status=='single'),]  
  Rajisingle <- Raji[which(Raji$status=='single'),]
  
  # free single E + free single T --> ET
  if(nrow(Jurkatsingle)> 0 & nrow(Rajisingle)> 0) {   
    
    for (j in 1:nrow(Jurkatsingle)){
      
      test12 <- runif(1) 
      
      if (test12 < Pse){  # cell encounter E+T
        
        Rrow <- Rajisingle[sample(nrow(Rajisingle),1,replace=TRUE), ] 
        xstart <- c(0,Jurkatsingle[j,3],Jurkatsingle[j,4],Rrow[,3],Rrow[,4])
        out<- as.data.frame(ode(times = timeseries,
                                y=xstart,
                                func = ode_binding, 
                                parms = params))
        
        bondtime <-round(runif(1,min=0.1,max=5),digits=1)
        outrow <- bondtime*10+1  
        AYB1 <- out[outrow,2]  
        AYB2 <- out[601,2] 
        NAYB1  = AYB1*Sc1  
        NAYB2  = AYB2*Sc2 
        Pb7 <- 1-exp(-beta*NAYB1) 
        
        test13 <- runif(1)
        if (test13 < Pb7){   # cell adhension E+T
          synapse[n,]<-c(Jurkatsingle[j,1],Rrow[,1],"TB",Jurkatsingle[j,3],Jurkatsingle[j,4],Rrow[,3],Rrow[,4],NA,NA,NA,NA,NA,NA,NA,NA) # to fill synapse df
          Jurkat$status <- as.character(Jurkat$status) 
          Jurkat[Jurkatsingle[j,]$ID,]$status<- "TB"  
          Raji$status <- as.character(Raji$status)
          Raji[as.numeric(rownames(Rrow)),]$status<-"TB" 
          bond[u,]<- c(Pb7,AYB1,AYB2,NAYB1,NAYB2,bondtime) # to fill bond df
          n=n+1
          u=u+1
          
        }   # Pb7
      }   # Pse
    }   # for j
  }   # if Jurkatsingle > 0
  
  NTB <- sum(synapse$status == "TB",na.rm=TRUE)
  NTBT <- sum(synapse$status == "TBT",na.rm=TRUE)
  NBTB <- sum(synapse$status == "BTB",na.rm=TRUE)
  NTBTT <- sum(synapse$status == "TBTT",na.rm=TRUE)
  NTBTB <- sum(synapse$status == "TBTB",na.rm=TRUE)
  NBTBT <- sum(synapse$status == "BTBT",na.rm=TRUE)
  NBTBB <- sum(synapse$status == "BTBB",na.rm=TRUE)
  
  Nmulti<- NTBT+NBTB+NTBTT+NTBTB+NBTBT+NBTBB
  Nsynapse <- NTB+Nmulti
  Tperc <- (Nsynapse+NTBT+NTBTT*2+NTBTB+NBTBT)*100/Nj  # E% involved (Jurkat)
  Bperc <- (Nsynapse+NBTB+NTBTB+NBTBT+NBTBB*2)*100/Nr  # T% involved (Raji)
  TBperc <- (Nmulti*4+NTB*2-NTBT-NBTB)*100/(Nr+Nj)
  Tmulti<-  (Nmulti+NTBT+NTBTT*2+NTBTB+NBTBT)*100/Nj
  Bmulti<-  (Nmulti+NBTB+NTBTB+NBTBT+NBTBB*2)*100/Nr
  TBmulti<-  (Nmulti*4-NTBT-NBTB)*100/(Nr+Nj)
  multiofTB<- (Nmulti*4-NTBT-NBTB)*100/(NTB*2+Nmulti*4-NTBT-NBTB)
  
  if (NTB>0|NTBT>0|NBTB>0|NTBTT>0|NTBTB>0|NBTBT>0|NBTBB>0)
    {TBratio <- (Nsynapse+NTBT+NTBTT*2+NTBTB+NBTBT)/(Nsynapse+NBTB+NTBTB+NBTBT+NBTBB*2)}
  else {TBratio <- NA}
  
  if(NTBT>0|NBTB>0|NTBTT>0|NTBTB>0|NBTBT>0|NBTBB>0)   
     {TBmultiratio <- (Nmulti+NTBT+NTBTT*2+NTBTB+NBTBT)/(Nmulti+NBTB+NTBTB+NBTBT+NBTBB*2)}
  else {TBmultiratio <- NA}
  
  Pb <- as.numeric(colMeans(bond["Pb"], na.rm = TRUE))
  bond5s<-  as.numeric(colMeans(bond["NAYB1"], na.rm = TRUE))
  bondmean <-  as.numeric(colMeans(bond["NAYB2"], na.rm = TRUE))
  bondsum <-  as.numeric(colSums(bond["NAYB2"], na.rm = TRUE))
  bondTime <-  as.numeric(colMeans(bond["bondtime"], na.rm = TRUE))
  MeanxA <-  as.numeric(colMeans(Jurkat["xA"], na.rm = TRUE))
  MeanxAY <-  as.numeric(colMeans(Jurkat["xAY"], na.rm = TRUE))
  MeanyB <-  as.numeric(colMeans(Raji["yB"], na.rm = TRUE))
  MeanyYB <-  as.numeric(colMeans(Raji["yYB"], na.rm = TRUE))
  
  output[t,] <- c(NTB,NTBT,NBTB,Nmulti,Nsynapse,Tperc,Bperc,TBperc,
                  Tmulti,Bmulti,TBmulti,multiofTB,TBratio,TBmultiratio,Pse,PmeTBT,PmeBTB,Pb,
                  bond5s,bondmean,bondsum,bondTime,MeanxA,MeanxAY,MeanyB,MeanyYB,NTBTT,NTBTB,NBTBT,NBTBB)  
  
}  # t
print (output)

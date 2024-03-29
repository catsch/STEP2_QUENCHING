########################################################################
## Quenching Correction based on Louis Terrats Work 
##
## Inputs : Dark corrected CHLA profiles 
## 
## Version : 13 January 2021
## Version : 22 November 2021
## 
## Bug Correction : light threshold and same value in the first 10m
########################################################################

library(ncdf4)
library(stringr)

source("./MLD_calc.R")

source("./IPAR_15_DEPTH.R")

source("./RunningFilter.R")

uf=commandArgs()

file_dark  <- uf[2]

### reading the dark value

DARK=read.table(as.character(file_dark),header=FALSE, as.is=TRUE)

DARK_CHLA=as.numeric(DARK$V1)

#### Creating the list of files for which we need to recompute  
liste_to_do=read.table("liste_all_B",header=FALSE, as.is=TRUE)

# List of the file to process
LIST_nc=liste_to_do$V1

for (IDnc in LIST_nc) {

### IDnc 

	print(IDnc)

### Initialisation FLAG_SHALLOW

	FLAG_SHALLOW=FALSE

### Initialisation FLAG_SHALLOW

	FLAG_NO_PAR=FALSE
 
#### Getting the name of the core file
	file_in_C=str_replace(IDnc,"/B","/")

# if B and C are not in the same mode 
	if (!file.exists(file_in_C)) file_in_C=str_replace(file_in_C,"profiles/R","profiles/D")
	if (!file.exists(file_in_C)) file_in_C=str_replace(file_in_C,"profiles/D","profiles/R")

#########################
# Open the nc File
#########################

# Open the C file
	filenc_C=nc_open(file_in_C,readunlim=FALSE,write=FALSE)

# Open the B file
	filenc_B=nc_open(IDnc,readunlim=FALSE,write=TRUE)

############################################################
## work on index 
#############################################################

#### Get the list of parameters in the profile

	STATION_PARAMETERS=ncvar_get(filenc_B,"STATION_PARAMETERS")

# Stations parameters has a fixed length 64 characters 

	CHLA_STRING=str_pad("CHLA",64,"right")

	PAR_STRING=str_pad("DOWNWELLING_PAR",64,"right")

# Find the profile containing CHLA/BBP and PAR  

	index_chla=which(STATION_PARAMETERS == CHLA_STRING, arr.ind=TRUE)

	index_par=which(STATION_PARAMETERS == PAR_STRING, arr.ind=TRUE)

	if (length(index_chla)==0) {

		next # jump to the next profile if no chla in the profile

	} else {

		iprof_chla=index_chla[,2]

	} 

	index_par=which(STATION_PARAMETERS == PAR_STRING, arr.ind=TRUE)	

	if (length(index_par)==0) { 

		FLAG_NO_PAR=TRUE

	} else {

		iprof_par=index_par[,2]
	}	

##########################################################
# MLD estimation
##########################################################

#### Read the C file to estimate the MLD

	TEMP_CTD=ncvar_get(filenc_C,"TEMP")

	PSAL_CTD=ncvar_get(filenc_C,"PSAL")

	PRES_CTD=ncvar_get(filenc_C,"PRES")

	MLD=MLD_calc(PRES_CTD, PSAL_CTD , TEMP_CTD)

###############################################################################
######## Jump to next profile or dont take PAR into account if PRES_QC is too bad
###############################################################################

	PROFILE_PRES_CTD_QC=strsplit(ncvar_get(filenc_C,"PROFILE_PRES_QC"),split="")

	if (PROFILE_PRES_CTD_QC[[1]][iprof_chla]=="E" | PROFILE_PRES_CTD_QC[[1]][iprof_chla]=="F") next

	if ( !FLAG_NO_PAR ) {
		if (PROFILE_PRES_CTD_QC[[1]][iprof_par]=="E" | PROFILE_PRES_CTD_QC[[1]][iprof_par]=="F") FLAG_NO_PAR=TRUE
	}

###########################################################
# Light Threshold
###########################################################

#### ipar_15_depth

	if ( !FLAG_NO_PAR ) {

		DOWNWELLING_PAR=ncvar_get(filenc_B,"DOWNWELLING_PAR")

##FEV22		Filtering the par to get the correct depth of ipar15

		FILTER_PAR=DOWNWELLING_PAR

		FILTER_PAR[,iprof_par]=RunningFilter(2,DOWNWELLING_PAR[,iprof_par],na.fill=T, ends.fill=T, Method="Median")
		
		if (length(which(!is.na(DOWNWELLING_PAR[,iprof_par]))) < 2 ) FLAG_NO_PAR=TRUE

		ipar_15_depth=IPAR_15_DEPTH(PRES_CTD,DOWNWELLING_PAR)
		print(ipar_15_depth)

		ipar_15_depth=IPAR_15_DEPTH(PRES_CTD,FILTER_PAR)

		print(ipar_15_depth)

	} else {

		ipar_15_depth=MLD

	}

###########################################################
# Testing shallow mixing 
###########################################################

	if (ipar_15_depth > MLD) FLAG_SHALLOW=TRUE

############################################################
# Get PAR and CHLA on the same levels 
############################################################

#### Working on PAR
	if (!FLAG_NO_PAR) {
 
		PAR_CHLA=approx(PRES_CTD[,iprof_par],FILTER_PAR[,iprof_par],PRES_CTD[,iprof_chla], rule=2)$y

#### Filtering a bit the PAR 
		MED_PAR_CHLA=RunningFilter(2,PAR_CHLA,na.fill=T, ends.fill=T, Method="Median")

	} 

################################################################
#### Getting BBP700 and CHLA 
################################################################ 

###  Getting CHLA

	CHLA=ncvar_get(filenc_B,"CHLA")

###  DARK correction with input from quentin's Estimation 
	CHLA=CHLA-DARK_CHLA

###  Get the CHLA_QC and change them into CHLA_ADJUSTED_QC
	CHLA_ADJUSTED_QC=ncvar_get(filenc_B,"CHLA_QC") 

###  PROFILE_CHLA

	PROFILE_CHLA_QC=str_squish(ncvar_get(filenc_B,"PROFILE_CHLA_QC"))

###  Initialising the CHLA derived variables without the quenching 

	MED_CHLA=CHLA # unspiked CHLA

	CHLA_NPQ=CHLA 	# final CHLA corrected from NPQ and SPIKED

	CHLA_NPQ_D=CHLA	# CHLA corrected from NPQ and despiked

###  Getting BBP 

	BBP700=ncvar_get(filenc_B,"BBP700")

	MED_BBP700=BBP700 # Filtered BBP



############################################################################################
# Filtering the signal To remove the spikes for BBP and CHLA (Briggs et al.) 
#	1. median filter 5 
#	2. Mean filter 7 
############################################################################################

### 1st median filter 5

	MED_CHLA[,iprof_chla]=RunningFilter(2,CHLA[,iprof_chla],na.fill=T, ends.fill=T, Method="Median")

	MED_BBP700[,iprof_chla]=RunningFilter(2,BBP700[,iprof_chla],na.fill=T, ends.fill=T, Method="Median")

### 2nd Mean filter 7 for BBP and CHLA 

	MED_CHLA[,iprof_chla]=RunningFilter(3,MED_CHLA[,iprof_chla],na.fill=T, ends.fill=T, Method="Mean")

	MED_BBP700[,iprof_chla]=RunningFilter(3,MED_BBP700[,iprof_chla],na.fill=T, ends.fill=T, Method="Mean")

### Isolating the CHLA spikes 

	SPIKE_CHLA=CHLA-MED_CHLA

### index of the MLD
	i_mld=which.min(abs(PRES_CTD[,iprof_chla]-MLD))

### index of ipar15
	i_ipar15=which.min(abs(PRES_CTD[,iprof_chla]-ipar_15_depth))

### inndex of i_10
	i_10=which.min(abs(PRES_CTD[,iprof_chla]-10))

### RAPPORT in the MLD with the light limitation
	RAPP_in_MLD=(MED_CHLA[1:i_ipar15,iprof_chla])/(MED_BBP700[1:i_ipar15,iprof_chla])

### Sackmann, max of the Rapp between CHLA and BBP

	i_sack=which.max(RAPP_in_MLD)

#######################################################################
### Estimating the NPQ correction on a smoothed CHLA
#######################################################################

	if (FLAG_SHALLOW & !FLAG_NO_PAR) {

# sigmoid 
#		CHLA_NPQ_D[(i_mld+1):i_ipar15,iprof_chla]=MED_CHLA[(i_mld+1):i_ipar15,iprof_chla]/(0.092+0.908/(1+(MED_PAR_CHLA[(i_mld+1):i_ipar15]/261)^2.2))

		CHLA_NPQ_D[i_mld:i_ipar15,iprof_chla]=MED_CHLA[i_mld:i_ipar15,iprof_chla]/(0.092+0.908/(1+(MED_PAR_CHLA[i_mld:i_ipar15]/261)^2.2))

		RAPP_at_MLD=(CHLA_NPQ_D[i_mld,iprof_chla])/(MED_BBP700[i_mld,iprof_chla])

# shallow than the MLD 
		CHLA_NPQ_D[1:i_mld,iprof_chla]=RAPP_at_MLD*MED_BBP700[1:i_mld,iprof_chla]

	} else {

		CHLA_NPQ_D[1:i_sack,iprof_chla]=MED_BBP700[1:i_sack,iprof_chla]*max(RAPP_in_MLD)

	} # end if FLAG_SHALLOW

######## CATHERINE #############################################################################
######## Double check that we won't miss a max in the MLD 

	if ( length(which(CHLA_NPQ_D[1:i_mld,iprof_chla]<CHLA[1:i_mld,iprof_chla])) > 5 & median(CHLA_NPQ_D[1:i_mld,iprof_chla]-CHLA[1:i_mld,iprof_chla])>0.05 ) {

		CHLA_NPQ_D=CHLA

		CHLA_NPQ_D[1:i_sack,iprof_chla]=MED_BBP700[1:i_sack,iprof_chla]*max(RAPP_in_MLD)

	}
######## CATHERINE ##############################################################################

#########################################################################################
####  Estimating the Q_NPQ (factor of correction to apply to correct of the quenching) till i_npq index of the pressure level
#########################################################################################

	i_npq=max(i_ipar15,i_sack)

	Q_NPQ=CHLA_NPQ_D[1:i_npq,iprof_chla]/MED_CHLA[1:i_npq,iprof_chla]

	Q_NPQ[which(abs(MED_CHLA[1:i_npq,iprof_chla])<0.007)]=NA ### no other idea for that point 

	#### Only make sense when Q_NPQ > 0 and apply a smooth filtering

	Q_NPQ[which(Q_NPQ<=0)]=NA

	#MED_Q_NPQ=RunningFilter(2,Q_NPQ,na.fill=T, ends.fill=T, Method="Median")

#########################################################################################
###	Correct the spikes of the quenching
#########################################################################################

	#SPIKE_CHLA_NPQ=SPIKE_CHLA[1:i_npq,iprof_chla]*MED_Q_NPQ

	#SPIKE_CHLA_NPQ[is.na(SPIKE_CHLA_NPQ)]=0.0

#########################################################################################
### 	Put the spikes back on the CHLA filtered corrected from the Quenching 
#########################################################################################

#	CHLA_NPQ[1:i_npq,iprof_chla]=CHLA_NPQ_D[1:i_npq,iprof_chla]+SPIKE_CHLA_NPQ

	CHLA_NPQ[1:i_npq,iprof_chla]=CHLA_NPQ_D[1:i_npq,iprof_chla]

#	Same value in the first 10m Try_cat (Xing 2018)
	print(CHLA_NPQ_D[i_10,iprof_chla])
	print(median(CHLA_NPQ[1:i_10,iprof_chla]))

 	if(i_sack>i_10)CHLA_NPQ[1:i_10,iprof_chla]=median(CHLA_NPQ_D[1:i_10,iprof_chla])


#########################################
# PLOT
#########################################

	MINDEPTH=max(PRES_CTD[,iprof_chla])
	MAXCHLA=max(CHLA_NPQ,na.rm=TRUE)
	MINDEPTH=MLD+50

	name_file=str_sub(IDnc,str_length(IDnc)-15,str_length(IDnc)-3)


	png(file=paste(name_file,".png",sep=""))

	matplot(CHLA_NPQ[,iprof_chla],PRES_CTD[,iprof_chla],col=2,lwd=2,type="l",ylab="Depth [m]",cex.lab=1.5,cex.axis=1.5,xlab=expression("Chlorophyll a [mg."*m ^ -3 * "]"),xlim=c(-0.2,MAXCHLA+0.5),ylim=rev(c(0, MINDEPTH)))
	matplot(CHLA[,iprof_chla],PRES_CTD[,iprof_chla],col=5,lwd=2,type="l",ylab="Depth [m]",cex.lab=1.5,cex.axis=1.5,xlab=expression("Chlorophyll a [mg."*m ^ -3 * "]"),xlim=c(-0.2,MAXCHLA+0.5),ylim=rev(c(0, MINDEPTH)),add=TRUE)

	legend("bottomright",c("CHLA_NPQ","CHLA"),pch=c(".","."),lwd=c(2,2),col=c(2,5),lty=c(1,1),cex=1.2)


#text(MAXJULD-100,MINGAIN-0.1,paste("mean GAIN m =",round(mean(GAIN_month,na.rm=TRUE),4)," error from sd =",err_month))
#text(MAXJULD-100,MINGAIN-0.15,paste("mean GAIN a =",round(mean(GAIN_ann,na.rm=TRUE),4)," error from sd =",err_ann)) 

	dev.off()

###############################################################
####    changing the NPQ_QC in CHLA_ADJUSTED_QC (to confirm) 
###############################################################

	N_QC=nchar(CHLA_ADJUSTED_QC[iprof_chla])

	QC_test=unlist(strsplit(CHLA_ADJUSTED_QC[iprof_chla],split="")) # to keep 4 set by GL

	if ( length(which(QC_test != "4" | QC_test != " ")) == 0 )  PROFILE_QC = "GL"

	index_qc=which(!is.na(CHLA[,iprof_chla]))

	for (i_qc in seq(1,length(index_qc))) {

		j=index_qc[i_qc] # to avoid _ filled_values 

		if ( j <= i_npq) { 

			substr(CHLA_ADJUSTED_QC[iprof_chla],j,j)<-as.character(5) ### Quenching QC

		} else { 

 			substr(CHLA_ADJUSTED_QC[iprof_chla],j,j)<-as.character(1)

		}

		#### and redo the Range Checking 		

		if ( CHLA_NPQ[j,iprof_chla] < -0.1 | CHLA_NPQ[j,iprof_chla] > 50. | PROFILE_CHLA_QC == "F" ) {

			substr(CHLA_ADJUSTED_QC[iprof_chla],j,j)<-as.character(4)

		}	
	
	}


###########################################
# Write in the working files 
###########################################

	ncvar_put(filenc_B,"CHLA_ADJUSTED",CHLA_NPQ)

	ncvar_put(filenc_B,"CHLA_ADJUSTED_QC",CHLA_ADJUSTED_QC)	

	
	nc_close(filenc_B)

	nc_close(filenc_C)

}# end loop on B files  

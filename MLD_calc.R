MLD_calc <- function ( PRES_CTD, PSAL_CTD , TEMP_CTD ) {

require(oce)

MLD_LIMIT=0.03

THETA=swTheta(PSAL_CTD,TEMP_CTD,PRES_CTD)

POTDENS=swSigmaTheta(PSAL_CTD,TEMP_CTD,PRES_CTD)

FLAG_BAD_POTDENS=rep(FALSE,length(PRES_CTD[,1]))

for(i in 2 : length(PRES_CTD[,1])) {

	TEST=POTDENS[i,1]-POTDENS[i-1,1]

	if(!is.na(TEST)) {

        if( TEST<=-0.03) FLAG_BAD_POTDENS[i-1]=TRUE

	}
}

POTDENS[which(FLAG_BAD_POTDENS==TRUE),1]=NA

PRES_CTD[is.na(POTDENS)]=NA

TEMP_CTD[is.na(POTDENS)]=NA

PSAL_CTD[is.na(POTDENS)]=NA

abs_pres_10=abs(PRES_CTD[,1]-10)

POTDENS_10=max(POTDENS[which(abs_pres_10==min(abs_pres_10,na.rm=TRUE)),1])

# index pour trouver la profondeur de la MLD

MLD_CALC=(POTDENS-POTDENS_10)

if (length(which((MLD_CALC>MLD_LIMIT)==TRUE)) > 0) {

	i_MLD=min(which(MLD_CALC>MLD_LIMIT))

	MLD=PRES_CTD[i_MLD]

} else {

# on initialise au max de profondeur au cas ou on ne trouverait pas de MLD

	MLD=max(PRES_CTD,na.rm=TRUE)

}

return(MLD)

}


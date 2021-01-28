IPAR_15_DEPTH <- function ( PRES_CTD, PAR ) {

PAR_LIMIT=15

if (length(which((PAR < PAR_LIMIT)==TRUE)) > 0) {

	i_par=min(which(PAR<PAR_LIMIT))

	ipar_15_depth=PRES_CTD[i_par]

} else {

	ipar_15_depth=max(PRES_CTD[which(!is.na(PAR))],na.rm=TRUE) ### bof bof 

}

return (ipar_15_depth)

}



!********************************************************************

	subroutine read_weutro_param(param_file)

	implicit none

	character*(*) param_file

	include 'weutro.h'

	real iavpar
	common /iavpar/ iavpar
	save /iavpar/

	namelist /nml_general_parameters/ graztype,NUTLIM,CCHL,CCHLX
	namelist /nml_phyto/ K1RT,K1RC,K1T,K1C,KMPHYT,KPZDC,KPZDT, &
     &		KMNG1,KMPG1,FOPO4,K1G,K1D,ZOO,ZOOSG,KGRZ, &
     &		KPHYZ,EFF,KDZ
	namelist /nml_phosphorus/ PCRB,FOP,K58T,K58C,KOPDC,KOPDT,KPO4
	namelist /nml_nitrogen/ NCRB,FON,K1320C,K1320T,K140C,K140T,KNIT, &
     &		KNO3,K1013C,K1013T,KONDC,KONDT
	namelist /nml_CBOD/ OCRB,KDC,KDT,KDSC,KDST,KBOD
	namelist /nml_dissoxyg/ WIND,AIRTMP,WTYPE,XICECVR,K2,SOD1D,SODTA
	namelist /nml_light_limitation/ LGHTSW,IS2
	namelist /nml_ditoro/ FDAY,IS1,IS1X
	namelist /nml_smith/ PHIMX,XKC,ITOT,iavpar,iav,DTDAY,NEWDAY
	namelist /nml_light_attenuation/ KE,KEFN,KESG

	integer iu,ios

	iu = 1
	open(iu,file=param_file,status='old')

	read(iu,nml=nml_general_parameters,iostat=ios)
	read(iu,nml=nml_phyto,iostat=ios)
	read(iu,nml=nml_phosphorus,iostat=ios)
	read(iu,nml=nml_nitrogen,iostat=ios)
	read(iu,nml=nml_CBOD,iostat=ios)
	read(iu,nml=nml_dissoxyg,iostat=ios)
	read(iu,nml=nml_light_limitation,iostat=ios)
	read(iu,nml=nml_ditoro,iostat=ios)
	read(iu,nml=nml_smith,iostat=ios)
	read(iu,nml=nml_light_attenuation,iostat=ios)

	close(iu)

	CCHLX(1) = CCHL
        IS1X(1) = IS1
        iav = iavpar*itot/FDAY

	end

!********************************************************************

	subroutine write_weutro_param(iunit)

	implicit none

	integer iu,iunit
	include 'weutro.h'

	real iavpar
	common /iavpar/ iavpar
	save /iavpar/

	iu = iunit
	if( iu == 0 ) iu = 6 

	write(iu,*) 'section = general_parameters'
	   write(iu,*) '   graztype   = ',graztype  
	   write(iu,*) '   NUTLIM     = ',NUTLIM    
	   write(iu,*) '   CCHL       = ',CCHL      
	   write(iu,*) '   CCHLX(1)   = ',CCHLX(1)  
	write(iu,*) 'section = phyto'
	   write(iu,*) '   K1RT       = ',K1RT      
	   write(iu,*) '   K1RC       = ',K1RC      
	   write(iu,*) '   K1T        = ',K1T       
	   write(iu,*) '   K1C        = ',K1C       
	   write(iu,*) '   KMPHYT     = ',KMPHYT    
	   write(iu,*) '   KPZDC      = ',KPZDC     
	   write(iu,*) '   KPZDT      = ',KPZDT     
	   write(iu,*) '   KMNG1      = ',KMNG1     
	   write(iu,*) '   KMPG1      = ',KMPG1     
	   write(iu,*) '   FOPO4      = ',FOPO4     
	   write(iu,*) '   K1G        = ',K1G       
	   write(iu,*) '   K1D        = ',K1D       
	   write(iu,*) '   ZOO        = ',ZOO       
	   write(iu,*) '   ZOOSG(1)   = ',ZOOSG(1)  
	   write(iu,*) '   KGRZ       = ',KGRZ      
	   write(iu,*) '   KPHYZ      = ',KPHYZ     
	   write(iu,*) '   EFF        = ',EFF       
	   write(iu,*) '   KDZ        = ',KDZ       
	write(iu,*) 'section = phosphorus'
	   write(iu,*) '   PCRB       = ',PCRB      
	   write(iu,*) '   FOP        = ',FOP       
	   write(iu,*) '   K58T       = ',K58T      
	   write(iu,*) '   K58C       = ',K58C      
	   write(iu,*) '   KOPDC      = ',KOPDC     
	   write(iu,*) '   KOPDT      = ',KOPDT     
	   write(iu,*) '   KPO4       = ',KPO4      
	write(iu,*) 'section = nitrogen'
	   write(iu,*) '   NCRB       = ',NCRB      
	   write(iu,*) '   FON        = ',FON       
	   write(iu,*) '   K1320C     = ',K1320C    
	   write(iu,*) '   K1320T     = ',K1320T    
	   write(iu,*) '   K140C      = ',K140C     
	   write(iu,*) '   K140T      = ',K140T     
	   write(iu,*) '   KNIT       = ',KNIT      
	   write(iu,*) '   KNO3       = ',KNO3      
	   write(iu,*) '   K1013C     = ',K1013C    
	   write(iu,*) '   K1013T     = ',K1013T    
	   write(iu,*) '   KONDC      = ',KONDC     
	   write(iu,*) '   KONDT      = ',KONDT     
	write(iu,*) 'section = CBOD'
	   write(iu,*) '   OCRB       = ',OCRB      
	   write(iu,*) '   KDC        = ',KDC       
	   write(iu,*) '   KDT        = ',KDT       
	   write(iu,*) '   KDSC       = ',KDSC      
	   write(iu,*) '   KDST       = ',KDST      
	   write(iu,*) '   KBOD       = ',KBOD      
	write(iu,*) 'section = dissoxyg'
	   write(iu,*) '   WIND       = ',WIND      
	   write(iu,*) '   AIRTMP     = ',AIRTMP    
	   write(iu,*) '   WTYPE      = ',WTYPE     
	   write(iu,*) '   XICECVR    = ',XICECVR   
	   write(iu,*) '   K2         = ',K2        
	   write(iu,*) '   SOD1D(1)   = ',SOD1D(1)  
	   write(iu,*) '   SODTA(1)   = ',SODTA(1)  
	write(iu,*) 'section = light_limitation'
	   write(iu,*) '   LGHTSW     = ',LGHTSW    
	   write(iu,*) '   IS2        = ',IS2       
	write(iu,*) 'section = ditoro'
	   write(iu,*) '   FDAY       = ',FDAY      
	   write(iu,*) '   IS1        = ',IS1       
	   write(iu,*) '   IS1X(1)    = ',IS1X(1)   
	write(iu,*) 'section = smith'
	   write(iu,*) '   PHIMX      = ',PHIMX     
	   write(iu,*) '   XKC        = ',XKC       
	   write(iu,*) '   ITOT       = ',ITOT      
	   write(iu,*) '   iavpar     = ',iavpar    
	   write(iu,*) '   iav        = ',iav       
	   write(iu,*) '   DTDAY      = ',DTDAY     
	   write(iu,*) '   NEWDAY     = ',NEWDAY    
	write(iu,*) 'section = light_attenuation'
	   write(iu,*) '   KE(1)      = ',KE(1)     
	   write(iu,*) '   KEFN(1)    = ',KEFN(1)   
	   write(iu,*) '   KESG(1)    = ',KESG(1)   

	end

!********************************************************************

!	program read_nml
!	implicit none
!	call read_weutro_param('weutro.nml')
!	call write_weutro_param(0)
!	end

!********************************************************************


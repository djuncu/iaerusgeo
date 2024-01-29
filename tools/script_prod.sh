#! /bin/bash

#------------------------------------------------------------------------------------------
# DEFINITIONS LIEES A L'ENVIRONNEMENT DE TRAVAIL

# ICARE - B. Six
#DATE_EXE=date
## Repertoire de base des donnees d'entree GeoL1B, Angles, AE1 Full disk
#FDSK_BASE_DIR=/DATA/LIENS/GEO/MSG+0000
## Repertoire de base des extractions AERONET des donnees d'entree GeoL1B, Angles, AE1 et AE2 
#ARNT_BASE_DIR=/work_dev/data/six/AERUS/SEV_AERUS.AERONET-2D
## Repertoire de base des produits d'entree AE2
##AER_BASE_DIR=/work_dev/data/six/AERUS
#AER_BASE_DIR=./OUTPUT
## Repertoire de base des produits de sortie AE2
##AER_OUT_DIR=/work_dev/data/six/AERUS
#AER_OUT_DIR=./OUTPUT

# CNRM - X. Ceamanos
DATE_EXE=date
# Repertoire de base des donnees d'entree GeoL1B, Angles, AE1 Full disk
FDSK_BASE_DIR=/cnrm/vegeo/SAT/DATA/AERUS/FULLDISK
# Repertoire de base des extractions AERONET des donnees d'entree GeoL1B, Angles, AE1 et AE2 
ARNT_BASE_DIR=/cnrm/vegeo/SAT/DATA/AERUS/AERONET
# Repertoire de base des produits d'entree AE2
AER_BASE_DIR=outputs
# Repertoire de base des produits de sortie AE2
AER_OUT_DIR=outputs

# MacBook - X. Ceamanos
#DATE_EXE=gdate
## Repertoire de base des donnees d'entree GeoL1B, Angles, AE1 Full disk
#FDSK_BASE_DIR=/Users/xavier.ceamanos/Documents/CNRM/Data/AERUS-GEO/input/FULLDISK
## Repertoire de base des extractions AERONET des donnees d'entree GeoL1B, Angles, AE1 et AE2 
#ARNT_BASE_DIR=/Users/xavier.ceamanos/Documents/CNRM/Data/AERUS-GEO/input/AERONET
## Repertoire de base des produits d'entree AE2
#AER_BASE_DIR=/Users/xavier.ceamanos/Documents/CNRM/Data/AERUS-GEO/output
## Repertoire de base des produits de sortie AE2
#AER_OUT_DIR=/Users/xavier.ceamanos/Documents/CNRM/Data/AERUS-GEO/output

#------------------------------------------------------------------------------------------

function usage {
    echo
    echo "Usage : `basename $0` -d <datedeb> [ -D <datefin> ] [ -m 0|1|2 ] [ -n ] [ -a ] [ -s ] [ -b ] [ -c ] [ -e ] [ -u ] [ -v ] [ -V] [ -t <tag> ] [ -x ] [ -F ]"
    echo " -d: processing beginning date: YYYYMMDD (daily mode) or YYYYMMDDhhmm (instantaneous mode)"
    echo " -D: processing end date (same format as beg. date); beginning date if missing"
    echo " -m: mode - 0:daily (default), 1:instantaneous; 2:daily + instantaneous"
    echo " -n: NRT mode (else reanalysis mode)"
    echo " -a: AERONET subsetting (else full disk mode)"
    echo " -s: starts a new processing period (else try to continue)"
    echo " -b: add browses to L3 product (unactivated if AERONET subsetting)"
    echo " -c: check L1 data disponibility for the day of run and L3 data for the eve"
    echo " -e: no L1 check at the end date"
    echo " -F: use FLOTSAM as RTM"
    echo " -u: use unbiased AOD inputs (default: standard AOD inputs)"
    echo " -v: verif only (do not lauch execution)"
    echo " -V: verbose (displays executed commands); automatic if -v"
    echo " -t: tag corresponding to the code version"
    echo " -x: execute en mode accelere (traitement de VIS06 uniquement)"
    echo
}

function getdate {
    dt=$1; inc=$2; daily=$3
    dd=${dt:0:8}
    if [ $daily -ne 0 ]; then
	$DATE_EXE -d "$dd ${inc}day" +"%Y%m%d" -u
    else
#	$DATE_EXE -d "${dd} ${dt:8} ${inc}hour" +"%Y%m%d%H%M" -u
	$DATE_EXE -d "${dd} ${dt:8} ${inc}minute" +"%Y%m%d%H%M" -u
    fi
}

function process {
    cmd="$1"
    if [ $VERBOSE -ne 0 ]; then echo $cmd; fi
    if [ $DEBUG -eq 0 ]; then eval $cmd; fi
}

START=0; BROWSE=0; DATEDEB=''; DATEFIN=''; ENDCHK=1; DEBUG=0; DAILY=1; NRT=0;
AERONET=0;EXECMODE=0;CHECK=0;VERBOSE=0;TAG='';UNBIAS=0;ACCELERATE=0;FLOTSAM_FLAG=0;
while getopts d:D:m:naxsbceuvVt:Fh c; do
    case $c in
	d) DATEDEB=$OPTARG;;
	D) DATEFIN=$OPTARG;;
	m) EXECMODE=$OPTARG;;
	n) NRT=1;;
	a) AERONET=1;;
	x) ACCELERATE=1;;
	s) START=1;;
	b) BROWSE=1;;
	c) CHECK=1;;
	e) ENDCHK=0;;
	u) UNBIAS=1;;
	v) DEBUG=1;;
	V) VERBOSE=1;;
	t) TAG=$OPTARG;;
        F) FLOTSAM_FLAG=1;;
	h) usage; exit;;
	\?) echo "Option invalide"; usage; exit;;
    esac
done
shift $(($OPTIND - 1))

if [ ${#DATEDEB} -eq 0 ]; then usage; exit; fi
if [ ${#DATEFIN} -eq 0 ]; then DATEFIN=${DATEDEB}; fi

if [ $DEBUG -ne 0 ]; then VERBOSE=1; fi

if [ $EXECMODE -eq 0 ]; then
    DAILY=1; dtinc=1; TRIHOR=0; DATEDEB=`getdate $DATEDEB 0 1`; DATEFIN=`getdate $DATEFIN 0 1`
else
#    DAILY=0; dtinc=3; TRIHOR=1; DATEDEB=`getdate $DATEDEB 0 0`; DATEFIN=`getdate $DATEFIN 21 0`
    DAILY=0; dtinc=15; TRIHOR=0; DATEDEB=`getdate $DATEDEB 0 0`; DATEFIN=`getdate $DATEFIN 0 0`
fi

# Repertoire principal de la chaîne de traitement.
MAIN_DIR=.

BASE_DIR=$FDSK_BASE_DIR

if   [ $NRT -ne 0 ]; then
    STRNRT='-NRT'
    pcf=tools/recursion_template_nrt.pcf
    AE1VRS=1.00; AE2VRS=1.02
else
    AE1VRS=1.03; AE2VRS=1.04
    STRNRT=''
    if [ $AERONET -ne 0 ]; then
	BASE_DIR=$ARNT_BASE_DIR
	#AER_BASE_DIR=$BASE_DIR
	#AER_OUT_DIR=$BASE_DIR
    	pcf=tools/recursion_template_aeronet.pcf 
	STRAERO='.AERONET'
    else
	pcf=tools/recursion_template_rea.pcf 
	STRAERO=''
    fi
fi

if [ $UNBIAS -ne 0 ]; then
    UNB_SFX='_u'
#    AE1_DIR=SEV_AERUS${STRNRT}-L1${STRAERO}.debias.v${AE1VRS}
    AE1_DIR=GEO_AERUS${STRNRT}-L1-MSG+0000${STRAERO}.v${AE1VRS}
    AE2_DIR=aerus-l3.1.1.9.beta_iAERUSGEO_v2.0.7.2_test_1.0${UNB_SFX}
else
    UNB_SFX='';
    AE1_DIR=SEV_AERUS${STRNRT}-L1${STRAERO}.v${AE1VRS}
    AE2_DIR=aerus-l3.1.1.9.beta
fi

PCF_BASE_DIR=$MAIN_DIR/PCF${UNB_SFX}
LOG_BASE_DIR=$MAIN_DIR/LOG${UNB_SFX}
TMP_DIR=$MAIN_DIR/TMP${UNB_SFX}$TAG

AE1_BASE_DIR=$BASE_DIR/$AE1_DIR
#AE2_BASE_DIR=$BASE_DIR/$AE2_DIR
AE2_BASE_DIR=$AER_OUT_DIR
AE2_OUT_DIR=$AER_OUT_DIR/$AE2_DIR

BRW_DIR=$MAIN_DIR/AERUS_Browse

process "mkdir -p $TMP_DIR"

dt=$DATEDEB;dfin=`getdate $DATEFIN $dtinc $DAILY`
jdeb=`getdate $dt 0 1`;ojdt=`getdate $dt -1 1`
verb=2
first_day=1

# Number of slots before if instantaneous mode
#NSLOTS=0  # 1 slot and 0 slots of prior
#NSLOTS=4  # 1 slot and 4 slots of prior (1 hour)
NSLOTS=8  # 1 slot and 8 slots of prior (2 hours)

# Number of history files J-1, J-2, ...
MAXHIST=3

while [ $dt != $dfin ]; do

    echo $dt

    year=${dt:0:4}
    day=${dt:0:4}_${dt:4:2}_${dt:6:2}
    PCF_DIR=$PCF_BASE_DIR$TAG/$year
    LOG_DIR=$LOG_BASE_DIR$TAG/$year
    OUT_DIR=$AE2_OUT_DIR$TAG/$year/$day

    jdt=${dt:0:8}

    if [ $jdt != $jdeb ]; then START=0;fi

    if [ $ojdt != $jdt ] && [ $CHECK -ne 0 ]; then

	if [ $DEBUG -ne 0 ]; then
	    if [ $dt != $DATEFIN ] || [ $ENDCHK -ne 0 ]; then
		echo "L1 current production checking"
	    else
		echo "No L1 current production checking"
	    fi
	    echo "L3 current production checking"
	else

            # Check L1 current production
	    if [ $dt != $DATEFIN ] || [ $ENDCHK -ne 0 ]; then
		dtchk1=`getdate $dt $dtinc $DAILY`
		ychk1=${dtchk1:0:4}; dchk1=${dtchk1:0:4}_${dtchk1:4:2}_${dtchk1:6:2}
		CHK1_DIR=$AE1_BASE_DIR/$ychk1/$dchk1
		chk1=`ls $CHK1_DIR/* &> /dev/null`
		if [ $? -ne 0 ]; then echo "No data in $CHK1_DIR. Sleeping 5mn"; sleep 300; continue; fi
	    fi

            # Check L3 current production
	    dtchk3=`getdate $dt -1 1`
	    ychk3=${dtchk3:0:4}; dchk3=${dtchk3:0:4}_${dtchk3:4:2}_${dtchk3:6:2}
	    CHK3_DIR=$AE2_BASE_DIR/$ychk3/$dchk3
	    chk3=`ls $CHK3_DIR/* &> /dev/null`
	    if [ $? -ne 0 ]; then
		echo $dt $dchk3
		if [ $START -eq 0 ]; then echo "No data in $CHK3_DIR. Cannot continue processing"; exit; fi
	    else
		if [ $START -ne 0 ]; then echo "$CHK3_DIR contains data. Cannot start"; exit; fi
	    fi
	fi
    fi

    if [ $ojdt != $jdt ] && [ $EXECMODE -eq 2 ]; then
	rdt=$jdt; trh=0
    else
	rdt=$dt ; trh=$TRIHOR
    fi

    LOGFIL=${LOG_DIR}/${rdt}.log
    PCFFIL=${PCF_DIR}/${rdt}.pcf

    process "mkdir -p $PCF_DIR $OUT_DIR $LOG_DIR"
    process "python tools/MakeAerusL3Pcf.py -v $AE2VRS -M $MAIN_DIR -O $OUT_DIR -T $TMP_DIR -B $AER_BASE_DIR -A $BASE_DIR -1 $AE1_DIR -2 $AE2_DIR -t $trh -s $START -x $ACCELERATE -d $rdt -V $verb -n $NRT -D 3 -N $NSLOTS -m $MAXHIST -F $FLOTSAM_FLAG $pcf &> $PCFFIL"

#    process "(time python src/FMK_AerusL3.py -Y 360 -g MSG+0000 $PCFFIL) &> $LOGFIL"
#    process "(time python src/FMK_AerusL3.py -Y 180,210,360,515,578 -g MSG+0000 $PCFFIL) &> $LOGFIL"
#    process "(time python src/FMK_AerusL3.py -Y 88,141,180,204,210,286,360,482,515,578 -g MSG+0000 $PCFFIL) &> $LOGFIL"
#    process "(time python src/FMK_AerusL3.py -Y 85,141,180,204,210,286,330,360,479,482,515,560,574,578 -g MSG+0000 $PCFFIL) &> $LOGFIL"
#    process "(time python src/FMK_AerusL3.py -Y 180,286,578 -g MSG+0000 $PCFFIL) &> $LOGFIL"
#    process "(time python src/FMK_AerusL3.py -Y 85,87,362 -g MSG+0000 $PCFFIL) &> $LOGFIL"
#    process "(time python src/FMK_AerusL3.py -Y 141,180,204,210,286,360,479,482,515,560,578 -g MSG+0000 $PCFFIL) &> $LOGFIL"
#    process "(time python src/FMK_AerusL3.py -Y 1,96,111,180,217,360,515,578,583 -g MSG+0000 $PCFFIL) &> $LOGFIL"
    #process "(time python src/FMK_AerusL3.py -g MSG+0000 $PCFFIL) &> $LOGFIL"
    #process "(time python src/FMK_AerusL3_flotsam.py -g MSG+0000 $PCFFIL) &> $LOGFIL"
    process "(time python src/FMK_AerusL3.py -g MSG+0000 $PCFFIL)"

    process "rm -f $TMP_DIR/*"

    if [ $BROWSE -ne 0 ] && [ $AERONET -eq 0 ]; then
	if [ $DEBUG -eq 0 ]; then echo "--------------------------------------------------------------------------------" >> $LOGFIL 2>&1; fi
	for f in $OUT_DIR/SEV_AERUS-AEROSOL*.h5 $OUT_DIR/SEV_AERUS-ALBEDO*.h5; do
	    process "(time python $BRW_DIR/src/Browse_AerusL3.py -v DEFAULT -O $OUT_DIR -M $BRW_DIR -i $f) >> $LOGFIL 2>&1"
	done
    fi

    if [ $EXECMODE -ne 2 ] || [ $ojdt == $jdt ] ; then dt=`getdate $dt $dtinc $DAILY`; fi
    ojdt=$jdt

done

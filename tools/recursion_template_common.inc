# Region geostationnaire a traiter (MSG+0000 ou MSG+0415).
GEOREG = $GEOREG

# Repertoires de base des produits AERUS-GEO
AE1BASEDIR = $GEOBASEDIR
AE2BASEDIR = $GEOBASEDIR

# Repertoire de base des produits GEO
GEOBASEDIR = $GEOBASEDIR

# Sous-repertoire des angles
ANGLEDIR = GEO_L1B-ANGLES.v1.00

# Repertoires de base des produits AERUS-GEO
#$AE1BASEDIR = $GEOBASEDIR
#$AE2BASEDIR = $GEOBASEDIR

# Sous-repertoire des produits AERUS-L1
AE1DIR = GEO_AERUS${NRTEXT}-L1-${GEOREG}.v${VERSL1}

# Sous-repertoire des produits "memoire" AERUS-L3
AE2DIR = SEV_AERUS-INTERNAL${NRTEXT}-D3.v${VERSL3}

# Version produit (DEFAULT -> gere par le code)
PROD_VER=$PROD_VER

# Identifiant Icare - pour production op. uniquement
ICARE_ID=RFU

# Recursion?
RECURSION = 1

# Composition?
COMPO = 0

# Premier jour de la periode de production?
START = 0

# NRT ?
NRT = $NRT

# Structure des arborescences AE1 et AE2
# Niveau "year" et "date" (YYYY_MM_DD) dans l'arborescence 
#     |  0  |  1  |  2  |  3  |
#-----+-----+-----+-----+-----+
# AE1 | NON | OUI | NON | OUI |
# AE2 | NON | NON | OUI | OUI |

YEARDIR = 3
DATEDIR = 1

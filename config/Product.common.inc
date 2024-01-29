# NAM_ALG_MODES
MODE        = 1
OCEAN_OK    = True
LAND_OK     = True
COAST_OK    = True

# Auxilliary input files pathnames
LAT         = ${REGION}.3km.lat
LON         = ${REGION}.3km.lon
CST_PIX     = ${REGION}.3km.cst.hdf

# NAM_STOP_STATUS_RETRIEVAL_FREQ
FREQINSECS  = 5
 
# Production mode : to produce outputs as only 2 (or 3) files if True
PROD_MODE   = True

# To output Albedos as an official product if True
ALB_OFF_PRD = True

# References
@References.cfg

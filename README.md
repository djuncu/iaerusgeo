CHAINE AERUS-L3
===============


DESCRIPTIF
----------

Cette chaîne constitute le second module AERUS. Elle a pour but la génération de cartes
journalière d'aérosols et d'albédos; elle utilise un processus à mémoire. Dans cette version,
seul MSG+0000 est supporté (MSG 1, 2 ou 3).

Deux flux de traitement sont prévus: un mode "NRT" exploitant les forecasts ECMWF et un mode
retraitement exploitant les réanalyses ERA-Interim.

Cette version unifie les différents modes de traitement posibles : daily, tri-horaire ou
instantané, avec choix de la région stationnaire à considérer. Le traitement multi-géos n'est
pas supporté dans cette version.

- Nom de l'exécutable: FMK_AerusL3.py
- Version: 1.1.9 - 12/06/2020
- Type d'exécutable: programme Python exécutant du code C et Fortran 90 compilé.
- Nombre de fichiers utilisés:
  * 288 fichiers produits de la chaîne AERUS-L1 (au maximum)
  *  96 fichiers d'angles géostationnaires (au maximum)
  *   1 fichier produit AERUS-L3 de la veille (sauf le premier jour de traitement)
  *   3 fichiers de sortie (AEROSOL, ALBEDO, INTERNAL)
  *   4 fichiers de données ancillaires statiques (paramètres techniques, latitudes)
- Temps d'exécution: de l'ordre de 50' CPU par jour
- Mémoire utilisée: environ 2 Go pour MSG
- Taille des fichiers de sortie: environ 600 Mo par jour


BIBLIOTHÈQUES/EXÉCUTABLES REQUIS
--------------------------------

La disponibilité de python (>= 2.6.4, < 3), numpy (>= 1.4), pyhdf (>= 0.8.3), h5py
(>= 1.3.0) et f2py est requise.


COMPILATION
-----------

Dans le répertoire principal de la chaîne :

    make all

tout simplement.

La compilation crée dans le répertoire src une bibliothèque dynamique AerusL3.so utilisée par
le code python.


EXÉCUTION
---------

Le code de la chaîne est FMK_AerusL3.py. La commande d'exécution est:

    python src/FMK_AerusL3.py [-options] [cmdfile]

les paramètres d'exécution pouvant être indifféremment passés en option ou dans un fichier de
commande.

L'exécution en mode DAILY produit dans le répertoire de sortie 3 fichiers format HDF5, de noms:

- SEV_AERUS-AEROSOL[-NRT]-D3_<DATE>_<VERSION>.h5  : Produit officiel public
- SEV_AERUS-ALBEDO[-NRT]-D3_<DATE>_<VERSION>.h5   : Produit officiel accès restreint
- SEV_AERUS-INTERNAL[-NRT]-D3_<DATE>_<VERSION>.h5 : Produit pour récursion uniquement

avec:

    <DATE>   : date du jour de traitement (YYYY-mm-dd)
    <VERSION>: version du produit
    [-NRT]   : le composant NRT n'est présent que dans le mode NRT.

Suivant le mode de traitement, les sorties peuvent avois des noms différents.


STRUCTURE DU FICHIER DE COMMANDE
--------------------------------

    # Identifiant ICARE. Usage interne (inutilisé ici car sortie non HDF).
    # Facultatif. Défaut: "RFU". Option -i.
    ICARE_ID =
    
    # Version du produit de sortie. Bien que requis (comme demandé par l'exploitation), ce
    # paramètre peut prendre la valeur DEFAULT, auquel cas la version produit sera celle
    # définie dans le code (les 2 premiers composants de la version du code).
    # Requis. Option -v.
    PROD_VER =
    
    # Chemin du répertoire principal de la chaîne.
    # Facultatif. Défaut: ".". Option -M.
    MAINDIR =
    
    # Chemin du répertoire de sortie du produit de la chaîne.
    # Facultatif. Défaut: ".". Option -O.
    OUTDIR =
    
    # Chemin du répertoire de travail temporaire (inutilisé mais imposé en production)
    # Facultatif. Défaut: "/tmp". Option -T.
    TMPDIR =
    
    # Niveau des affichages (0: pas d'affichage, 1: affichage restreint; >1: affichage complet)
    # Facultatif. Défaut: 1. Option -V.
    VERBOSE = 
    
    # Date de traitement.
    # Format "YYYYMMDD" en mode daily, format "YYYYMMDDhhmm" en modes tri-horaire ou instantané.
    # Requis. Option -d.
    DATE =

    # Région géostationnaire à traiter (MSG+0000 ou MSG+0415 supportés)
    # Requis. Option -g.
    GEOREG =

    # Exécution en mode quasi temps réel ou non (0: non, 1: oui)
    # Requis. Option -n.
    NRT = 

    # 1er jour de traitement? 0: NON, <>0: OUI.
    # Facultatif. Défaut: 0. Option -s
    START =
    
    # Run daily (1) ou non (0 -> tri-horaire ou instantané).
    # Facultatif. Défaut: 1. Option -D.
    DAILY_MODE =

    # Run tri-horaire (1) ou instantané (0)?.
    # Facultatif. Défaut: 0. Option -t.
    TRIHOR =

    # Fichiers "BRDF Model et covariance matrix" du jour précédant la date du Run.
    # Produit "mémoire" SEV_AERUS-INTERNAL.
    # Requis sauf si START <> 0. Option -c.
    CK_K012_IN =
    
    # Liste des fichiers "BRDF Model et covariance matrix" du(des) jour(s) précédant la veille
    # de la date du Run. Produit "mémoire" SEV_AERUS-INTERNAL.
    # Le nombre de jours dépend du type de run: 30 en 1ère passe NRT, 0 dans les autres modes.
    # Séparateur: espace. # Facultatif. Défaut: "". Option -C.
    ALT_CK_K012_IN =
    
    # Nombre maximal de jours toléré entre le produit "mémoire" et le jour du Run.
    # Dépend du type de run: > 1 (variable) en 1ère passe NRT, 1 dans les autres modes.
    # Facultatif. Défaut: 1. Option -m.
    MAXHIST =
    
    # Information à priori à prendre dans les daily plutôt que dans les 3h. Sans objet en run
    # daily ou instantané.
    # Facultatif. Défaut: 1. Option -h.
    HIST_DAY =
    
    # Liste des fichiers des angles géostationnaires à traiter (produit GeoL1B-Angles)
    # Séparateur: espace. Requis. Option -a.
    ANGLES = 
    
    # Liste des fichiers produits AERUS-L1 à traiter pour la bande VIS06
    # Séparateur: espace. Requis. Option -6.
    VIS06 =  
    
    # Liste des fichiers produits AERUS-L1 à traiter pour la bande VIS08
    # Séparateur: espace. Requis. Option -8.
    VIS08 =  
    
    # Liste des fichiers produits AERUS-L1 à traiter pour la bande IR016
    # Séparateur: espace. Requis. Option -1.
    IR016 =  
    
    # Nom du fichier HDF modèle de subsetting (doit exister dans le répertoire ancillary).
    # Facultatif. Défaut: "" (pas de subsetting). Option -S.
    SUBSET =
    
    # Extension à ajouter en fin de nom de fichier (utile par exemple en cas de subsetting).
    # Facultatif. Défaut: "" (pas d'extension). Option -e.
    EXTENSION =
    	
    # Sélection d'un sous-ensemble d'indices de colonnes du subset (typiquement, par exemple,
    # certaines stations AERONET sélectionnées dans le modèle de subsetting AERONET).
    # Format: i1,i2,..,in, sans espace. Sans effet si SUBSET = "".
    # Séparateur: ','. Facultatif. Défaut: "". Option -Y.
    SUBSUBSET =

    # Transposition du modèle de subsetting (<>0: transposition, 0: pas de transposition)
    # Sans effet si SUBSET = "".
    # Facultatif. Défaut: 0. Option -X.
    SUBTRANSP =
    
    # Sorties des fichiers instantanés de contrôle pour chaque slot.
    # Facultatif. Défaut: 0. Option -I.
    INST_OUT =


OPTIONS D'EXÉCUTION
-------------------

Tous les paramètres du fichier de commande peuvent être passés indifféremment en option. On
peut avoir seulement un fichier de commande, seulement des options ou bien les deux. Si un même
paramètre est passé à la fois dans le fichier de commande et en option, c'est celui passé en
option qui sera pris en compte.


TESTS:
-----

Dans le répertoire test, 3 répertoires:
- pcf   : contient les PCFs de test;
- input : contient les entrées nécessaires;
- output: contient les produits de sortie associés aux PCFs et les logs d'exécution.

Par ailleurs, le répertoire tools contient un générateur de PCF avec quelques templates.


GÉNÉRATION DES FICHIERS DE COMMANDE
-----------------------------------

Le répertoire tools contient un code python MakeAerusL3Pcf.py permettant de générer les
fichiers de commandes pour FMK_AerusL3.py, et en particulier la création des liste de fichiers
de données en fonction de la date et du mode de traitement.

Pour plus de détails sur les paramètres de ce code, voir par exemple le fichier de paramètres
générique recursion_template_common.inc inclus dans tous les autres.


RÈGLES DE PRODUCTION:
--------------------

Contrairement aux versions précédentes, le premier run NRT aura aussi un paramètre START=0
(sauf indication contraire spécifique). Le paramètre CK_K012_IN ne pouvant alors être
renseigné, un fichier de départ (fakestart) est fourni dans le répertoire ancillary.  Cela
permet une grande amélioration des produits pendant la période de démarrage (une vingtaine de
jours).

En mode NRT (paramètre NRT = 1), on démarre le 01/04/2014 et on traite au fil de l'eau dès
que les différentes entrées sont disponibles:
- 1ère passe avec contrainte relâchée sur le fichier "mémoire" (MAXHIST > 1)
- 2ème passe différée avec contrainte stricte sur le fichier "mémoire" (MAXHIST = 1)

En mode retraitement (Réanalyses ERA-Interim, paramètre NRT = 0), on remonte aussi loin que
possible dans le passé pour toute la période pour laquelle ERA-Interim est disponible.

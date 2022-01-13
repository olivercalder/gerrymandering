# Shapefile info:
SHAPEFILE_PATH = '../data/mn_mggg/MN16/MN_precincts16.shp'  # Path to shapefile
POPULATION_COL = 'TOTPOP'       # Column where total population is stored
VAP_COL = 'VAP'                 # Column where voting age population is stored
BVAP_COL = 'BVAP'               # Column where black voting age population is stored
HVAP_COL = 'HVAP'               # Column where hispanic voting age population is stored
DISTRICT_COL = 'CONGDIST'       # Column where district assignment is stored
COUNTY_COL = 'COUNTYNAME'       # Column where county is stored
PAD_DISTRICT_NUMBERS = False    # Pad district numbers with leading 0s

# Output info:
TOTAL_STEPS = 100               # Total number of accepted maps to generate
OUTPUT_PATH = 'MN_results.csv'  # File path for output CSV

# Election info:
ELECTION_NAME = 'PRES16'        # Name of the election used for partisan metrics
ELECTION_DEM_COL = 'PRES16D'    # Name of the column where democratic votes are stored
ELECTION_REP_COL = 'PRES16R'    # Name of the column where republican votes are stored

# Districts:
TOTAL_DISTRICTS = 8             # Total number of districts in the given shapefile

# VRA Compliance:
OPPORTUNITY_THRESHOLD = 0.3     # Percent minority population to be considered "minority opportunity"
BLACK_OPP_TARGET = 2            # Target number of black opportunity districts
BLACK_OPP_WEIGHT = 1            # The weight of black opportunity district VRA compliance
HISPANIC_OPP_TARGET = 0         # Target number of hispanic opportunity districts
HISPANIC_OPP_WEIGHT = 0         # The weight of hispanic opportunity district VRA compliance

# County Split Scoring:
COUNTY_SPLIT_COEFFICIENT = 10   # (this)^k is the coefficient of the w_{2+k} score

# Score function weights:
POPULATION_SCORE_WEIGHT = 1     # Weight for population score
COMPACTNESS_SCORE_WEIGHT = 1    # Weight for compactness score
COUNTY_SCORE_WEIGHT = 1         # Weight for county split score
VRA_SCORE_WEIGHT = 1            # Weight for VRA score

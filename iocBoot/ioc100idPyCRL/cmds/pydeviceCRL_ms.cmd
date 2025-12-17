################################################################################
# PyDevice Support for Transfocator
################################################################################

################################################################################
# Local Definitions

epicsEnvSet("PY_OBJECT", "msCRL")
epicsEnvSet("SYS_ID", "MS")
epicsEnvSet("STACK_SUBS_FILE", "substitutions/pyDevCRL_stacks_ms.substitutions")
epicsEnvSet("SYSTEM_SUBS_FILE", "substitutions/pyDevCRL_system_ms.substitutions")
epicsEnvSet("TOML_FILE", "toml/crl_setup_ms.toml")

# Set CRL numbers
epicsEnvSet("_STACKS1","10")  # Number of stacks
epicsEnvSet("_STACKS2","10")  # Number of stacks
epicsEnvSet("_CONFIGS","1024")  # Possible configurations: 2^(stacks1) 

# Creating beam energy PVs for testing
epicsEnvSet("MONOE","testMonoE") # for testing -- replace with real mono energy PV 
epicsEnvSet("IDENERGY","testIDE") # for testing -- replace with real ID energy PV 

# Setting Mono energy PV
#epicsEnvSet("BLE","$(PREFIX)$(MONOE)")	# Beam energy PV at CRL (testing uses MONOE defined earlier)
epicsEnvSet("BLE","100idPySBL:BraggERdbkAO")

################################################################################
# Load DBs and python code

# Testing tools for energy
dbLoadRecords("${TOP}/db/energyTestTools.db","P=$(PREFIX), MONOE=$(MONOE), IDENERGY=$(IDENERGY)")

# CRL DBs and defining substitution file to get stack properties
dbLoadTemplate("$(STACK_SUBS_FILE)","P=$(PREFIX),SYSID=$(SYS_ID)")
pydev("stack_subFile = '$(STACK_SUBS_FILE)'")

# Add elements (including simulated/real motors for z-position)
dbLoadTemplate("$(SYSTEM_SUBS_FILE)","P=$(PREFIX),SYSID=$(SYS_ID),OBJ=$(PY_OBJECT), ELEM=$(_CONFIGS)")

# Import Transfocator class
pydev("from pyCRL_system import focusingSystem")

# Create Transfocator object
pydev("$(PY_OBJECT) = focusingSystem(crl_setup = '$(TOML_FILE)')")

# DB file for system controls and sample position(s)
dbLoadRecords("${TOP}/db/pyDevCRL_general.db","P=$(PREFIX), SYSID=$(SYS_ID), OBJ=$(PY_OBJECT), KEV=$(BLE), ELEM=$(_CONFIGS)")
dbLoadRecords("${TOP}/db/pyDevCRL_2systems.db","P=$(PREFIX), SYSID=$(SYS_ID), OBJ=$(PY_OBJECT), SYSA=B, SYSB=C")
# For 28 (3 CRLs):
# dbLoadRecords("${TOP}/db/pyDevCRL_3systems.db","P=$(PREFIX), SYSID=$(SYS_ID), OBJ=$(PY_OBJECT), SYSA=B, SYSB=C, SYSC=D")

################################################################################
# Initial setting of some PVs

# Verboseness turned on for debugging -- set to 0 to turn off
doAfterIocInit("dbpf('$(PREFIX)$(SYS_ID):verbosity','1')")

# After iocInit, need to setup focal size lookup table -- because setting slit
# sizes and energy causes the table to be recalculated, will first disable 
# lookup table calc PV, set the slit sizes and energy, calculate the table 
# directly via python call, then enable calc PV so future slit/energy changes
# update the lookup table

# Disable calc of lookup table
doAfterIocInit("dbpf('$(PREFIX)$(SYS_ID):recalc_enable','0')") # Disable table calc

# TODO can 1, 2 be replaced programmatically?
# Process current slit and energy PV settings
doAfterIocInit("dbpf('$(PREFIX)$(SYS_ID):1:slitSize_H.PROC','1')") # PROC H slit size read and update CRL1 object
doAfterIocInit("dbpf('$(PREFIX)$(SYS_ID):1:slitSize_V.PROC','1')") # PROC V slit size read and update CRL1 object
doAfterIocInit("dbpf('$(PREFIX)$(SYS_ID):2:slitSize_H.PROC','1')") # PROC H slit size read and update CRL1 object
doAfterIocInit("dbpf('$(PREFIX)$(SYS_ID):2:slitSize_V.PROC','1')") # PROC V slit size read and update CRL1 object
doAfterIocInit("dbpf('$(PREFIX)$(SYS_ID):energy.PROC','1')") # PROC energy read and update CRL1 object

# Check initial energy and slit settings
doAfterIocInit("pydev('print($(PY_OBJECT).energy)')")
doAfterIocInit("pydev('print($(PY_OBJECT).slits)')")

# Calc lookup table via python; later re-calculations started via EPICS (changes to energy, slits)
# 	stack_subFile - string holding name/rel. location of substitutions file loaded
#                     loaded in  dbLoadTemplate
doAfterIocInit("pydev('$(PY_OBJECT).setupLookupTable(stack_subFile)')") 

# Enable lookup table calc PV
doAfterIocInit("dbpf('$(PREFIX)$(SYS_ID):recalc_enable','1')") # Enable table calc
################################################################################

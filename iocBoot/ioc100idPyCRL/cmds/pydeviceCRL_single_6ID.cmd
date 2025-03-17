################################################################################
# PyDevice Support for Transfocator
################################################################################

################################################################################
# Local Definitions

epicsEnvSet("PY_OBJECT", "CRL6ID")
epicsEnvSet("SYS_ID", "CRL6ID")
epicsEnvSet("SUBS_FILE", "substitutions/pyDevCRL_stacks_single_6ID.substitutions")
epicsEnvSet("TOML_FILE", "toml/crl_setup_single_6ID.toml")

# Set CRL numbers
epicsEnvSet("_STACKS1","9")  # Number of stacks
epicsEnvSet("_CONFIGS","512")  # Possible configurations: 2^(min(stacks1, stacks2)) 

# Creating beam energy PVs for testing
epicsEnvSet("MONOE","testMonoE") # for testing -- replace with real mono energy PV 
epicsEnvSet("IDENERGY","testIDE") # for testing -- replace with real ID energy PV 

# Setting slit PVs
epicsEnvSet('_SLIT1H',"$(PREFIX)testSSH1")	# Horizontal size of slit PV before CRL 1 (testing)
epicsEnvSet('_SLIT1V',"$(PREFIX)testSSV1")	# Vertical size of slit PV before CRL 1 (testing)
#epicsEnvSet('_SLIT1H',"")	# Horizontal size of slit PV before CRL 1
#epicsEnvSet('_SLIT1V',"")	# Vertical size of slit PV before CRL 1

# Setting CRL Z-translation PVs (if z-translation part of CRL system)
#epicsEnvSet('_OEPOS1',"$(PREFIX)testCRL1z")	# Z-motion of CRL 1 (testing)
epicsEnvSet('_OEPOS1',"6idbSoft:TRANS:m6")	# Z-motion of CRL 1

# Setting Sample Z-translation PVs
epicsEnvSet('_SAMPOS',"$(PREFIX)testSAMz")	# Z-motion of sample (testing or unneeded)
#epicsEnvSet('_SAMPOS',"")	# Z-motion of sample

# Setting Mono energy PV
epicsEnvSet("BLE","$(PREFIX)$(MONOE)")	# Beam energy PV at CRL (testing uses MONOE defined earlier)
#epicsEnvSet("BLE","")

################################################################################
# Load DBs and python code

# Next several lines set up some testing tools for energy, slits and z-positions
dbLoadRecords("${TOP}/db/energyTestTools.db","P=$(PREFIX), MONOE=$(MONOE), IDENERGY=$(IDENERGY)")
dbLoadRecords("${TOP}/db/slitTestTools.db","P=$(PREFIX), SLITH=testSSH1, SLITV=testSSV1")
dbLoadRecords("${TOP}/db/posTestTools.db","P=$(PREFIX), ZELEM=testCRL1z")
dbLoadRecords("${TOP}/db/posTestTools.db","P=$(PREFIX), ZELEM=testSAMz")

# CRL DBs and defining substitution file to get stack properties
dbLoadTemplate("$(SUBS_FILE)","P=$(PREFIX),SYSID=$(SYS_ID)")
pydev("stack_subFile = '$(SUBS_FILE)'")

# Add elements
dbLoadRecords("${TOP}/db/pyDevCRL_elem.db","P=$(PREFIX),SYSID=$(SYS_ID),OBJ=$(PY_OBJECT),OE=1,ELEM=$(_CONFIGS),OEPOS=$(_OEPOS1)")
#dbLoadRecords("${TOP}/db/pyDevCRL_elem.db","P=$(PREFIX),SYSID=$(SYS_ID),OBJ=$(PY_OBJECT),OE=2,ELEM=$(_CONFIGS),OEPOS=$(_OEPOS2)")

# Add slits for each element
# Transforcators are numbered (1, 2); KB has string identifier ('kb') for OE ID
dbLoadRecords("${TOP}/db/pyDevCRL_slits.db","P=$(PREFIX),SYSID=$(SYS_ID),OBJ=$(PY_OBJECT),OE=1,SLITH=$(_SLIT1H),SLITV=$(_SLIT1V)")
#dbLoadRecords("${TOP}/db/pyDevCRL_slits.db","P=$(PREFIX),SYSID=$(SYS_ID),OBJ=$(PY_OBJECT),OE=2,SLITH=$(_SLIT2H),SLITV=$(_SLIT2V)")
#dbLoadRecords("${TOP}/db/pyDevCRL_slits.db","P=$(PREFIX),SYSID=$(SYS_ID),OBJ=$(PY_OBJECT),OE=kb,SLITH=$(_SLITKBH),SLITV=$(_SLITKBV)")

# Import Transfocator class
pydev("from pyCRL_system import focusingSystem")

# Create Transfocator object
pydev("$(PY_OBJECT) = focusingSystem(crl_setup = '$(TOML_FILE)')")

# DB file for system controls
dbLoadRecords("${TOP}/db/pyDevCRL_general.db","P=$(PREFIX), SYSID=$(SYS_ID), OBJ=$(PY_OBJECT), KEV=$(BLE), ELEM=$(_CONFIGS),SAMPOS=$(_SAMPOS)")

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

# Process current slit and energy PV settings
doAfterIocInit("dbpf('$(PREFIX)$(SYS_ID):1:slitSize_H.PROC','1')") # PROC H slit size read and update CRL1 object
doAfterIocInit("dbpf('$(PREFIX)$(SYS_ID):1:slitSize_V.PROC','1')") # PROC V slit size read and update CRL1 object
doAfterIocInit("dbpf('$(PREFIX)$(SYS_ID):energy.PROC','1')") # PROC energy read and update CRL1 object

# Check initial energy and slit settings
doAfterIocInit("pydev('print($(PY_OBJECT).energy)')")
doAfterIocInit("pydev('print($(PY_OBJECT).slits)')")

# Calc lookup table via python; later re-calculations started via EPICS (changes to energy, slits)
# 	stack_subFile - string holding name/rel. location of substitutions file loaded
#                     loaded in  dbLoadTemplate
#   $(_STACKS1)  - number of stacks of lenses in first transfocator
#   $(_STACKS2)  - number of stacks of lenses in second transfocator
doAfterIocInit("pydev('$(PY_OBJECT).setupLookupTable(stack_subFile, [$(_STACKS1)])')") 

# Enable lookup table calc PV
doAfterIocInit("dbpf('$(PREFIX)$(SYS_ID):recalc_enable','1')") # Enable table calc
################################################################################

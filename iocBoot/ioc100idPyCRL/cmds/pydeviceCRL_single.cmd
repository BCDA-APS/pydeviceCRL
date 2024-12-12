################################################################################
# PyDevice Support for Transfocator
################################################################################
epicsEnvSet("PY_OBJECT", "CRL")
epicsEnvSet("SYS_ID", "CRL")

# Set CRL numbers
epicsEnvSet("_STACKS1","10")  # Number of stacks
epicsEnvSet("_CONFIGS","1024")  # Possible configurations: 2^(min(stacks1, stacks2)) 

epicsEnvSet("MONOE","testMonoE") # for testing -- replace with real mono energy PV 
epicsEnvSet("IDENERGY","testIDE") # for testing -- replace with real ID energy PV 
# Next two lines set up some testing tools for energy and slits
dbLoadRecords("${TOP}/db/energyTestTools.db","P=$(PREFIX), MONOE=$(MONOE), IDENERGY=$(IDENERGY)")
dbLoadRecords("${TOP}/db/slitTestTools.db","P=$(PREFIX), SLITH=testSSH1, SLITV=testSSV1")

# Setting slit PVs
epicsEnvSet('_SLIT1H',"$(PREFIX)testSSH1")	# Horizontal size of slit PV before CRL 1 (testing)
epicsEnvSet('_SLIT1V',"$(PREFIX)testSSV1")	# Vertical size of slit PV before CRL 1 (testing)
#epicsEnvSet('_SLIT1H',"")	# Horizontal size of slit PV before CRL 1
#epicsEnvSet('_SLIT1V',"")	# Vertical size of slit PV before CRL 1

#epicsEnvSet('_SLIT2H',"$(PREFIX)testSSH2")	# Horizontal size of slit PV before CRL 2 (testing)
#epicsEnvSet('_SLIT2V',"$(PREFIX)testSSV2")	# Vertical size of slit PV before CRL 2 (testing)
#epicsEnvSet('_SLIT2H',"")	# Horizontal size of slit PV before CRL 2
#epicsEnvSet('_SLIT2V',"")	# Vertical size of slit PV before CRL 2

#epicsEnvSet('_SLITKBH',"$(PREFIX)testSSHKB")	# Horizontal size of slit PV before KB (testing)
#epicsEnvSet('_SLITKBV',"$(PREFIX)testSSVKB")	# Vertical size of slit PV before KB (testing)
#epicsEnvSet('_SLITKBH',"")	# Horizontal size of slit PV before KB
#epicsEnvSet('_SLITKBV',"")	# Vertical size of slit PV before KB


# Setting Mono energy PV
epicsEnvSet("BLE","$(PREFIX)$(MONOE)")	# Beam energy PV at CRL (testing)


# CRL DBs and defining substitution file to get stack properties
dbLoadTemplate("substitutions/pyDevCRL_stacks_single.substitutions","P=$(PREFIX),SYSID=$(SYS_ID)")
pydev("stack_subFile = 'substitutions/pyDevCRL_stacks_single.substitutions'")

# Add elements
dbLoadRecords("${TOP}/db/pyDevCRL_elem.db","P=$(PREFIX),SYSID=$(SYS_ID),OBJ=$(PY_OBJECT),OE=1,ELEM=$(_CONFIGS)")
#dbLoadRecords("${TOP}/db/pyDevCRL_elem.db","P=$(PREFIX),SYSID=$(SYS_ID),OBJ=$(PY_OBJECT),OE=2,ELEM=$(_CONFIGS)")

# Add slits for each element
# Transforcators are numbered (1, 2); KB has string identifier ('kb') for OE ID
dbLoadRecords("${TOP}/db/pyDevCRL_slits.db","P=$(PREFIX),SYSID=$(SYS_ID),OBJ=$(PY_OBJECT),OE=1,SLITH=$(_SLIT1H),SLITV=$(_SLIT1V)")
#dbLoadRecords("${TOP}/db/pyDevCRL_slits.db","P=$(PREFIX),SYSID=$(SYS_ID),OBJ=$(PY_OBJECT),OE=2,SLITH=$(_SLIT2H),SLITV=$(_SLIT2V)")
#dbLoadRecords("${TOP}/db/pyDevCRL_slits.db","P=$(PREFIX),SYSID=$(SYS_ID),OBJ=$(PY_OBJECT),OE=kb,SLITH=$(_SLITKBH),SLITV=$(_SLITKBV)")

# Import Transfocator class
pydev("from pyCRL_system import focusingSystem")

# Create Transfocator object
pydev("$(PY_OBJECT) = focusingSystem(crl_setup = 'crl_setup_single.toml')")

# DB file for system controls
dbLoadRecords("${TOP}/db/pyDevCRL_general.db","P=$(PREFIX), SYSID=$(SYS_ID), OBJ=$(PY_OBJECT), KEV=$(BLE), ELEM=$(_CONFIGS)")

################################################################################
# Verboseness turned on for debugging -- set to 0 to turn off
doAfterIocInit("dbpf('$(PREFIX)$(SYS_ID):verbosity','1')")

# After iocInit, need to setup focal size lookup table -- because setting slit
# sizes and energy causes the table to be recalculated, will first disable 
# lookup table calc PV, set the slit sizes and energy, calculate the table 
# directly via python call, then enable calc PV so future slit/energy changes
# update the lookup table

# Disable calc of lookup table
doAfterIocInit("dbpf('$(PREFIX)CRL$(_N):recalc_enable','0')") # Disable table calc

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

# 100idPyCRL
EPICS IOC using PyDevice module to integrate python code for transfocator calculations with lens actuation.  Uses xraylib 
python module

The IOC is based off of the xxx template from APS BCDA synApps collection. For 
more information on synAps, see
   http://www.aps.anl.gov/bcda/synApps

More information on the PyDevice module can be found here:
   https://github.com/klemenv/PyDevice

## Setting up

### Python Environment

Need to create conda environment from which to build and run the IOC (should be kept external to IOC file structure)
Create a source script for activating the environment (e.g. 100idPyFilter_startup_env)


### Startup script changes
 
xxx.pl --> 100idPyCRL.pl

Added definitions: 
PYDEV_STARTUP —> file to be sourced that includes conda activation and updates to LD_LIBRARY_PATH
PYDEV_IOC_CMD —> combined source and IOC_CMD for use in screen call

Changes to command files in iocBoot/ioc100idPyCRL/softioc/commands:
run.pl — had to combine source command and IOC startup command into one line:
```
system("source ${PYDEV_STARTUP} && ${IOC_CMD}");
```
 
start.pl — And in 6-3 looks like:
```
system("$SCREEN -dm -S $IOC_NAME -h 5000 -L -Logfile $LOG_FILE bash -c \"$PYDEV_IOC_CMD\"");
```

### caQtDM script changes

start_caQtDM_100idPyCRL and start_caQtDM_100idPyCRL_double bring up the main screen
for their respective systems (not the default xxx tabbed menu screen). To customize for 
a new IOC, the UI_FILE_MACROS line needs updating, e.g. for transfocators, the original line is:
```
export UI_FILE_MACROS=${2:-"P=100idPyCRL:,CRL=CRL"}
```
P and CRL need to be updated for the new IOC (i.e. match those set in the settings.iocsh and st.cmd).
 
## Running
Should be able run like any other xxx-based (synApps) IOC i.e. by running the xxx.pl file (or as in one the new examples in this IOC 100idPyCRL.pl)


+-----------------------------+
| msotlib.itheramedical.com   |
+-----------------------------+



MATLAB Interfacing Library to the ViewMSOT Proprietary File Format





1. INSTALL
--------------------------------------------------------------------------------------

In order to be able to use loadMSOT, please add the following to your startup.m file:
  
  javaaddpath <DIRECTORY>\MSOTBeans\xbean.jar
  javaaddpath <DIRECTORY>\MSOTBeans\msotbeans.jar

Your startup.m file should reside in your MATLAB Startup directory (find out using
the MATLAB command "userpath". If startup.m does not exist in this directory, please
create it.






2. CONTENTS
--------------------------------------------------------------------------------------

Listing of functions:
- listMSOT	    Lists contents of a study folder
- loadMSOT          Loads MSOT META information
- loadMSOTRecon     Loads MSOT Reconstructions
- loadMSOTMsp       Loads MSOT MSP (multispectrally processed) data
- loadMSOTSignals   Loads optoacoustic Signals as acquired by the transducers






3. USAGE
--------------------------------------------------------------------------------------

For usage please refer to the function documentation using help <function>
# Microsoft Developer Studio Project File - Name="rtheta" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=rtheta - Win32 uBLAS Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "rtheta.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "rtheta.mak" CFG="rtheta - Win32 uBLAS Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "rtheta - Win32 uBLAS Debug" (based on "Win32 (x86) Console Application")
!MESSAGE "rtheta - Win32 uBLAS Release" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""$/BayesFilter/rtheta", UBDAAAAA"
# PROP Scc_LocalPath "."
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "rtheta - Win32 uBLAS Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "rtheta___Win32_uBLAS_Debug"
# PROP BASE Intermediate_Dir "rtheta___Win32_uBLAS_Debug"
# PROP BASE Ignore_Export_Lib 0
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "rtheta_uBLAS_Debug"
# PROP Intermediate_Dir "rtheta_uBLAS_Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /MDd /W3 /Gm /Gi /GX /ZI /Od /I "../MTLmatSup" /I ".." /D "_CONSOLE" /D "_DEBUG" /D "WIN32" /D "_MBCS" /D "MTL_EXCEPTIONS" /FD /GZ /c
# SUBTRACT BASE CPP /YX
# ADD CPP /nologo /MDd /W3 /Gm /Gi /GX /ZI /Od /I "../uBLASmatSup" /I ".." /D "_CONSOLE" /D "_DEBUG" /D "WIN32" /D "_MBCS" /FD /GZ /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0xc09 /d "_DEBUG"
# ADD RSC /l 0xc09 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib comdlg32.lib advapi32.lib shell32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib uuid.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib comdlg32.lib advapi32.lib shell32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib uuid.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept

!ELSEIF  "$(CFG)" == "rtheta - Win32 uBLAS Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "rtheta___Win32_uBLAS_Release"
# PROP BASE Intermediate_Dir "rtheta___Win32_uBLAS_Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "rtheta_uBLAS_Release"
# PROP Intermediate_Dir "rtheta_uBLAS_Release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /MD /W3 /GX /O2 /Op /I "../MTLmatSup" /I ".." /D "_CONSOLE" /D "NDEBUG" /D "WIN32" /D "_MBCS" /D "BAYESFILTER_MTL" /FD /c
# SUBTRACT BASE CPP /YX
# ADD CPP /nologo /MD /W3 /GX /O2 /Op /I "../uBLASmatSup" /I ".." /D "_CONSOLE" /D "NDEBUG" /D "WIN32" /D "_MBCS" /D "BAYESFILTER_MTL" /FD /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0xc09 /d "NDEBUG"
# ADD RSC /l 0xc09 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386

!ENDIF 

# Begin Target

# Name "rtheta - Win32 uBLAS Debug"
# Name "rtheta - Win32 uBLAS Release"
# Begin Source File

SOURCE=.\rtheta.cpp
# End Source File
# Begin Source File

SOURCE=.\TestResult.txt
# End Source File
# Begin Source File

SOURCE=.\tests.cpp
# End Source File
# Begin Source File

SOURCE=..\..\CLAPACK\Debug\clapack_6d.lib

!IF  "$(CFG)" == "rtheta - Win32 uBLAS Debug"

!ELSEIF  "$(CFG)" == "rtheta - Win32 uBLAS Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\CLAPACK\Release\clapack_6.lib

!IF  "$(CFG)" == "rtheta - Win32 uBLAS Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "rtheta - Win32 uBLAS Release"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\VClib\BayesFilter_uBLAS_MD6.lib

!IF  "$(CFG)" == "rtheta - Win32 uBLAS Debug"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "rtheta - Win32 uBLAS Release"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\VClib\BayesFilter_uBLAS_MD6d.lib

!IF  "$(CFG)" == "rtheta - Win32 uBLAS Debug"

!ELSEIF  "$(CFG)" == "rtheta - Win32 uBLAS Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# End Target
# End Project

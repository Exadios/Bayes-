# Microsoft Developer Studio Project File - Name="VClib" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=VClib - Win32 uBLAS Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "VClib.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "VClib.mak" CFG="VClib - Win32 uBLAS Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "VClib - Win32 uBLAS Debug" (based on "Win32 (x86) Static Library")
!MESSAGE "VClib - Win32 uBLAS Release" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""$/BayesFilter/VClib", VADAAAAA"
# PROP Scc_LocalPath "."
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "VClib - Win32 uBLAS Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "VClib___Win32_uBLAS_Debug"
# PROP BASE Intermediate_Dir "VClib___Win32_uBLAS_Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "VClib_uBLAS_Debug"
# PROP Intermediate_Dir "VClib_uBLAS_Debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /MDd /W3 /Gm /GX /Zi /Od /Gf /I "../MTLmatSup" /D "_DEBUG" /D "_LIB" /D "WIN32" /D "_MBCS" /D for=if(1)for /D "MTL_EXCEPTIONS" /Fd"BayesFilter_MTL_MDd.pdb" /FD /GZ /c
# SUBTRACT BASE CPP /Gy /YX
# ADD CPP /nologo /MDd /W3 /Gm /GX /Zi /Od /Gf /I "../uBLASmatSup" /I ".." /D "_LIB" /D for=if(1)for /D "_DEBUG" /D "WIN32" /D "_MBCS" /Fd"BayesFilter_MTL_MDd.pdb" /FD /GZ /Zm150 /c
# ADD BASE RSC /l 0xc09 /d "_DEBUG"
# ADD RSC /l 0xc09 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo /out:"BayesFilter_MTL_MDd.lib"
# ADD LIB32 /nologo /out:"BayesFilter_uBLAS_MD6d.lib"

!ELSEIF  "$(CFG)" == "VClib - Win32 uBLAS Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "VClib___Win32_uBLAS_Release"
# PROP BASE Intermediate_Dir "VClib___Win32_uBLAS_Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "VClib_uBLAS_Release"
# PROP Intermediate_Dir "VClib_uBLAS_Release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /MD /W3 /GX /Ox /Ot /Oa /Og /Oi /Op /Ob2 /Gf /I "../MTLmatSup" /D "NDEBUG" /D "_LIB" /D "WIN32" /D "_MBCS" /D for=if(1)for /FD /c
# SUBTRACT BASE CPP /Gy /YX
# ADD CPP /nologo /MD /W3 /GX /Ox /Ot /Oa /Og /Oi /Op /Ob2 /Gf /I "../uBLASmatSup" /I ".." /D "NDEBUG" /D "_LIB" /D "WIN32" /D "_MBCS" /D for=if(1)for /FD /c
# SUBTRACT CPP /Gy /YX
# ADD BASE RSC /l 0xc09 /d "NDEBUG"
# ADD RSC /l 0xc09 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo /out:"BayesFilter_MTL_MD.lib"
# ADD LIB32 /nologo /out:"BayesFilter_uBLAS_MD6.lib"

!ENDIF 

# Begin Target

# Name "VClib - Win32 uBLAS Debug"
# Name "VClib - Win32 uBLAS Release"
# Begin Group "BayesFilter Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\BayesFilter\bayesFlt.cpp
# End Source File
# Begin Source File

SOURCE=..\BayesFilter\bayesFltAlg.cpp
# End Source File
# Begin Source File

SOURCE=..\BayesFilter\bayesFltDummy.cpp
# PROP BASE Exclude_From_Build 1
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=..\BayesFilter\covFlt.cpp
# End Source File
# Begin Source File

SOURCE=..\BayesFilter\infFlt.cpp
# End Source File
# Begin Source File

SOURCE=..\BayesFilter\infRtFlt.cpp
# End Source File
# Begin Source File

SOURCE=..\BayesFilter\itrFlt.cpp
# End Source File
# Begin Source File

SOURCE=..\BayesFilter\matSup.cpp
# End Source File
# Begin Source File

SOURCE=..\BayesFilter\SIRFlt.cpp
# End Source File
# Begin Source File

SOURCE=..\BayesFilter\UDFlt.cpp
# End Source File
# Begin Source File

SOURCE=..\BayesFilter\UdU.cpp
# End Source File
# Begin Source File

SOURCE=..\BayesFilter\unsFlt.cpp
# End Source File
# End Group
# Begin Group "BayesFilter Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\BayesFilter\allFlt.h
# End Source File
# Begin Source File

SOURCE=..\BayesFilter\bayesFlt.h
# End Source File
# Begin Source File

SOURCE=..\BayesFilter\compatibility.h
# End Source File
# Begin Source File

SOURCE=..\BayesFilter\covFlt.h
# End Source File
# Begin Source File

SOURCE=..\BayesFilter\infFlt.h
# End Source File
# Begin Source File

SOURCE=..\BayesFilter\infRtFlt.h
# End Source File
# Begin Source File

SOURCE=..\BayesFilter\itrFlt.h
# End Source File
# Begin Source File

SOURCE=..\BayesFilter\matSup.h
# End Source File
# Begin Source File

SOURCE=..\BayesFilter\models.h
# End Source File
# Begin Source File

SOURCE=..\BayesFilter\schemeFlt.h
# End Source File
# Begin Source File

SOURCE=..\BayesFilter\SIRFlt.h
# End Source File
# Begin Source File

SOURCE=..\BayesFilter\UDFlt.h
# End Source File
# Begin Source File

SOURCE=..\BayesFilter\unsFlt.h
# End Source File
# End Group
# Begin Group "Matrix Support"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\fastMatSup\matSupSub.h
# End Source File
# Begin Source File

SOURCE=..\uBLASmatSup\uBLASmatrix.h
# End Source File
# End Group
# Begin Group "Filters"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\BayesFilter\filters\average1.h
# End Source File
# Begin Source File

SOURCE=..\BayesFilter\filters\indirect.h
# End Source File
# End Group
# Begin Source File

SOURCE="..\BayesFilter TODO.txt"
# End Source File
# End Target
# End Project

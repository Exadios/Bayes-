# Microsoft Developer Studio Project File - Name="MuPAD" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) External Target" 0x0106

CFG=MuPAD - Win32 MTL Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "MuPAD.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "MuPAD.mak" CFG="MuPAD - Win32 MTL Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "MuPAD - Win32 MTL Debug" (based on "Win32 (x86) External Target")
!MESSAGE "MuPAD - Win32 MTL Release" (based on "Win32 (x86) External Target")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""$/BayesFilter/MuPAD", CCGAAAAA"
# PROP Scc_LocalPath "."

!IF  "$(CFG)" == "MuPAD - Win32 MTL Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "MuPAD___Win32_MTL_Debug"
# PROP BASE Intermediate_Dir "MuPAD___Win32_MTL_Debug"
# PROP BASE Cmd_Line "build"
# PROP BASE Rebuild_Opt "/a"
# PROP BASE Target_File "bfilter.mdm"
# PROP BASE Bsc_Name ""
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "MuPAD___Win32_MTL_Debug"
# PROP Intermediate_Dir "MuPAD___Win32_MTL_Debug"
# PROP Cmd_Line "nmake debug"
# PROP Rebuild_Opt "/a"
# PROP Target_File "bfilterD.mdm"
# PROP Bsc_Name ""
# PROP Target_Dir ""

!ELSEIF  "$(CFG)" == "MuPAD - Win32 MTL Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "MuPAD___Win32_MTL_Release"
# PROP BASE Intermediate_Dir "MuPAD___Win32_MTL_Release"
# PROP BASE Cmd_Line "NMAKE /f MuPAD.mak"
# PROP BASE Rebuild_Opt "/a"
# PROP BASE Target_File "MuPAD.exe"
# PROP BASE Bsc_Name "MuPAD.bsc"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "MuPAD___Win32_MTL_Release"
# PROP Intermediate_Dir "MuPAD___Win32_MTL_Release"
# PROP Cmd_Line "nmake release"
# PROP Rebuild_Opt "/a"
# PROP Target_File "bfilterR.mdm"
# PROP Bsc_Name ""
# PROP Target_Dir ""

!ENDIF 

# Begin Target

# Name "MuPAD - Win32 MTL Debug"
# Name "MuPAD - Win32 MTL Release"

!IF  "$(CFG)" == "MuPAD - Win32 MTL Debug"

!ELSEIF  "$(CFG)" == "MuPAD - Win32 MTL Release"

!ENDIF 

# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=.\bfilter.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=.\MuPadConvert.h
# End Source File
# End Group
# Begin Source File

SOURCE=.\Makefile
# End Source File
# End Target
# End Project

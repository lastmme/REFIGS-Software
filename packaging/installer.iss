#ifndef AppVersion
  #define AppVersion "1.0.0"
#endif
[Setup]
AppId={{DB3BCAD2-A1BB-40E7-A2ED-7F57F76EE932}
AppName=RE-FIGS
AppVersion={#AppVersion}
AppPublisher=RE-FIGS
DefaultDirName={localappdata}\Programs\RE-FIGS
DefaultGroupName=RE-FIGS
PrivilegesRequired=lowest
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
OutputDir=..\out
OutputBaseFilename=RE-FIGS-{#AppVersion}-windows-x64-setup
SetupIconFile=..\REFIGS.ico
UninstallDisplayIcon={app}\RE-FIGS.exe
Compression=lzma2
SolidCompression=yes
WizardStyle=modern
[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"
[Tasks]
Name: "desktopicon"; Description: "Create a desktop shortcut"; Flags: unchecked
[Files]
Source: "..\out\dist\RE-FIGS\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs
[Icons]
Name: "{group}\RE-FIGS"; Filename: "{app}\RE-FIGS.exe"
Name: "{autodesktop}\RE-FIGS"; Filename: "{app}\RE-FIGS.exe"; Tasks: desktopicon
[Run]
Filename: "{app}\RE-FIGS.exe"; Description: "Launch RE-FIGS"; Flags: nowait postinstall skipifsilent

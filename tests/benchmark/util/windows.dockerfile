# escape=`

FROM mcr.microsoft.com/dotnet/framework/sdk:4.8-windowsservercore-2004

# Restore the default Windows shell for correct batch processing.
SHELL ["cmd", "/S", "/C"]

# Disable UAC - this is needed to be able to create symlinks.
RUN reg add HKLM\Software\Microsoft\Windows\CurrentVersion\policies\system /v EnableLUA /t Reg_DWORD /d 0 /f

# Download the Build Tools
ADD https://github.com/git-for-windows/git/releases/download/v2.28.0.windows.1/Git-2.28.0-64-bit.exe C:\TEMP\git-installer.exe
ADD https://aka.ms/vs/16/release/vs_buildtools.exe C:\TEMP\vs_buildtools.exe
ADD https://www.python.org/ftp/python/3.8.5/python-3.8.5-amd64.exe C:\TEMP\python-installer.exe

# Install git
RUN C:\TEMP\git-installer.exe /VERYSILENT /NORESTART /NOCANCEL /SP- /CLOSEAPPLICATIONS /RESTARTAPPLICATIONS /COMPONENTS=""

# Add git tools to the PATH
RUN setx /M PATH "%PATH%;C:\Program Files\Git\usr\bin"

# Associate .sh with bash.exe
RUN ftype Bash.Script="C:\Program Files\Git\usr\bin\bash.exe" %L %* && assoc .sh=Bash.Script

# Install Python
RUN C:\TEMP\python-installer.exe /quiet InstallAllUsers=1 TargetDir=C:\Python PrependPath=1 Include_test=0 Include_doc=0

# Install MSVC Build Tools
RUN C:\TEMP\vs_buildtools.exe --quiet --wait --norestart --nocache `
    --installPath C:\BuildTools `
    --add Microsoft.VisualStudio.Workload.VCTools `
    --add Microsoft.VisualStudio.Component.VC.CMake.Project `
    --add Microsoft.VisualStudio.Component.VC.Tools.x86.x64 `
    --add Microsoft.VisualStudio.Component.Windows10SDK.18362 `
    --remove Microsoft.VisualStudio.Component.Windows10SDK.10240 `
    --remove Microsoft.VisualStudio.Component.Windows10SDK.10586 `
    --remove Microsoft.VisualStudio.Component.Windows10SDK.14393 `
    --remove Microsoft.VisualStudio.Component.Windows81SDK `
    || IF "%ERRORLEVEL%"=="3010" EXIT 0

# Define the entry point for the docker container.
ENV BUILD_ARCH=x64
ENTRYPOINT ["C:\\BuildTools\\Common7\\Tools\\VsDevCmd.bat", "-arch=%BUILD_ARCH%", "-host_arch=amd64", "&&", "cmd.exe"]
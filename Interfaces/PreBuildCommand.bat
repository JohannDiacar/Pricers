echo Parsing of Interface Functions
date /t

SET interfaceGen=%1

del report_file.txt

for %%v in (InterfaceFunctions\*.h) do (%InterfaceGen% %%v >> report_file.txt)
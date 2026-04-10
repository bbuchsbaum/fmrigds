@echo off
R --no-save --no-restore -s -e "fmrigds:::fmrigds_cli_exec()" --args %*


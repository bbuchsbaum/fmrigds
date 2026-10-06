@echo off
R --no-save --no-restore -s -e "quit(save = 'no', status = fmrigds:::fmrigds_cli_exec(), runLast = FALSE)" --args %*


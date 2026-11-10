# scripts

The command-line tools comRAD installs, unchanged, under their usual names
(`trim2bRAD.py`, `pickQ.py`, `mind`, ...). The list is `tools/SCRIPTS`;
`python3 -m pip install` puts each one on `PATH` (`~/.local/bin` with `--user`).

To add or change one: edit it here, add its name to `tools/SCRIPTS` if it is new,
commit, and reinstall (`python3 -m pip install --user .` from a clone, or the
`git+https://` install after pushing).

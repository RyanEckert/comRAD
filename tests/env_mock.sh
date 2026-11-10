# source this to use an existing mock HPC (after setup_mock.sh) without wiping it
T=${T:-${TMPDIR:-/tmp}/comrad-mock}
export PATH=$T/bin:${COMRAD_VENV:+$COMRAD_VENV/bin:}$PATH koko=$T/koko EMAIL=test@example.org USER=me TMPDIR=$T/tmp RTJOBS_POLL=1
P=$koko/2bRAD/proj
python3 -c 'import moments, demes' 2>/dev/null || { mkdir -p $T/py; touch $T/py/moments.py $T/py/demes.py; export PYTHONPATH=$T/py; }

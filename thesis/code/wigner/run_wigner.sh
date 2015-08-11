#!/usr/bin/env bash

nucleus="Be"
path="./${nucleus}"
path="${path}/"
path="$path$nucleus"

echo "Calculating Wigner"

python wigner_0s.py "$nucleus" 0 3 0.02 &
python wigner.py "$nucleus" 0 3 0.02 0 1 0 &
python wigner.py "$nucleus" 0 3 0.02 0 1 1 

wait

echo "done Be"

nucleus="Al"
path="./${nucleus}"
path="${path}/"
path="$path$nucleus"

echo "Calculating Wigner"

python wigner_0s.py "$nucleus" 0 3 0.02 &
python wigner.py "$nucleus" 0 3 0.02 0 1 0 &
python wigner.py "$nucleus" 0 3 0.02 0 1 1 &
python wigner.py "$nucleus" 0 3 0.02 0 2 0 &
python wigner.py "$nucleus" 0 3 0.02 0 2 1 &
python wigner.py "$nucleus" 0 3 0.02 0 2 2 &

wait

echo "done Al"

nucleus="Fe"
path="./${nucleus}"
path="${path}/"
path="$path$nucleus"

echo "Calculating Wigner"

python wigner_0s.py "$nucleus" 0 3 0.02 &
python wigner.py "$nucleus" 0 3 0.02 0 1 0 &
python wigner.py "$nucleus" 0 3 0.02 0 1 1 &
python wigner.py "$nucleus" 0 3 0.02 0 2 0 &
python wigner.py "$nucleus" 0 3 0.02 0 2 1 &
python wigner.py "$nucleus" 0 3 0.02 0 2 2 &
python wigner.py "$nucleus" 0 3 0.02 1 0 0 &

wait

python wigner.py "$nucleus" 0 3 0.02 0 3 0 &
python wigner.py "$nucleus" 0 3 0.02 0 3 1 &
python wigner.py "$nucleus" 0 3 0.02 0 3 2 &
python wigner.py "$nucleus" 0 3 0.02 0 3 3 &
python wigner.py "$nucleus" 0 3 0.02 1 1 0 &
python wigner.py "$nucleus" 0 3 0.02 1 1 1 &

wait 

echo "done Fe"


nucleus="Ca48"
path="./${nucleus}"
path="${path}/"
path="$path$nucleus"

echo "Calculating Wigner"

python wigner_0s.py "$nucleus" 0 3 0.02 &
python wigner.py "$nucleus" 0 3 0.02 0 1 0 &
python wigner.py "$nucleus" 0 3 0.02 0 1 1 &
python wigner.py "$nucleus" 0 3 0.02 0 2 0 &
python wigner.py "$nucleus" 0 3 0.02 0 2 1 &
python wigner.py "$nucleus" 0 3 0.02 0 2 2 &
python wigner.py "$nucleus" 0 3 0.02 1 0 0 &

wait

python wigner.py "$nucleus" 0 3 0.02 0 3 0 &
python wigner.py "$nucleus" 0 3 0.02 0 3 1 &
python wigner.py "$nucleus" 0 3 0.02 0 3 2 &
python wigner.py "$nucleus" 0 3 0.02 0 3 3 &

wait 

echo "done Fe"







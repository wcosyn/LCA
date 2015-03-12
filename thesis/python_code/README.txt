#uitleg bij bestandjes:
1. onebmomentum.py berekent one-body momentum distributie voor verschillende kernen #werkt prima
2. moshinsky.py berekent clebsch-gordan coefficienten en moshinsky brackets voor de two-body momentum distributie #probleem
3. twobody.py berekent two-body relative momentum distr volgens formule (89) uit men voorlopige thesis (Nuclearmomentumdistributions.pdf) #probleem
4. He_twobody.py calculates two-body cm momentum distribution following thesis M. Vanhalst, p74-75 #werkt prima
5. full-shell.py berekent two-body total momentum distributie voor kernen met volledig gevulde schillen (N) (wood-saxon opvulling met Spin-orbit = HO opvulling) (matrices.py, read.py zijn hulp-scripts)


(2.) en (3.) werken niet dus heb ik eerst voor eenvoudigere gevallen geprobeerd (4. en 5.). Bij (5.) geeft hetzelfde als (4.) voor Helium maar voor zwaardere kernen is er nog een factor 3.299 in de normalisatie
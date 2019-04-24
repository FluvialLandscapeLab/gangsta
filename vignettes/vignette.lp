/******************** OPTIMIZATION ********************/
/***** Optimization function (Exprsn. 1) *****/
MAX: Total_Biomass.finalMolecules;
/***** The Total_Biomass.finalMolecules equals the sum of the Organism.finalMolecules for each organism type (Exprsn. 2) *****/
Total_Biomass.finalMolecules = Aut.finalMolecules + Met.finalMolecules;
/******************** INITIALIZATION ********************/
/***** Set FiniteCompound.initialMolecules & InfiniteCompound.initialMolecules *****/
CO2.initialMolecules = 0.2;
CH4.initialMolecules = 0.1;
NH4.initialMolecules = 0.3;
NO3.initialMolecules = 0.1;
O2.initialMolecules = 0.9;
DOM.initialMolecules = 0.3;
Ox.initialMolecules = 0;
Aut.initialMolecules = 1;
Met.initialMolecules = 1;
/***** For each InfiniteCompound, remove constraints on InfiniteCompound.finalMolecules *****/
-Inf < Ox.finalMolecules < +Inf;
/***** For each Pool associated with an InfiniteCompound, remove constraints on Pool.finalAtoms *****/
-Inf < Ox_O.finalAtoms < +Inf;
-Inf < Ox_O.initialAtoms < +Inf;
/***** For each Process, the Process.netEnergy is unconstrained *****/
-Inf < AutNitrif.netEnergy < +Inf;
-Inf < MetMethaneOxid.netEnergy < +Inf;
-Inf < AutAssimCO2.netEnergy < +Inf;
-Inf < MetAssimCH4.netEnergy < +Inf;
-Inf < AutAssimNO3.netEnergy < +Inf;
-Inf < MetAssimNO3.netEnergy < +Inf;
-Inf < AutAssimNH4.netEnergy < +Inf;
-Inf < MetAssimNH4.netEnergy < +Inf;
-Inf < AutDecay.netEnergy < +Inf;
-Inf < MetDecay.netEnergy < +Inf;
/***** For each Organism type, the Organism.respirationEnergy is unconstrained *****/
-Inf < Aut.respirationEnergy < +Inf;
-Inf < Met.respirationEnergy < +Inf;
/******************** STOICHIOMETRY ********************/
/***** For each elemental Pool, the Pool.initialAtoms must conform to Compound stoichiometry (Exprsn. 3) *****/
CO2_C.initialAtoms = 1 CO2.initialMolecules;
CO2_O.initialAtoms = 2 CO2.initialMolecules;
CH4_C.initialAtoms = 1 CH4.initialMolecules;
NH4_N.initialAtoms = 1 NH4.initialMolecules;
NO3_N.initialAtoms = 1 NO3.initialMolecules;
NO3_O.initialAtoms = 3 NO3.initialMolecules;
O2_O.initialAtoms = 2 O2.initialMolecules;
DOM_C.initialAtoms = 1 DOM.initialMolecules;
DOM_N.initialAtoms = 0.150943396226415 DOM.initialMolecules;
DOM_O.initialAtoms = 1.0377358490566 DOM.initialMolecules;
Ox_O.initialAtoms = 1 Ox.initialMolecules;
Aut_C.initialAtoms = 1 Aut.initialMolecules;
Aut_N.initialAtoms = 0.150943396226415 Aut.initialMolecules;
Aut_O.initialAtoms = 1.0377358490566 Aut.initialMolecules;
Met_C.initialAtoms = 1 Met.initialMolecules;
Met_N.initialAtoms = 0.150943396226415 Met.initialMolecules;
Met_O.initialAtoms = 1.0377358490566 Met.initialMolecules;
/***** For each elemental Pool, the Pool.finalAtoms must conform to Compound stoichiometry (Exprsn. 4) *****/
CO2_C.finalAtoms = 1 CO2.finalMolecules;
CO2_O.finalAtoms = 2 CO2.finalMolecules;
CH4_C.finalAtoms = 1 CH4.finalMolecules;
NH4_N.finalAtoms = 1 NH4.finalMolecules;
NO3_N.finalAtoms = 1 NO3.finalMolecules;
NO3_O.finalAtoms = 3 NO3.finalMolecules;
O2_O.finalAtoms = 2 O2.finalMolecules;
DOM_C.finalAtoms = 1 DOM.finalMolecules;
DOM_N.finalAtoms = 0.150943396226415 DOM.finalMolecules;
DOM_O.finalAtoms = 1.0377358490566 DOM.finalMolecules;
Ox_O.finalAtoms = 1 Ox.finalMolecules;
Aut_C.finalAtoms = 1 Aut.finalMolecules;
Aut_N.finalAtoms = 0.150943396226415 Aut.finalMolecules;
Aut_O.finalAtoms = 1.0377358490566 Aut.finalMolecules;
Met_C.finalAtoms = 1 Met.finalMolecules;
Met_N.finalAtoms = 0.150943396226415 Met.finalMolecules;
Met_O.finalAtoms = 1.0377358490566 Met.finalMolecules;
/******************** MOLAR TRANSFER ********************/
/***** The number of atoms transferred during each Process must be in accordance with the stoichiometry and chemical affinity of the Process (Exprsn. 5) *****/
AutNitrif.netEnergy = 0.0003485 AutNitrif_NH4_N_NO3_N.atoms;
AutNitrif.netEnergy = 0.000116166666666667 AutNitrif_O2_O_NO3_O.atoms;
AutNitrif.netEnergy = 0.0003485 AutNitrif_O2_O_Ox_O.atoms;
MetMethaneOxid.netEnergy = 0.000818 MetMethaneOxid_CH4_C_CO2_C.atoms;
MetMethaneOxid.netEnergy = 0.000409 MetMethaneOxid_O2_O_CO2_O.atoms;
MetMethaneOxid.netEnergy = 0.000409 MetMethaneOxid_O2_O_Ox_O.atoms;
AutAssimCO2.netEnergy = -0.0035 AutAssimCO2_CO2_C_Aut_C.atoms;
AutAssimCO2.netEnergy = -0.00175 AutAssimCO2_CO2_O_Aut_O.atoms + -0.00175 AutAssimCO2_CO2_O_Ox_O.atoms;
MetAssimCH4.netEnergy = -0.00109 MetAssimCH4_CH4_C_Met_C.atoms;
AutAssimNO3.netEnergy = -0.000155 AutAssimNO3_NO3_N_Aut_N.atoms;
AutAssimNO3.netEnergy = -5.16666666666667e-05 AutAssimNO3_NO3_O_Aut_O.atoms + -5.16666666666667e-05 AutAssimNO3_NO3_O_Ox_O.atoms;
MetAssimNO3.netEnergy = -0.000155 MetAssimNO3_NO3_N_Met_N.atoms;
MetAssimNO3.netEnergy = -5.16666666666667e-05 MetAssimNO3_NO3_O_Met_O.atoms + -5.16666666666667e-05 MetAssimNO3_NO3_O_Ox_O.atoms;
AutAssimNH4.netEnergy = -3.18e-05 AutAssimNH4_NH4_N_Aut_N.atoms;
MetAssimNH4.netEnergy = -3.18e-05 MetAssimNH4_NH4_N_Met_N.atoms;
/******************** TURNOVER ********************/
/***** For each organism type, decay must be at least equal to turnover rate *****/
AutDecay_Aut_C_DOM_C.atoms >= 0.05 Aut_C.initialAtoms;
AutDecay_Aut_N_DOM_N.atoms >= 0.05 Aut_N.initialAtoms;
AutDecay_Aut_O_DOM_O.atoms >= 0.05 Aut_O.initialAtoms;
MetDecay_Met_C_DOM_C.atoms >= 0.05 Met_C.initialAtoms;
MetDecay_Met_N_DOM_N.atoms >= 0.05 Met_N.initialAtoms;
MetDecay_Met_O_DOM_O.atoms >= 0.05 Met_O.initialAtoms;
/******************** MAXIMUM GROWTH RATE ********************/
/***** For each organism type, growth may not exceed a maximum rate *****/
Aut_C.finalAtoms <= 1.5 Aut_C.initialAtoms;
Aut_N.finalAtoms <= 1.5 Aut_N.initialAtoms;
Aut_O.finalAtoms <= 1.5 Aut_O.initialAtoms;
Met_C.finalAtoms <= 1.5 Met_C.initialAtoms;
Met_N.finalAtoms <= 1.5 Met_N.initialAtoms;
Met_O.finalAtoms <= 1.5 Met_O.initialAtoms;
/******************** ELEMENTAL MOLAR BALANCE ********************/
/***** Elemental molar balance must be conserved (Exprsn. 6) *****/
CO2_C.finalAtoms = CO2_C.initialAtoms + MetMethaneOxid_CH4_C_CO2_C.atoms - AutAssimCO2_CO2_C_Aut_C.atoms;
CO2_O.finalAtoms = CO2_O.initialAtoms + MetMethaneOxid_O2_O_CO2_O.atoms - AutAssimCO2_CO2_O_Aut_O.atoms - AutAssimCO2_CO2_O_Ox_O.atoms;
CH4_C.finalAtoms = CH4_C.initialAtoms - MetMethaneOxid_CH4_C_CO2_C.atoms - MetAssimCH4_CH4_C_Met_C.atoms;
NH4_N.finalAtoms = NH4_N.initialAtoms - AutNitrif_NH4_N_NO3_N.atoms - AutAssimNH4_NH4_N_Aut_N.atoms - MetAssimNH4_NH4_N_Met_N.atoms;
NO3_N.finalAtoms = NO3_N.initialAtoms + AutNitrif_NH4_N_NO3_N.atoms - AutAssimNO3_NO3_N_Aut_N.atoms - MetAssimNO3_NO3_N_Met_N.atoms;
NO3_O.finalAtoms = NO3_O.initialAtoms + AutNitrif_O2_O_NO3_O.atoms - AutAssimNO3_NO3_O_Aut_O.atoms - AutAssimNO3_NO3_O_Ox_O.atoms - MetAssimNO3_NO3_O_Met_O.atoms - MetAssimNO3_NO3_O_Ox_O.atoms;
O2_O.finalAtoms = O2_O.initialAtoms - AutNitrif_O2_O_NO3_O.atoms - AutNitrif_O2_O_Ox_O.atoms - MetMethaneOxid_O2_O_CO2_O.atoms - MetMethaneOxid_O2_O_Ox_O.atoms;
DOM_C.finalAtoms = DOM_C.initialAtoms + AutDecay_Aut_C_DOM_C.atoms + MetDecay_Met_C_DOM_C.atoms;
DOM_N.finalAtoms = DOM_N.initialAtoms + AutDecay_Aut_N_DOM_N.atoms + MetDecay_Met_N_DOM_N.atoms;
DOM_O.finalAtoms = DOM_O.initialAtoms + AutDecay_Aut_O_DOM_O.atoms + MetDecay_Met_O_DOM_O.atoms;
Ox_O.finalAtoms = Ox_O.initialAtoms + AutNitrif_O2_O_Ox_O.atoms + MetMethaneOxid_O2_O_Ox_O.atoms + AutAssimCO2_CO2_O_Ox_O.atoms + AutAssimNO3_NO3_O_Ox_O.atoms + MetAssimNO3_NO3_O_Ox_O.atoms;
Aut_C.finalAtoms = Aut_C.initialAtoms + AutAssimCO2_CO2_C_Aut_C.atoms - AutDecay_Aut_C_DOM_C.atoms;
Aut_N.finalAtoms = Aut_N.initialAtoms + AutAssimNO3_NO3_N_Aut_N.atoms + AutAssimNH4_NH4_N_Aut_N.atoms - AutDecay_Aut_N_DOM_N.atoms;
Aut_O.finalAtoms = Aut_O.initialAtoms + AutAssimCO2_CO2_O_Aut_O.atoms + AutAssimNO3_NO3_O_Aut_O.atoms - AutDecay_Aut_O_DOM_O.atoms;
Met_C.finalAtoms = Met_C.initialAtoms + MetAssimCH4_CH4_C_Met_C.atoms - MetDecay_Met_C_DOM_C.atoms;
Met_N.finalAtoms = Met_N.initialAtoms + MetAssimNO3_NO3_N_Met_N.atoms + MetAssimNH4_NH4_N_Met_N.atoms - MetDecay_Met_N_DOM_N.atoms;
Met_O.finalAtoms = Met_O.initialAtoms + MetAssimNO3_NO3_O_Met_O.atoms - MetDecay_Met_O_DOM_O.atoms;
/***** Generally, Transfers out of each Pool must be less than the Pool.initialAtoms (Exprsn. 7) *****/
CO2_C.initialAtoms >= AutAssimCO2_CO2_C_Aut_C.atoms;
CO2_O.initialAtoms >= AutAssimCO2_CO2_O_Aut_O.atoms + AutAssimCO2_CO2_O_Ox_O.atoms;
CH4_C.initialAtoms >= MetMethaneOxid_CH4_C_CO2_C.atoms + MetAssimCH4_CH4_C_Met_C.atoms;
NH4_N.initialAtoms >= AutNitrif_NH4_N_NO3_N.atoms + AutAssimNH4_NH4_N_Aut_N.atoms + MetAssimNH4_NH4_N_Met_N.atoms;
NO3_N.initialAtoms >= AutAssimNO3_NO3_N_Aut_N.atoms + MetAssimNO3_NO3_N_Met_N.atoms;
NO3_O.initialAtoms >= AutAssimNO3_NO3_O_Aut_O.atoms + AutAssimNO3_NO3_O_Ox_O.atoms + MetAssimNO3_NO3_O_Met_O.atoms + MetAssimNO3_NO3_O_Ox_O.atoms;
O2_O.initialAtoms >= AutNitrif_O2_O_NO3_O.atoms + AutNitrif_O2_O_Ox_O.atoms + MetMethaneOxid_O2_O_CO2_O.atoms + MetMethaneOxid_O2_O_Ox_O.atoms;
Aut_C.initialAtoms >= AutDecay_Aut_C_DOM_C.atoms;
Aut_N.initialAtoms >= AutDecay_Aut_N_DOM_N.atoms;
Aut_O.initialAtoms >= AutDecay_Aut_O_DOM_O.atoms;
Met_C.initialAtoms >= MetDecay_Met_C_DOM_C.atoms;
Met_N.initialAtoms >= MetDecay_Met_N_DOM_N.atoms;
Met_O.initialAtoms >= MetDecay_Met_O_DOM_O.atoms;
/******************** ENERGETIC CONSTRAINTS ********************/
/***** Maintenance respiration energy is proportional to the biomass of each organism type (Exprsn. 8) *****/
Aut.respirationEnergy = -6.792e-05 Aut.finalMolecules;
Met.respirationEnergy = -6.792e-05 Met.finalMolecules;
/***** For each organism type, the total energy associated with assimilatory and dissimilatory Processes must equal the sum of the Process.netEnergy for all Processes (Exprsn. 9) *****/
Aut.totalProcessEnergy = AutNitrif.netEnergy + AutAssimCO2.netEnergy + AutAssimNO3.netEnergy + AutAssimNH4.netEnergy;
Met.totalProcessEnergy = MetMethaneOxid.netEnergy + MetAssimCH4.netEnergy + MetAssimNO3.netEnergy + MetAssimNH4.netEnergy;
/***** Organisms must meet their respiratory energy demands (Exprsn. 10) *****/
0 = Aut.respirationEnergy + Aut.totalProcessEnergy;
0 = Met.respirationEnergy + Met.totalProcessEnergy;

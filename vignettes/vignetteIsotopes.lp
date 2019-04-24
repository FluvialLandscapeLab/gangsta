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
-Inf < Ox_O16.finalAtoms < +Inf;
-Inf < Ox_O17.finalAtoms < +Inf;
-Inf < Ox_O18.finalAtoms < +Inf;
-Inf < Ox_O16.initialAtoms < +Inf;
-Inf < Ox_O17.initialAtoms < +Inf;
-Inf < Ox_O18.initialAtoms < +Inf;
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
NH4_N.initialAtoms = 1 NH4.initialMolecules;
NO3_N.initialAtoms = 1 NO3.initialMolecules;
DOM_N.initialAtoms = 0.150943396226415 DOM.initialMolecules;
Aut_N.initialAtoms = 0.150943396226415 Aut.initialMolecules;
Met_N.initialAtoms = 0.150943396226415 Met.initialMolecules;
CO2_C12.initialAtoms + CO2_C13.initialAtoms = 1 CO2.initialMolecules;
CH4_C12.initialAtoms + CH4_C13.initialAtoms = 1 CH4.initialMolecules;
DOM_C12.initialAtoms + DOM_C13.initialAtoms = 1 DOM.initialMolecules;
Aut_C12.initialAtoms + Aut_C13.initialAtoms = 1 Aut.initialMolecules;
Met_C12.initialAtoms + Met_C13.initialAtoms = 1 Met.initialMolecules;
CO2_O16.initialAtoms + CO2_O17.initialAtoms + CO2_O18.initialAtoms = 2 CO2.initialMolecules;
NO3_O16.initialAtoms + NO3_O17.initialAtoms + NO3_O18.initialAtoms = 3 NO3.initialMolecules;
O2_O16.initialAtoms + O2_O17.initialAtoms + O2_O18.initialAtoms = 2 O2.initialMolecules;
DOM_O16.initialAtoms + DOM_O17.initialAtoms + DOM_O18.initialAtoms = 1.0377358490566 DOM.initialMolecules;
Ox_O16.initialAtoms + Ox_O17.initialAtoms + Ox_O18.initialAtoms = 1 Ox.initialMolecules;
Aut_O16.initialAtoms + Aut_O17.initialAtoms + Aut_O18.initialAtoms = 1.0377358490566 Aut.initialMolecules;
Met_O16.initialAtoms + Met_O17.initialAtoms + Met_O18.initialAtoms = 1.0377358490566 Met.initialMolecules;
/***** For each elemental Pool, the Pool.finalAtoms must conform to Compound stoichiometry (Exprsn. 4) *****/
NH4_N.finalAtoms = 1 NH4.finalMolecules;
NO3_N.finalAtoms = 1 NO3.finalMolecules;
DOM_N.finalAtoms = 0.150943396226415 DOM.finalMolecules;
Aut_N.finalAtoms = 0.150943396226415 Aut.finalMolecules;
Met_N.finalAtoms = 0.150943396226415 Met.finalMolecules;
CO2_C12.finalAtoms + CO2_C13.finalAtoms = 1 CO2.finalMolecules;
CH4_C12.finalAtoms + CH4_C13.finalAtoms = 1 CH4.finalMolecules;
DOM_C12.finalAtoms + DOM_C13.finalAtoms = 1 DOM.finalMolecules;
Aut_C12.finalAtoms + Aut_C13.finalAtoms = 1 Aut.finalMolecules;
Met_C12.finalAtoms + Met_C13.finalAtoms = 1 Met.finalMolecules;
CO2_O16.finalAtoms + CO2_O17.finalAtoms + CO2_O18.finalAtoms = 2 CO2.finalMolecules;
NO3_O16.finalAtoms + NO3_O17.finalAtoms + NO3_O18.finalAtoms = 3 NO3.finalMolecules;
O2_O16.finalAtoms + O2_O17.finalAtoms + O2_O18.finalAtoms = 2 O2.finalMolecules;
DOM_O16.finalAtoms + DOM_O17.finalAtoms + DOM_O18.finalAtoms = 1.0377358490566 DOM.finalMolecules;
Ox_O16.finalAtoms + Ox_O17.finalAtoms + Ox_O18.finalAtoms = 1 Ox.finalMolecules;
Aut_O16.finalAtoms + Aut_O17.finalAtoms + Aut_O18.finalAtoms = 1.0377358490566 Aut.finalMolecules;
Met_O16.finalAtoms + Met_O17.finalAtoms + Met_O18.finalAtoms = 1.0377358490566 Met.finalMolecules;
/******************** MOLAR TRANSFER ********************/
/***** The number of atoms transferred during each Process must be in accordance with the stoichiometry and chemical affinity of the Process (Exprsn. 5) *****/
AutNitrif.netEnergy = 0.0003485 AutNitrif_NH4_N_NO3_N.atoms;
AutNitrif.netEnergy = 0.000116166666666667 AutNitrif_O2_O16_NO3_O16.atoms;
AutNitrif.netEnergy = 0.000116166666666667 AutNitrif_O2_O17_NO3_O17.atoms;
MetMethaneOxid.netEnergy = 0.000818 MetMethaneOxid_CH4_C12_CO2_C12.atoms;
MetMethaneOxid.netEnergy = 0.000818 MetMethaneOxid_CH4_C13_CO2_C13.atoms;
MetMethaneOxid.netEnergy = 0.000409 MetMethaneOxid_O2_O16_CO2_O16.atoms;
AutAssimCO2.netEnergy = -0.0035 AutAssimCO2_CO2_C12_Aut_C12.atoms;
AutAssimCO2.netEnergy = -0.0035 AutAssimCO2_CO2_C13_Aut_C13.atoms + -0.00175 AutAssimCO2_CO2_O16_Aut_O16.atoms;
MetAssimCH4.netEnergy = -0.00109 MetAssimCH4_CH4_C12_Met_C12.atoms;
AutAssimNO3.netEnergy = -0.000155 AutAssimNO3_NO3_N_Aut_N.atoms;
AutAssimNO3.netEnergy = -5.16666666666667e-05 AutAssimNO3_NO3_O16_Aut_O16.atoms + -5.16666666666667e-05 AutAssimNO3_NO3_O17_Aut_O17.atoms;
MetAssimNO3.netEnergy = -0.000155 MetAssimNO3_NO3_N_Met_N.atoms;
MetAssimNO3.netEnergy = -5.16666666666667e-05 MetAssimNO3_NO3_O16_Met_O16.atoms + -5.16666666666667e-05 MetAssimNO3_NO3_O17_Met_O17.atoms;
AutAssimNH4.netEnergy = -3.18e-05 AutAssimNH4_NH4_N_Aut_N.atoms;
MetAssimNH4.netEnergy = -3.18e-05 MetAssimNH4_NH4_N_Met_N.atoms;
/******************** ELEMENTAL MOLAR BALANCE ********************/
/***** Elemental molar balance must be conserved (Exprsn. 6) *****/
NH4_N.finalAtoms = NH4_N.initialAtoms - AutNitrif_NH4_N_NO3_N.atoms - AutAssimNH4_NH4_N_Aut_N.atoms - MetAssimNH4_NH4_N_Met_N.atoms;
NO3_N.finalAtoms = NO3_N.initialAtoms + AutNitrif_NH4_N_NO3_N.atoms - AutAssimNO3_NO3_N_Aut_N.atoms - MetAssimNO3_NO3_N_Met_N.atoms;
DOM_N.finalAtoms = DOM_N.initialAtoms + AutDecay_Aut_N_DOM_N.atoms + MetDecay_Met_N_DOM_N.atoms;
Aut_N.finalAtoms = Aut_N.initialAtoms + AutAssimNO3_NO3_N_Aut_N.atoms + AutAssimNH4_NH4_N_Aut_N.atoms - AutDecay_Aut_N_DOM_N.atoms;
Met_N.finalAtoms = Met_N.initialAtoms + MetAssimNO3_NO3_N_Met_N.atoms + MetAssimNH4_NH4_N_Met_N.atoms - MetDecay_Met_N_DOM_N.atoms;
CO2_C12.finalAtoms = CO2_C12.initialAtoms + MetMethaneOxid_CH4_C12_CO2_C12.atoms - AutAssimCO2_CO2_C12_Aut_C12.atoms;
CO2_C13.finalAtoms = CO2_C13.initialAtoms + MetMethaneOxid_CH4_C13_CO2_C13.atoms - AutAssimCO2_CO2_C13_Aut_C13.atoms;
CH4_C12.finalAtoms = CH4_C12.initialAtoms - MetMethaneOxid_CH4_C12_CO2_C12.atoms - MetAssimCH4_CH4_C12_Met_C12.atoms;
CH4_C13.finalAtoms = CH4_C13.initialAtoms - MetMethaneOxid_CH4_C13_CO2_C13.atoms - MetAssimCH4_CH4_C13_Met_C13.atoms;
DOM_C12.finalAtoms = DOM_C12.initialAtoms + AutDecay_Aut_C12_DOM_C12.atoms + MetDecay_Met_C12_DOM_C12.atoms;
DOM_C13.finalAtoms = DOM_C13.initialAtoms + AutDecay_Aut_C13_DOM_C13.atoms + MetDecay_Met_C13_DOM_C13.atoms;
Aut_C12.finalAtoms = Aut_C12.initialAtoms + AutAssimCO2_CO2_C12_Aut_C12.atoms - AutDecay_Aut_C12_DOM_C12.atoms;
Aut_C13.finalAtoms = Aut_C13.initialAtoms + AutAssimCO2_CO2_C13_Aut_C13.atoms - AutDecay_Aut_C13_DOM_C13.atoms;
Met_C12.finalAtoms = Met_C12.initialAtoms + MetAssimCH4_CH4_C12_Met_C12.atoms - MetDecay_Met_C12_DOM_C12.atoms;
Met_C13.finalAtoms = Met_C13.initialAtoms + MetAssimCH4_CH4_C13_Met_C13.atoms - MetDecay_Met_C13_DOM_C13.atoms;
CO2_O16.finalAtoms = CO2_O16.initialAtoms + MetMethaneOxid_O2_O16_CO2_O16.atoms - AutAssimCO2_CO2_O16_Aut_O16.atoms - AutAssimCO2_CO2_O16_Ox_O16.atoms;
CO2_O17.finalAtoms = CO2_O17.initialAtoms + MetMethaneOxid_O2_O17_CO2_O17.atoms - AutAssimCO2_CO2_O17_Aut_O17.atoms - AutAssimCO2_CO2_O17_Ox_O17.atoms;
CO2_O18.finalAtoms = CO2_O18.initialAtoms + MetMethaneOxid_O2_O18_CO2_O18.atoms - AutAssimCO2_CO2_O18_Aut_O18.atoms - AutAssimCO2_CO2_O18_Ox_O18.atoms;
NO3_O16.finalAtoms = NO3_O16.initialAtoms + AutNitrif_O2_O16_NO3_O16.atoms - AutAssimNO3_NO3_O16_Aut_O16.atoms - AutAssimNO3_NO3_O16_Ox_O16.atoms - MetAssimNO3_NO3_O16_Met_O16.atoms - MetAssimNO3_NO3_O16_Ox_O16.atoms;
NO3_O17.finalAtoms = NO3_O17.initialAtoms + AutNitrif_O2_O17_NO3_O17.atoms - AutAssimNO3_NO3_O17_Aut_O17.atoms - AutAssimNO3_NO3_O17_Ox_O17.atoms - MetAssimNO3_NO3_O17_Met_O17.atoms - MetAssimNO3_NO3_O17_Ox_O17.atoms;
NO3_O18.finalAtoms = NO3_O18.initialAtoms + AutNitrif_O2_O18_NO3_O18.atoms - AutAssimNO3_NO3_O18_Aut_O18.atoms - AutAssimNO3_NO3_O18_Ox_O18.atoms - MetAssimNO3_NO3_O18_Met_O18.atoms - MetAssimNO3_NO3_O18_Ox_O18.atoms;
O2_O16.finalAtoms = O2_O16.initialAtoms - AutNitrif_O2_O16_NO3_O16.atoms - AutNitrif_O2_O16_Ox_O16.atoms - MetMethaneOxid_O2_O16_CO2_O16.atoms - MetMethaneOxid_O2_O16_Ox_O16.atoms;
O2_O17.finalAtoms = O2_O17.initialAtoms - AutNitrif_O2_O17_NO3_O17.atoms - AutNitrif_O2_O17_Ox_O17.atoms - MetMethaneOxid_O2_O17_CO2_O17.atoms - MetMethaneOxid_O2_O17_Ox_O17.atoms;
O2_O18.finalAtoms = O2_O18.initialAtoms - AutNitrif_O2_O18_NO3_O18.atoms - AutNitrif_O2_O18_Ox_O18.atoms - MetMethaneOxid_O2_O18_CO2_O18.atoms - MetMethaneOxid_O2_O18_Ox_O18.atoms;
DOM_O16.finalAtoms = DOM_O16.initialAtoms + AutDecay_Aut_O16_DOM_O16.atoms + MetDecay_Met_O16_DOM_O16.atoms;
DOM_O17.finalAtoms = DOM_O17.initialAtoms + AutDecay_Aut_O17_DOM_O17.atoms + MetDecay_Met_O17_DOM_O17.atoms;
DOM_O18.finalAtoms = DOM_O18.initialAtoms + AutDecay_Aut_O18_DOM_O18.atoms + MetDecay_Met_O18_DOM_O18.atoms;
Ox_O16.finalAtoms = Ox_O16.initialAtoms + AutNitrif_O2_O16_Ox_O16.atoms + MetMethaneOxid_O2_O16_Ox_O16.atoms + AutAssimCO2_CO2_O16_Ox_O16.atoms + AutAssimNO3_NO3_O16_Ox_O16.atoms + MetAssimNO3_NO3_O16_Ox_O16.atoms;
Ox_O17.finalAtoms = Ox_O17.initialAtoms + AutNitrif_O2_O17_Ox_O17.atoms + MetMethaneOxid_O2_O17_Ox_O17.atoms + AutAssimCO2_CO2_O17_Ox_O17.atoms + AutAssimNO3_NO3_O17_Ox_O17.atoms + MetAssimNO3_NO3_O17_Ox_O17.atoms;
Ox_O18.finalAtoms = Ox_O18.initialAtoms + AutNitrif_O2_O18_Ox_O18.atoms + MetMethaneOxid_O2_O18_Ox_O18.atoms + AutAssimCO2_CO2_O18_Ox_O18.atoms + AutAssimNO3_NO3_O18_Ox_O18.atoms + MetAssimNO3_NO3_O18_Ox_O18.atoms;
Aut_O16.finalAtoms = Aut_O16.initialAtoms + AutAssimCO2_CO2_O16_Aut_O16.atoms + AutAssimNO3_NO3_O16_Aut_O16.atoms - AutDecay_Aut_O16_DOM_O16.atoms;
Aut_O17.finalAtoms = Aut_O17.initialAtoms + AutAssimCO2_CO2_O17_Aut_O17.atoms + AutAssimNO3_NO3_O17_Aut_O17.atoms - AutDecay_Aut_O17_DOM_O17.atoms;
Aut_O18.finalAtoms = Aut_O18.initialAtoms + AutAssimCO2_CO2_O18_Aut_O18.atoms + AutAssimNO3_NO3_O18_Aut_O18.atoms - AutDecay_Aut_O18_DOM_O18.atoms;
Met_O16.finalAtoms = Met_O16.initialAtoms + MetAssimNO3_NO3_O16_Met_O16.atoms - MetDecay_Met_O16_DOM_O16.atoms;
Met_O17.finalAtoms = Met_O17.initialAtoms + MetAssimNO3_NO3_O17_Met_O17.atoms - MetDecay_Met_O17_DOM_O17.atoms;
Met_O18.finalAtoms = Met_O18.initialAtoms + MetAssimNO3_NO3_O18_Met_O18.atoms - MetDecay_Met_O18_DOM_O18.atoms;
/***** Generally, Transfers out of each Pool must be less than the Pool.initialAtoms (Exprsn. 7) *****/
CO2_C12.initialAtoms >= AutAssimCO2_CO2_C12_Aut_C12.atoms;
CO2_C13.initialAtoms >= AutAssimCO2_CO2_C13_Aut_C13.atoms;
CO2_O16.initialAtoms >= AutAssimCO2_CO2_O16_Aut_O16.atoms + AutAssimCO2_CO2_O16_Ox_O16.atoms;
CO2_O17.initialAtoms >= AutAssimCO2_CO2_O17_Aut_O17.atoms + AutAssimCO2_CO2_O17_Ox_O17.atoms;
CO2_O18.initialAtoms >= AutAssimCO2_CO2_O18_Aut_O18.atoms + AutAssimCO2_CO2_O18_Ox_O18.atoms;
CH4_C12.initialAtoms >= MetMethaneOxid_CH4_C12_CO2_C12.atoms + MetAssimCH4_CH4_C12_Met_C12.atoms;
CH4_C13.initialAtoms >= MetMethaneOxid_CH4_C13_CO2_C13.atoms + MetAssimCH4_CH4_C13_Met_C13.atoms;
NH4_N.initialAtoms >= AutNitrif_NH4_N_NO3_N.atoms + AutAssimNH4_NH4_N_Aut_N.atoms + MetAssimNH4_NH4_N_Met_N.atoms;
NO3_N.initialAtoms >= AutAssimNO3_NO3_N_Aut_N.atoms + MetAssimNO3_NO3_N_Met_N.atoms;
NO3_O16.initialAtoms >= AutAssimNO3_NO3_O16_Aut_O16.atoms + AutAssimNO3_NO3_O16_Ox_O16.atoms + MetAssimNO3_NO3_O16_Met_O16.atoms + MetAssimNO3_NO3_O16_Ox_O16.atoms;
NO3_O17.initialAtoms >= AutAssimNO3_NO3_O17_Aut_O17.atoms + AutAssimNO3_NO3_O17_Ox_O17.atoms + MetAssimNO3_NO3_O17_Met_O17.atoms + MetAssimNO3_NO3_O17_Ox_O17.atoms;
NO3_O18.initialAtoms >= AutAssimNO3_NO3_O18_Aut_O18.atoms + AutAssimNO3_NO3_O18_Ox_O18.atoms + MetAssimNO3_NO3_O18_Met_O18.atoms + MetAssimNO3_NO3_O18_Ox_O18.atoms;
O2_O16.initialAtoms >= AutNitrif_O2_O16_NO3_O16.atoms + AutNitrif_O2_O16_Ox_O16.atoms + MetMethaneOxid_O2_O16_CO2_O16.atoms + MetMethaneOxid_O2_O16_Ox_O16.atoms;
O2_O17.initialAtoms >= AutNitrif_O2_O17_NO3_O17.atoms + AutNitrif_O2_O17_Ox_O17.atoms + MetMethaneOxid_O2_O17_CO2_O17.atoms + MetMethaneOxid_O2_O17_Ox_O17.atoms;
O2_O18.initialAtoms >= AutNitrif_O2_O18_NO3_O18.atoms + AutNitrif_O2_O18_Ox_O18.atoms + MetMethaneOxid_O2_O18_CO2_O18.atoms + MetMethaneOxid_O2_O18_Ox_O18.atoms;
Aut_N.initialAtoms >= AutDecay_Aut_N_DOM_N.atoms;
Aut_C12.initialAtoms >= AutDecay_Aut_C12_DOM_C12.atoms;
Aut_C13.initialAtoms >= AutDecay_Aut_C13_DOM_C13.atoms;
Aut_O16.initialAtoms >= AutDecay_Aut_O16_DOM_O16.atoms;
Aut_O17.initialAtoms >= AutDecay_Aut_O17_DOM_O17.atoms;
Aut_O18.initialAtoms >= AutDecay_Aut_O18_DOM_O18.atoms;
Met_N.initialAtoms >= MetDecay_Met_N_DOM_N.atoms;
Met_C12.initialAtoms >= MetDecay_Met_C12_DOM_C12.atoms;
Met_C13.initialAtoms >= MetDecay_Met_C13_DOM_C13.atoms;
Met_O16.initialAtoms >= MetDecay_Met_O16_DOM_O16.atoms;
Met_O17.initialAtoms >= MetDecay_Met_O17_DOM_O17.atoms;
Met_O18.initialAtoms >= MetDecay_Met_O18_DOM_O18.atoms;
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

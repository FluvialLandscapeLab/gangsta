/******************** OPTIMIZATION ********************/
/***** Optimization function (Exprsn. 1) *****/
MAX: Total_Biomass.finalMolecules;
/***** The Total_Biomass.finalMolecules equals the sum of the Organism.finalMolecules for each organism type (Exprsn. 2) *****/
Total_Biomass.finalMolecules = Het.finalMolecules + Aut.finalMolecules + Met.finalMolecules;
/******************** INITIALIZATION ********************/
/***** Set FiniteCompound.initialMolecules & InfiniteCompound.initialMolecules *****/
Het.initialMolecules = 0;
Aut.initialMolecules = 0;
Met.initialMolecules = 0;
DOM.initialMolecules = 0;
CH4.initialMolecules = 0;
NH4.initialMolecules = 0;
NO3.initialMolecules = 0;
CO2.initialMolecules = 0;
N2.initialMolecules = 0;
/***** For each InfiniteCompound, remove constraints on InfiniteCompound.finalMolecules *****/
/***** For each Pool associated with an InfiniteCompound, remove constraints on Pool.finalAtoms *****/
/***** For each Process, the Process.netEnergy is unconstrained *****/
-Inf < HetAssimDOM.netEnergy < +Inf;
-Inf < AutAssimCO2.netEnergy < +Inf;
-Inf < MetAssimCH4.netEnergy < +Inf;
-Inf < HetAssimNO3.netEnergy < +Inf;
-Inf < AutAssimNO3.netEnergy < +Inf;
-Inf < MetAssimNO3.netEnergy < +Inf;
-Inf < HetAssimNH4.netEnergy < +Inf;
-Inf < AutAssimNH4.netEnergy < +Inf;
-Inf < MetAssimNH4.netEnergy < +Inf;
-Inf < HetAerobicofBiomass.netEnergy < +Inf;
-Inf < HetAerobicofDOM.netEnergy < +Inf;
-Inf < HetDenitofBiomass.netEnergy < +Inf;
-Inf < HetDenitofDOM.netEnergy < +Inf;
-Inf < HetMethanogenesisofBiomass.netEnergy < +Inf;
-Inf < HetMethanogenesisofDOM.netEnergy < +Inf;
-Inf < AutNitrif.netEnergy < +Inf;
-Inf < MetMethaneOxid.netEnergy < +Inf;
-Inf < HetDecay.netEnergy < +Inf;
-Inf < AutDecay.netEnergy < +Inf;
-Inf < MetDecay.netEnergy < +Inf;
/***** For each Organism type, the Organism.respirationEnergy is unconstrained *****/
-Inf < Het.respirationEnergy < +Inf;
-Inf < Aut.respirationEnergy < +Inf;
-Inf < Met.respirationEnergy < +Inf;
/******************** STOICHIOMETRY ********************/
/***** For each elemental Pool, the Pool.initialAtoms must conform to Compound stoichiometry (Exprsn. 3) *****/
Het_C.initialAtoms = 1 Het.initialMolecules;
Het_N.initialAtoms = 0.150943396226415 Het.initialMolecules;
Aut_C.initialAtoms = 1 Aut.initialMolecules;
Aut_N.initialAtoms = 0.150943396226415 Aut.initialMolecules;
Met_C.initialAtoms = 1 Met.initialMolecules;
Met_N.initialAtoms = 0.150943396226415 Met.initialMolecules;
DOM_C.initialAtoms = 1 DOM.initialMolecules;
DOM_N.initialAtoms = 0.150943396226415 DOM.initialMolecules;
CH4_C.initialAtoms = 1 CH4.initialMolecules;
NH4_N.initialAtoms = 1 NH4.initialMolecules;
NO3_N.initialAtoms = 1 NO3.initialMolecules;
CO2_C.initialAtoms = 1 CO2.initialMolecules;
N2_N.initialAtoms = 2 N2.initialMolecules;
/***** For each elemental Pool, the Pool.finalAtoms must conform to Compound stoichiometry (Exprsn. 4) *****/
Het_C.finalAtoms = 1 Het.finalMolecules;
Het_N.finalAtoms = 0.150943396226415 Het.finalMolecules;
Aut_C.finalAtoms = 1 Aut.finalMolecules;
Aut_N.finalAtoms = 0.150943396226415 Aut.finalMolecules;
Met_C.finalAtoms = 1 Met.finalMolecules;
Met_N.finalAtoms = 0.150943396226415 Met.finalMolecules;
DOM_C.finalAtoms = 1 DOM.finalMolecules;
DOM_N.finalAtoms = 0.150943396226415 DOM.finalMolecules;
CH4_C.finalAtoms = 1 CH4.finalMolecules;
NH4_N.finalAtoms = 1 NH4.finalMolecules;
NO3_N.finalAtoms = 1 NO3.finalMolecules;
CO2_C.finalAtoms = 1 CO2.finalMolecules;
N2_N.finalAtoms = 2 N2.finalMolecules;
/******************** MOLAR TRANSFER ********************/
/***** The number of atoms transferred during each Process must be in accordance with the stoichiometry and chemical affinity of the Process (Exprsn. 5) *****/
HetAssimDOM.netEnergy = -0.000432 HetAssimDOM_DOM_C_Het_C.atoms;
HetAssimDOM.netEnergy = -0.002862 HetAssimDOM_DOM_N_Het_N.atoms;
AutAssimCO2.netEnergy = -0.0035 AutAssimCO2_CO2_C_Aut_C.atoms;
MetAssimCH4.netEnergy = -0.00109 MetAssimCH4_CH4_C_Met_C.atoms;
HetAssimNO3.netEnergy = -0.000155 HetAssimNO3_NO3_N_Het_N.atoms;
AutAssimNO3.netEnergy = -0.000155 AutAssimNO3_NO3_N_Aut_N.atoms;
MetAssimNO3.netEnergy = -0.000155 MetAssimNO3_NO3_N_Met_N.atoms;
HetAssimNH4.netEnergy = -3.18e-05 HetAssimNH4_NH4_N_Het_N.atoms;
AutAssimNH4.netEnergy = -3.18e-05 AutAssimNH4_NH4_N_Aut_N.atoms;
MetAssimNH4.netEnergy = -3.18e-05 MetAssimNH4_NH4_N_Met_N.atoms;
HetAerobicofBiomass.netEnergy = 0.000437 HetAerobicofBiomass_Het_C_CO2_C.atoms;
HetAerobicofBiomass.netEnergy = 0.002895125 HetAerobicofBiomass_Het_N_NH4_N.atoms;
HetAerobicofDOM.netEnergy = 0.000437 HetAerobicofDOM_DOM_C_CO2_C.atoms;
HetAerobicofDOM.netEnergy = 0.002895125 HetAerobicofDOM_DOM_N_NH4_N.atoms;
HetDenitofBiomass.netEnergy = 0.0004102 HetDenitofBiomass_Het_C_CO2_C.atoms;
HetDenitofBiomass.netEnergy = 0.002717575 HetDenitofBiomass_Het_N_NH4_N.atoms;
HetDenitofBiomass.netEnergy = 0.00051275 HetDenitofBiomass_NO3_N_N2_N.atoms;
HetDenitofDOM.netEnergy = 0.0004102 HetDenitofDOM_DOM_C_CO2_C.atoms;
HetDenitofDOM.netEnergy = 0.002717575 HetDenitofDOM_DOM_N_NH4_N.atoms;
HetDenitofDOM.netEnergy = 0.00051275 HetDenitofDOM_NO3_N_N2_N.atoms;
HetMethanogenesisofBiomass.netEnergy = 5.6e-05 HetMethanogenesisofBiomass_Het_C_CO2_C.atoms;
HetMethanogenesisofBiomass.netEnergy = 5.6e-05 HetMethanogenesisofBiomass_Het_C_CH4_C.atoms;
HetMethanogenesisofBiomass.netEnergy = 0.0001855 HetMethanogenesisofBiomass_Het_N_NH4_N.atoms;
HetMethanogenesisofDOM.netEnergy = 5.6e-05 HetMethanogenesisofDOM_DOM_C_CO2_C.atoms;
HetMethanogenesisofDOM.netEnergy = 5.6e-05 HetMethanogenesisofDOM_DOM_C_CH4_C.atoms;
HetMethanogenesisofDOM.netEnergy = 0.0001855 HetMethanogenesisofDOM_DOM_N_NH4_N.atoms;
AutNitrif.netEnergy = 0.0003485 AutNitrif_NH4_N_NO3_N.atoms;
MetMethaneOxid.netEnergy = 0.000818 MetMethaneOxid_CH4_C_CO2_C.atoms;
/******************** ELEMENTAL MOLAR BALANCE ********************/
/***** Elemental molar balance must be conserved (Exprsn. 6) *****/
Het_C.finalAtoms = Het_C.initialAtoms + HetAssimDOM_DOM_C_Het_C.atoms - HetAerobicofBiomass_Het_C_CO2_C.atoms - HetDenitofBiomass_Het_C_CO2_C.atoms - HetMethanogenesisofBiomass_Het_C_CO2_C.atoms - HetMethanogenesisofBiomass_Het_C_CH4_C.atoms - HetDecay_Het_C_DOM_C.atoms;
Het_N.finalAtoms = Het_N.initialAtoms + HetAssimDOM_DOM_N_Het_N.atoms + HetAssimNO3_NO3_N_Het_N.atoms + HetAssimNH4_NH4_N_Het_N.atoms - HetAerobicofBiomass_Het_N_NH4_N.atoms - HetDenitofBiomass_Het_N_NH4_N.atoms - HetMethanogenesisofBiomass_Het_N_NH4_N.atoms - HetDecay_Het_N_DOM_N.atoms;
Aut_C.finalAtoms = Aut_C.initialAtoms + AutAssimCO2_CO2_C_Aut_C.atoms - AutDecay_Aut_C_DOM_C.atoms;
Aut_N.finalAtoms = Aut_N.initialAtoms + AutAssimNO3_NO3_N_Aut_N.atoms + AutAssimNH4_NH4_N_Aut_N.atoms - AutDecay_Aut_N_DOM_N.atoms;
Met_C.finalAtoms = Met_C.initialAtoms + MetAssimCH4_CH4_C_Met_C.atoms - MetDecay_Met_C_DOM_C.atoms;
Met_N.finalAtoms = Met_N.initialAtoms + MetAssimNO3_NO3_N_Met_N.atoms + MetAssimNH4_NH4_N_Met_N.atoms - MetDecay_Met_N_DOM_N.atoms;
DOM_C.finalAtoms = DOM_C.initialAtoms + HetDecay_Het_C_DOM_C.atoms + AutDecay_Aut_C_DOM_C.atoms + MetDecay_Met_C_DOM_C.atoms - HetAssimDOM_DOM_C_Het_C.atoms - HetAerobicofDOM_DOM_C_CO2_C.atoms - HetDenitofDOM_DOM_C_CO2_C.atoms - HetMethanogenesisofDOM_DOM_C_CO2_C.atoms - HetMethanogenesisofDOM_DOM_C_CH4_C.atoms;
DOM_N.finalAtoms = DOM_N.initialAtoms + HetDecay_Het_N_DOM_N.atoms + AutDecay_Aut_N_DOM_N.atoms + MetDecay_Met_N_DOM_N.atoms - HetAssimDOM_DOM_N_Het_N.atoms - HetAerobicofDOM_DOM_N_NH4_N.atoms - HetDenitofDOM_DOM_N_NH4_N.atoms - HetMethanogenesisofDOM_DOM_N_NH4_N.atoms;
CH4_C.finalAtoms = CH4_C.initialAtoms + HetMethanogenesisofBiomass_Het_C_CH4_C.atoms + HetMethanogenesisofDOM_DOM_C_CH4_C.atoms - MetAssimCH4_CH4_C_Met_C.atoms - MetMethaneOxid_CH4_C_CO2_C.atoms;
NH4_N.finalAtoms = NH4_N.initialAtoms + HetAerobicofBiomass_Het_N_NH4_N.atoms + HetAerobicofDOM_DOM_N_NH4_N.atoms + HetDenitofBiomass_Het_N_NH4_N.atoms + HetDenitofDOM_DOM_N_NH4_N.atoms + HetMethanogenesisofBiomass_Het_N_NH4_N.atoms + HetMethanogenesisofDOM_DOM_N_NH4_N.atoms - HetAssimNH4_NH4_N_Het_N.atoms - AutAssimNH4_NH4_N_Aut_N.atoms - MetAssimNH4_NH4_N_Met_N.atoms - AutNitrif_NH4_N_NO3_N.atoms;
NO3_N.finalAtoms = NO3_N.initialAtoms + AutNitrif_NH4_N_NO3_N.atoms - HetAssimNO3_NO3_N_Het_N.atoms - AutAssimNO3_NO3_N_Aut_N.atoms - MetAssimNO3_NO3_N_Met_N.atoms - HetDenitofBiomass_NO3_N_N2_N.atoms - HetDenitofDOM_NO3_N_N2_N.atoms;
CO2_C.finalAtoms = CO2_C.initialAtoms + HetAerobicofBiomass_Het_C_CO2_C.atoms + HetAerobicofDOM_DOM_C_CO2_C.atoms + HetDenitofBiomass_Het_C_CO2_C.atoms + HetDenitofDOM_DOM_C_CO2_C.atoms + HetMethanogenesisofBiomass_Het_C_CO2_C.atoms + HetMethanogenesisofDOM_DOM_C_CO2_C.atoms + MetMethaneOxid_CH4_C_CO2_C.atoms - AutAssimCO2_CO2_C_Aut_C.atoms;
N2_N.finalAtoms = N2_N.initialAtoms + HetDenitofBiomass_NO3_N_N2_N.atoms + HetDenitofDOM_NO3_N_N2_N.atoms;
/***** Generally, Transfers out of each Pool must be less than the Pool.initialAtoms (Exprsn. 7) *****/
Het_C.initialAtoms >= HetAerobicofBiomass_Het_C_CO2_C.atoms + HetDenitofBiomass_Het_C_CO2_C.atoms + HetMethanogenesisofBiomass_Het_C_CO2_C.atoms + HetMethanogenesisofBiomass_Het_C_CH4_C.atoms + HetDecay_Het_C_DOM_C.atoms;
Het_N.initialAtoms >= HetAerobicofBiomass_Het_N_NH4_N.atoms + HetDenitofBiomass_Het_N_NH4_N.atoms + HetMethanogenesisofBiomass_Het_N_NH4_N.atoms + HetDecay_Het_N_DOM_N.atoms;
Aut_C.initialAtoms >= AutDecay_Aut_C_DOM_C.atoms;
Aut_N.initialAtoms >= AutDecay_Aut_N_DOM_N.atoms;
Met_C.initialAtoms >= MetDecay_Met_C_DOM_C.atoms;
Met_N.initialAtoms >= MetDecay_Met_N_DOM_N.atoms;
DOM_C.initialAtoms >= HetAssimDOM_DOM_C_Het_C.atoms + HetAerobicofDOM_DOM_C_CO2_C.atoms + HetDenitofDOM_DOM_C_CO2_C.atoms + HetMethanogenesisofDOM_DOM_C_CO2_C.atoms + HetMethanogenesisofDOM_DOM_C_CH4_C.atoms;
DOM_N.initialAtoms >= HetAssimDOM_DOM_N_Het_N.atoms + HetAerobicofDOM_DOM_N_NH4_N.atoms + HetDenitofDOM_DOM_N_NH4_N.atoms + HetMethanogenesisofDOM_DOM_N_NH4_N.atoms;
CH4_C.initialAtoms >= MetAssimCH4_CH4_C_Met_C.atoms + MetMethaneOxid_CH4_C_CO2_C.atoms;
NH4_N.initialAtoms >= AutAssimNH4_NH4_N_Aut_N.atoms + MetAssimNH4_NH4_N_Met_N.atoms + AutNitrif_NH4_N_NO3_N.atoms;
NO3_N.initialAtoms >= HetAssimNO3_NO3_N_Het_N.atoms + AutAssimNO3_NO3_N_Aut_N.atoms + MetAssimNO3_NO3_N_Met_N.atoms + HetDenitofBiomass_NO3_N_N2_N.atoms + HetDenitofDOM_NO3_N_N2_N.atoms;
CO2_C.initialAtoms >= AutAssimCO2_CO2_C_Aut_C.atoms;
/******************** ENERGETIC CONSTRAINTS ********************/
/***** Maintenance respiration energy is proportional to the biomass of each organism type (Exprsn. 8) *****/
Het.respirationEnergy = -6.792e-05 Het.finalMolecules;
Aut.respirationEnergy = -6.792e-05 Aut.finalMolecules;
Met.respirationEnergy = -6.792e-05 Met.finalMolecules;
/***** For each organism type, the total energy associated with assimilatory and dissimilatory Processes must equal the sum of the Process.netEnergy for all Processes (Exprsn. 9) *****/
Het.totalProcessEnergy = HetAssimDOM.netEnergy + HetAssimNO3.netEnergy + HetAssimNH4.netEnergy + HetAerobicofBiomass.netEnergy + HetAerobicofDOM.netEnergy + HetDenitofBiomass.netEnergy + HetDenitofDOM.netEnergy + HetMethanogenesisofBiomass.netEnergy + HetMethanogenesisofDOM.netEnergy;
Aut.totalProcessEnergy = AutAssimCO2.netEnergy + AutAssimNO3.netEnergy + AutAssimNH4.netEnergy + AutNitrif.netEnergy;
Met.totalProcessEnergy = MetAssimCH4.netEnergy + MetAssimNO3.netEnergy + MetAssimNH4.netEnergy + MetMethaneOxid.netEnergy;
/***** Organisms must meet their respiratory energy demands (Exprsn. 10) *****/
0 = Het.respirationEnergy + Het.totalProcessEnergy;
0 = Aut.respirationEnergy + Aut.totalProcessEnergy;
0 = Met.respirationEnergy + Met.totalProcessEnergy;

/******************** OPTIMIZATION ********************/
/***** Optimization function (Exprsn. 1) *****/
MAX: Total_Biomass.finalAmountMols;
/***** The mols of total biomass equals the sum of the final mols of biomass in each organism type (Exprsn. 2) *****/
Total_Biomass.finalAmountMols = Het.finalAmountMols + Aut.finalAmountMols + Met.finalAmountMols;
/******************** INITIALIZATION ********************/
/***** Set initial mols of compounds (Exprsn. 3) & Initial mols of source/sink compounds must equal 0 (Exprsn. 4) *****/
Het.initialAmountMols = 0;
Aut.initialAmountMols = 0;
Met.initialAmountMols = 0;
DOM.initialAmountMols = 0;
CH4.initialAmountMols = 0;
NH4.initialAmountMols = 0;
NO3.initialAmountMols = 0;
CO2.initialAmountMols = 0;
N2.initialAmountMols = 0;
/***** Remove constraints on final mols of source/sink compounds (Exprsn. 5) *****/
/***** Remove constraints on final mols of source/sink pools (Exprsn. 6) *****/
/***** The energies associated with decay, assimilatory, and dissimilatory processes are unconstrained here, but depend on parameterization of process energies (Exprsn. 7) *****/
-Inf < HetAssimDOM.energyJoules < +Inf;
-Inf < AutAssimCO2.energyJoules < +Inf;
-Inf < MetAssimCH4.energyJoules < +Inf;
-Inf < HetAssimNO3.energyJoules < +Inf;
-Inf < AutAssimNO3.energyJoules < +Inf;
-Inf < MetAssimNO3.energyJoules < +Inf;
-Inf < HetAssimNH4.energyJoules < +Inf;
-Inf < AutAssimNH4.energyJoules < +Inf;
-Inf < MetAssimNH4.energyJoules < +Inf;
-Inf < HetDenitofBiomass.energyJoules < +Inf;
-Inf < HetDenitofDOM.energyJoules < +Inf;
-Inf < HetMethanogenesisofBiomass.energyJoules < +Inf;
-Inf < HetMethanogenesisofDOM.energyJoules < +Inf;
-Inf < HetDecay.energyJoules < +Inf;
-Inf < AutDecay.energyJoules < +Inf;
-Inf < MetDecay.energyJoules < +Inf;
/***** The respiration rate associated with each organism type is a user-defined constant; therefore, the respiration energy associated with each organism type must be unconstrained. (Exprsn. 8) *****/
-Inf < Het.respirationEnergyJoules < +Inf;
-Inf < Aut.respirationEnergyJoules < +Inf;
-Inf < Met.respirationEnergyJoules < +Inf;
/******************** STOICHIOMETRY ********************/
/***** The initial mols of elemental pools must conform to compound stoichiometry (Exprsn. 9) *****/
Het_C.initialAmountMols = 1 Het.initialAmountMols;
Het_N.initialAmountMols = 0.150943396226415 Het.initialAmountMols;
Aut_C.initialAmountMols = 1 Aut.initialAmountMols;
Aut_N.initialAmountMols = 0.150943396226415 Aut.initialAmountMols;
Met_C.initialAmountMols = 1 Met.initialAmountMols;
Met_N.initialAmountMols = 0.150943396226415 Met.initialAmountMols;
DOM_C.initialAmountMols = 1 DOM.initialAmountMols;
DOM_N.initialAmountMols = 0.150943396226415 DOM.initialAmountMols;
CH4_C.initialAmountMols = 1 CH4.initialAmountMols;
NH4_N.initialAmountMols = 1 NH4.initialAmountMols;
NO3_N.initialAmountMols = 1 NO3.initialAmountMols;
CO2_C.initialAmountMols = 1 CO2.initialAmountMols;
N2_N.initialAmountMols = 2 N2.initialAmountMols;
/***** The final mols of elemental pools must conform to compound stoichiometry (Exprsn. 10) *****/
Het_C.finalAmountMols = 1 Het.finalAmountMols;
Het_N.finalAmountMols = 0.150943396226415 Het.finalAmountMols;
Aut_C.finalAmountMols = 1 Aut.finalAmountMols;
Aut_N.finalAmountMols = 0.150943396226415 Aut.finalAmountMols;
Met_C.finalAmountMols = 1 Met.finalAmountMols;
Met_N.finalAmountMols = 0.150943396226415 Met.finalAmountMols;
DOM_C.finalAmountMols = 1 DOM.finalAmountMols;
DOM_N.finalAmountMols = 0.150943396226415 DOM.finalAmountMols;
CH4_C.finalAmountMols = 1 CH4.finalAmountMols;
NH4_N.finalAmountMols = 1 NH4.finalAmountMols;
NO3_N.finalAmountMols = 1 NO3.finalAmountMols;
CO2_C.finalAmountMols = 1 CO2.finalAmountMols;
N2_N.finalAmountMols = 2 N2.finalAmountMols;
/******************** MOLAR TRANSFER ********************/
/***** The energy and molar transfer for each process must be in accordance with the stoichiometry and free energy yield of the process (Exprsn. 11) *****/
HetAssimDOM.energyJoules = -0.000432 HetAssimDOM_DOM_C_Het_C.amountMolsTransfer;
HetAssimDOM.energyJoules = -0.002862 HetAssimDOM_DOM_N_Het_N.amountMolsTransfer;
AutAssimCO2.energyJoules = -0.035 AutAssimCO2_CO2_C_Aut_C.amountMolsTransfer;
MetAssimCH4.energyJoules = -0.00109 MetAssimCH4_CH4_C_Met_C.amountMolsTransfer;
HetAssimNO3.energyJoules = -0.000155 HetAssimNO3_NO3_N_Het_N.amountMolsTransfer;
AutAssimNO3.energyJoules = -0.000155 AutAssimNO3_NO3_N_Aut_N.amountMolsTransfer;
MetAssimNO3.energyJoules = -0.000155 MetAssimNO3_NO3_N_Met_N.amountMolsTransfer;
HetAssimNH4.energyJoules = -3.18e-05 HetAssimNH4_NH4_N_Het_N.amountMolsTransfer;
AutAssimNH4.energyJoules = -3.18e-05 AutAssimNH4_NH4_N_Aut_N.amountMolsTransfer;
MetAssimNH4.energyJoules = -3.18e-05 MetAssimNH4_NH4_N_Met_N.amountMolsTransfer;
HetDenitofBiomass.energyJoules = 0.0004102 HetDenitofBiomass_Het_C_CO2_C.amountMolsTransfer;
HetDenitofBiomass.energyJoules = 0.002717575 HetDenitofBiomass_Het_N_NH4_N.amountMolsTransfer;
HetDenitofBiomass.energyJoules = 0.00051275 HetDenitofBiomass_NO3_N_N2_N.amountMolsTransfer;
HetDenitofDOM.energyJoules = 0.0004102 HetDenitofDOM_DOM_C_CO2_C.amountMolsTransfer;
HetDenitofDOM.energyJoules = 0.002717575 HetDenitofDOM_DOM_N_NH4_N.amountMolsTransfer;
HetDenitofDOM.energyJoules = 0.00051275 HetDenitofDOM_NO3_N_N2_N.amountMolsTransfer;
HetMethanogenesisofBiomass.energyJoules = 5.6e-05 HetMethanogenesisofBiomass_Het_C_CO2_C.amountMolsTransfer;
HetMethanogenesisofBiomass.energyJoules = 5.6e-05 HetMethanogenesisofBiomass_Het_C_CH4_C.amountMolsTransfer;
HetMethanogenesisofBiomass.energyJoules = 0.0001855 HetMethanogenesisofBiomass_Het_N_NH4_N.amountMolsTransfer;
HetMethanogenesisofDOM.energyJoules = 5.6e-05 HetMethanogenesisofDOM_DOM_C_CO2_C.amountMolsTransfer;
HetMethanogenesisofDOM.energyJoules = 5.6e-05 HetMethanogenesisofDOM_DOM_C_CH4_C.amountMolsTransfer;
HetMethanogenesisofDOM.energyJoules = 0.0001855 HetMethanogenesisofDOM_DOM_N_NH4_N.amountMolsTransfer;
/******************** ELEMENTAL MOLAR BALANCE ********************/
/***** The mols of each elemental pool must be conserved (Exprsn. 12) *****/
Het_C.finalAmountMols = Het_C.initialAmountMols + HetAssimDOM_DOM_C_Het_C.amountMolsTransfer - HetDenitofBiomass_Het_C_CO2_C.amountMolsTransfer - HetMethanogenesisofBiomass_Het_C_CO2_C.amountMolsTransfer - HetMethanogenesisofBiomass_Het_C_CH4_C.amountMolsTransfer - HetDecay_Het_C_DOM_C.amountMolsTransfer;
Het_N.finalAmountMols = Het_N.initialAmountMols + HetAssimDOM_DOM_N_Het_N.amountMolsTransfer + HetAssimNO3_NO3_N_Het_N.amountMolsTransfer + HetAssimNH4_NH4_N_Het_N.amountMolsTransfer - HetDenitofBiomass_Het_N_NH4_N.amountMolsTransfer - HetMethanogenesisofBiomass_Het_N_NH4_N.amountMolsTransfer - HetDecay_Het_N_DOM_N.amountMolsTransfer;
Aut_C.finalAmountMols = Aut_C.initialAmountMols + AutAssimCO2_CO2_C_Aut_C.amountMolsTransfer - AutDecay_Aut_C_DOM_C.amountMolsTransfer;
Aut_N.finalAmountMols = Aut_N.initialAmountMols + AutAssimNO3_NO3_N_Aut_N.amountMolsTransfer + AutAssimNH4_NH4_N_Aut_N.amountMolsTransfer - AutDecay_Aut_N_DOM_N.amountMolsTransfer;
Met_C.finalAmountMols = Met_C.initialAmountMols + MetAssimCH4_CH4_C_Met_C.amountMolsTransfer - MetDecay_Met_C_DOM_C.amountMolsTransfer;
Met_N.finalAmountMols = Met_N.initialAmountMols + MetAssimNO3_NO3_N_Met_N.amountMolsTransfer + MetAssimNH4_NH4_N_Met_N.amountMolsTransfer - MetDecay_Met_N_DOM_N.amountMolsTransfer;
DOM_C.finalAmountMols = DOM_C.initialAmountMols + HetDecay_Het_C_DOM_C.amountMolsTransfer + AutDecay_Aut_C_DOM_C.amountMolsTransfer + MetDecay_Met_C_DOM_C.amountMolsTransfer - HetAssimDOM_DOM_C_Het_C.amountMolsTransfer - HetDenitofDOM_DOM_C_CO2_C.amountMolsTransfer - HetMethanogenesisofDOM_DOM_C_CO2_C.amountMolsTransfer - HetMethanogenesisofDOM_DOM_C_CH4_C.amountMolsTransfer;
DOM_N.finalAmountMols = DOM_N.initialAmountMols + HetDecay_Het_N_DOM_N.amountMolsTransfer + AutDecay_Aut_N_DOM_N.amountMolsTransfer + MetDecay_Met_N_DOM_N.amountMolsTransfer - HetAssimDOM_DOM_N_Het_N.amountMolsTransfer - HetDenitofDOM_DOM_N_NH4_N.amountMolsTransfer - HetMethanogenesisofDOM_DOM_N_NH4_N.amountMolsTransfer;
CH4_C.finalAmountMols = CH4_C.initialAmountMols + HetMethanogenesisofBiomass_Het_C_CH4_C.amountMolsTransfer + HetMethanogenesisofDOM_DOM_C_CH4_C.amountMolsTransfer - MetAssimCH4_CH4_C_Met_C.amountMolsTransfer;
NH4_N.finalAmountMols = NH4_N.initialAmountMols + HetDenitofBiomass_Het_N_NH4_N.amountMolsTransfer + HetDenitofDOM_DOM_N_NH4_N.amountMolsTransfer + HetMethanogenesisofBiomass_Het_N_NH4_N.amountMolsTransfer + HetMethanogenesisofDOM_DOM_N_NH4_N.amountMolsTransfer - HetAssimNH4_NH4_N_Het_N.amountMolsTransfer - AutAssimNH4_NH4_N_Aut_N.amountMolsTransfer - MetAssimNH4_NH4_N_Met_N.amountMolsTransfer;
NO3_N.finalAmountMols = NO3_N.initialAmountMols - HetAssimNO3_NO3_N_Het_N.amountMolsTransfer - AutAssimNO3_NO3_N_Aut_N.amountMolsTransfer - MetAssimNO3_NO3_N_Met_N.amountMolsTransfer - HetDenitofBiomass_NO3_N_N2_N.amountMolsTransfer - HetDenitofDOM_NO3_N_N2_N.amountMolsTransfer;
CO2_C.finalAmountMols = CO2_C.initialAmountMols + HetDenitofBiomass_Het_C_CO2_C.amountMolsTransfer + HetDenitofDOM_DOM_C_CO2_C.amountMolsTransfer + HetMethanogenesisofBiomass_Het_C_CO2_C.amountMolsTransfer + HetMethanogenesisofDOM_DOM_C_CO2_C.amountMolsTransfer - AutAssimCO2_CO2_C_Aut_C.amountMolsTransfer;
N2_N.finalAmountMols = N2_N.initialAmountMols + HetDenitofBiomass_NO3_N_N2_N.amountMolsTransfer + HetDenitofDOM_NO3_N_N2_N.amountMolsTransfer;
/***** Generally, molar transfers from each pool must be less than the initial mols of the pool (Exprsn. 13) *****/
Het_C.initialAmountMols >= HetDenitofBiomass_Het_C_CO2_C.amountMolsTransfer + HetMethanogenesisofBiomass_Het_C_CO2_C.amountMolsTransfer + HetMethanogenesisofBiomass_Het_C_CH4_C.amountMolsTransfer + HetDecay_Het_C_DOM_C.amountMolsTransfer;
Het_N.initialAmountMols >= HetDenitofBiomass_Het_N_NH4_N.amountMolsTransfer + HetMethanogenesisofBiomass_Het_N_NH4_N.amountMolsTransfer + HetDecay_Het_N_DOM_N.amountMolsTransfer;
Aut_C.initialAmountMols >= AutDecay_Aut_C_DOM_C.amountMolsTransfer;
Aut_N.initialAmountMols >= AutDecay_Aut_N_DOM_N.amountMolsTransfer;
Met_C.initialAmountMols >= MetDecay_Met_C_DOM_C.amountMolsTransfer;
Met_N.initialAmountMols >= MetDecay_Met_N_DOM_N.amountMolsTransfer;
DOM_C.initialAmountMols >= HetAssimDOM_DOM_C_Het_C.amountMolsTransfer + HetDenitofDOM_DOM_C_CO2_C.amountMolsTransfer + HetMethanogenesisofDOM_DOM_C_CO2_C.amountMolsTransfer + HetMethanogenesisofDOM_DOM_C_CH4_C.amountMolsTransfer;
DOM_N.initialAmountMols >= HetAssimDOM_DOM_N_Het_N.amountMolsTransfer + HetDenitofDOM_DOM_N_NH4_N.amountMolsTransfer + HetMethanogenesisofDOM_DOM_N_NH4_N.amountMolsTransfer;
CH4_C.initialAmountMols >= MetAssimCH4_CH4_C_Met_C.amountMolsTransfer;
NH4_N.initialAmountMols >= AutAssimNH4_NH4_N_Aut_N.amountMolsTransfer + MetAssimNH4_NH4_N_Met_N.amountMolsTransfer;
NO3_N.initialAmountMols >= HetAssimNO3_NO3_N_Het_N.amountMolsTransfer + AutAssimNO3_NO3_N_Aut_N.amountMolsTransfer + MetAssimNO3_NO3_N_Met_N.amountMolsTransfer + HetDenitofBiomass_NO3_N_N2_N.amountMolsTransfer + HetDenitofDOM_NO3_N_N2_N.amountMolsTransfer;
CO2_C.initialAmountMols >= AutAssimCO2_CO2_C_Aut_C.amountMolsTransfer;
/******************** ENERGETIC CONSTRAINTS ********************/
/***** Maintenance respiration energy must be proportional to mols of biomass (Exprsn. 14) *****/
Het.respirationEnergyJoules = -6.792e-05 Het.finalAmountMols;
Aut.respirationEnergyJoules = -6.792e-05 Aut.finalAmountMols;
Met.respirationEnergyJoules = -6.792e-05 Met.finalAmountMols;
/***** The net energy associated with assimilatory and dissimilatory processes must equal the summed energy of each assimilatory and dissimilatory process  (Exprsn. 15) *****/
Het.energyJoules = HetAssimDOM.energyJoules + HetAssimNO3.energyJoules + HetAssimNH4.energyJoules + HetDenitofBiomass.energyJoules + HetDenitofDOM.energyJoules + HetMethanogenesisofBiomass.energyJoules + HetMethanogenesisofDOM.energyJoules;
Aut.energyJoules = AutAssimCO2.energyJoules + AutAssimNO3.energyJoules + AutAssimNH4.energyJoules;
Met.energyJoules = MetAssimCH4.energyJoules + MetAssimNO3.energyJoules + MetAssimNH4.energyJoules;
/***** Organisms must meet their respiratory energy demands (Exprsn. 16) *****/
0 = Het.respirationEnergyJoules + Het.energyJoules;
0 = Aut.respirationEnergyJoules + Aut.energyJoules;
0 = Met.respirationEnergyJoules + Met.energyJoules;

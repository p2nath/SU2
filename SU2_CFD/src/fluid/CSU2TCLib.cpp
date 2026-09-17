/*!
 * \file CSU2TCLib.cpp
 * \brief Source of user defined 2T nonequilibrium gas model.
 * \author C. Garbacz, W. Maier, S. R. Copeland, J. Needels
 * \version 8.1.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2024, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../../include/fluid/CSU2TCLib.hpp"
#include "../../../Common/include/option_structure.hpp"

CSU2TCLib::CSU2TCLib(const CConfig* config, unsigned short val_nDim, bool viscous): CNEMOGas(config, val_nDim){

  unsigned short maxEl = 0;
  su2double mf = 0.0;

  const auto MassFrac_Freestream = config->GetGas_Composition();

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    mf += MassFrac_Freestream[iSpecies];

  /*--- Allocate vectors for gas properties ---*/
  nElStates.resize(nSpecies,0);
  CharVibTemp.resize(nSpecies,0.0);
  RotationModes.resize(nSpecies,0.0);
  Diss.resize(nSpecies,0.0);
  A.resize(5,0.0);
  Omega11.resize(nSpecies,nSpecies,4,0.0);
  Omega22.resize(nSpecies,nSpecies,4,0.0);
  RxnConstantTable.resize(6,5) = su2double(0.0);
  CatRecombTable.resize(nSpecies,2) = 0;
  Blottner.resize(nSpecies,3)  = su2double(0.0);
  taus.resize(nSpecies,0.0);
  eve_eq.resize(nSpecies,0.0);
  eve.resize(nSpecies,0.0);

  if (viscous) {
    MolarFracWBE.resize(nSpecies,0.0);
    phis.resize(nSpecies,0.0);
    mus.resize(nSpecies,0.0);
  }

  if (gas_model =="ARGON"){
    if (nSpecies != 1) {
      SU2_MPI::Error("CONFIG ERROR: nSpecies mismatch between gas model & gas composition", CURRENT_FUNCTION);
    }
    mf = 0.0;
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      mf += MassFrac_Freestream[iSpecies];
    if (mf != 1.0) {
      SU2_MPI::Error("CONFIG ERROR: Intial gas mass fractions do not sum to 1!", CURRENT_FUNCTION);
    }

    /*--- Define parameters of the gas model ---*/
    gamma       = 1.667;
    nReactions  = 0;

    // Molar mass [kg/kmol]
    MolarMass[0] = 39.948;
    // Rotational modes of energy storage
    RotationModes[0] = 0.0;
    // Characteristic vibrational temperatures
    CharVibTemp[0] = 0.0;

    Enthalpy_Formation[0] = 0.0;
    Ref_Temperature[0] = 0.0;
    nElStates[0] = 7;

    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      maxEl = max(maxEl, nElStates[iSpecies]);

    /*--- Allocate and initialize electron data arrays ---*/
    CharElTemp.resize(nSpecies,maxEl) = su2double(0.0);
    ElDegeneracy.resize(nSpecies,maxEl) = su2double(0.0);

    /*--- AR: Blottner coefficients. ---*/
    Blottner(0,0) = 3.83444322E-03;   Blottner(0,1) = 6.74718764E-01;   Blottner(0,2) = -1.24290388E+01;

    /*--- AR: 7 states ---*/
    CharElTemp(0,0) = 0.000000000000000E+00;
    CharElTemp(0,1) = 1.611135736988230E+05;
    CharElTemp(0,2) = 1.625833076870950E+05;
    CharElTemp(0,3) = 1.636126382960720E+05;
    CharElTemp(0,4) = 1.642329518358000E+05;
    CharElTemp(0,5) = 1.649426852542080E+05;
    CharElTemp(0,6) = 1.653517702884570E+05;
    ElDegeneracy(0,0) = 1;
    ElDegeneracy(0,1) = 9;
    ElDegeneracy(0,2) = 21;
    ElDegeneracy(0,3) = 7;
    ElDegeneracy(0,4) = 3;
    ElDegeneracy(0,5) = 5;
    ElDegeneracy(0,6) = 15;

    /*--- Catalytic wall table ---*/
    // Creation/Destruction (+1/-1), Index of monoatomic reactants.
    // Argon not used.
    CatRecombTable(0,0) = 0; CatRecombTable(0,1) = 0;

    /*--- Values used in the Sutherland's formula. ---*/
    if (viscous) {
      //F.M. White, Viscous Fluid Flow, 3rd ed., McGraw-Hill, 2006.
      mu_ref[0] = 2.125E-5;
      k_ref[0] = 0.0163;
      Sm_ref[0] = 114.0;
      Sk_ref[0] = 170;
    }

  } else if (gas_model == "N2"){
    /*--- Check for errors in the initialization ---*/
    if (nSpecies != 2) {
      SU2_MPI::Error("CONFIG ERROR: nSpecies mismatch between gas model & gas composition", CURRENT_FUNCTION);
    }
    mf = 0.0;
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      mf += MassFrac_Freestream[iSpecies];
    if (mf != 1.0) {
      SU2_MPI::Error("CONFIG ERROR: Intial gas mass fractions do not sum to 1!", CURRENT_FUNCTION);
    }

    /*--- Define parameters of the gas model ---*/
    gamma       = 1.4;
    nReactions  = 2;

    Reactions.resize(nReactions,2,6,0.0);
    ArrheniusCoefficient.resize(nReactions,0.0);
    ArrheniusEta.resize(nReactions,0.0);
    ArrheniusTheta.resize(nReactions,0.0);
    Tcf_a.resize(nReactions,0.0);
    Tcf_b.resize(nReactions,0.0);
    Tcb_a.resize(nReactions,0.0);
    Tcb_b.resize(nReactions,0.0);

    /*--- Assign gas properties ---*/
    // Rotational modes of energy storage
    RotationModes[0] = 2.0;
    RotationModes[1] = 0.0;
    // Molar mass [kg/kmol]
    MolarMass[0] = 2.0*14.0067;
    MolarMass[1] = 14.0067;
    // Characteristic vibrational temperatures
    CharVibTemp[0] = 3395.0;
    CharVibTemp[1] = 0.0;
    // Formation enthalpy: (JANAF values [KJ/Kmol])
    // J/kg - from Scalabrin
    Enthalpy_Formation[0] = 0.0;      //N2
    Enthalpy_Formation[1] = 3.36E7;   //N
    // Reference temperature (JANAF values, [K])
    Ref_Temperature[0] = 0.0;
    Ref_Temperature[1] = 0.0;
    // Blottner viscosity coefficients
    // A                       // B                       // C
    Blottner(0,0) = 2.68E-2;   Blottner(0,1) = 3.18E-1;   Blottner(0,2) = -1.13E1;  // N2
    Blottner(1,0) = 1.16E-2;   Blottner(1,1) = 6.03E-1;   Blottner(1,2) = -1.24E1;  // N
    // Number of electron states
    nElStates[0] = 15;                    // N2
    nElStates[1] = 3;                     // N
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      maxEl = max(maxEl, nElStates[iSpecies]);

    /*--- Allocate and initialize electron data arrays ---*/
    CharElTemp.resize(nSpecies,maxEl) = su2double(0.0);
    ElDegeneracy.resize(nSpecies,maxEl) = su2double(0.0);

    /*--- Assign values to data structures ---*/
    // N2: 15 states
    CharElTemp(0,0)  = 0.000000000000000E+00;
    CharElTemp(0,1)  = 7.223156514095200E+04;
    CharElTemp(0,2)  = 8.577862640384000E+04;
    CharElTemp(0,3)  = 8.605026716160000E+04;
    CharElTemp(0,4)  = 9.535118627874400E+04;
    CharElTemp(0,5)  = 9.805635702203200E+04;
    CharElTemp(0,6)  = 9.968267656935200E+04;
    CharElTemp(0,7)  = 1.048976467715200E+05;
    CharElTemp(0,8)  = 1.116489555200000E+05;
    CharElTemp(0,9)  = 1.225836470400000E+05;
    CharElTemp(0,10) = 1.248856873600000E+05;
    CharElTemp(0,11) = 1.282476158188320E+05;
    CharElTemp(0,12) = 1.338060936000000E+05;
    CharElTemp(0,13) = 1.404296391107200E+05;
    CharElTemp(0,14) = 1.504958859200000E+05;
    ElDegeneracy(0,0)  = 1;
    ElDegeneracy(0,1)  = 3;
    ElDegeneracy(0,2)  = 6;
    ElDegeneracy(0,3)  = 6;
    ElDegeneracy(0,4)  = 3;
    ElDegeneracy(0,5)  = 1;
    ElDegeneracy(0,6)  = 2;
    ElDegeneracy(0,7)  = 2;
    ElDegeneracy(0,8)  = 5;
    ElDegeneracy(0,9)  = 1;
    ElDegeneracy(0,10) = 6;
    ElDegeneracy(0,11) = 6;
    ElDegeneracy(0,12) = 10;
    ElDegeneracy(0,13) = 6;
    ElDegeneracy(0,14) = 6;
    // N: 3 states
    CharElTemp(1,0) = 0.000000000000000E+00;
    CharElTemp(1,1) = 2.766469645581980E+04;
    CharElTemp(1,2) = 4.149309313560210E+04;
    ElDegeneracy(1,0) = 4;
    ElDegeneracy(1,1) = 10;
    ElDegeneracy(1,2) = 6;
    /*--- Set Arrhenius coefficients for chemical reactions ---*/
    // Note: Data lists coefficients in (cm^3/mol-s) units, need to convert
    //       to (m^3/kmol-s) to be consistent with the rest of the code
    // Pre-exponential factor
    ArrheniusCoefficient[0]  = 7.0E21;
    ArrheniusCoefficient[1]  = 3.0E22;
    // Rate-controlling temperature exponent
    ArrheniusEta[0]  = -1.60;
    ArrheniusEta[1]  = -1.60;
    // Characteristic temperature
    ArrheniusTheta[0] = 113200.0;
    ArrheniusTheta[1] = 113200.0;
    /*--- Set reaction maps ---*/
    // N2 + N2 -> 2N + N2
    Reactions(0,0,0)=0;   Reactions(0,0,1)=0;   Reactions(0,0,2)=nSpecies;
    Reactions(0,1,0)=1;   Reactions(0,1,1)=1;   Reactions(0,1,2) =0;
    // N2 + N -> 2N + N
    Reactions(1,0,0)=0;   Reactions(1,0,1)=1;   Reactions(1,0,2)=nSpecies;
    Reactions(1,1,0)=1;   Reactions(1,1,1)=1;   Reactions(1,1,2)=1;
    /*--- Set rate-controlling temperature exponents ---*/
    //  -----------  Tc = Ttr^a * Tve^b  -----------
    //
    // Forward Reactions
    //   Dissociation:      a = 0.5, b = 0.5  (OR a = 0.7, b =0.3)
    //   Exchange:          a = 1,   b = 0
    //   Impact ionization: a = 0,   b = 1
    //
    // Backward Reactions
    //   Recomb ionization:      a = 0, b = 1
    //   Impact ionization:      a = 0, b = 1
    //   N2 impact dissociation: a = 0, b = 1
    //   Others:                 a = 1, b = 0
    Tcf_a[0] = 0.5; Tcf_b[0] = 0.5; Tcb_a[0] = 1;  Tcb_b[0] = 0;
    Tcf_a[1] = 0.5; Tcf_b[1] = 0.5; Tcb_a[1] = 1;  Tcb_b[1] = 0;

    /*--- Dissociation potential [KJ/kg] ---*/
    Diss[0] = 3.36E4;
    Diss[1] = 0.0;

    /*--- Collision integral data ---*/
    // Index 1: collider
    // Index 2: partner
    // Index 3: A1, A2, A3
    Omega11(0,0,0) = -6.0614558E-03;  Omega11(0,0,1) = 1.2689102E-01;   Omega11(0,0,2) = -1.0616948E+00;  Omega11(0,0,3) = 8.0955466E+02;
    Omega11(0,1,0) = -1.0796249E-02;  Omega11(0,1,1) = 2.2656509E-01;   Omega11(0,1,2) = -1.7910602E+00;  Omega11(0,1,3) = 4.0455218E+03;
    Omega11(1,0,0) = -1.0796249E-02;  Omega11(1,0,1) = 2.2656509E-01;   Omega11(1,0,2) = -1.7910602E+00;  Omega11(1,0,3) = 4.0455218E+03;
    Omega11(1,1,0) = -9.6083779E-03;  Omega11(1,1,1) = 2.0938971E-01;   Omega11(1,1,2) = -1.7386904E+00;  Omega11(1,1,3) = 3.3587983E+03;
    Omega22(0,0,0) = -7.6303990E-03;  Omega22(0,0,1) = 1.6878089E-01;   Omega22(0,0,2) = -1.4004234E+00;  Omega22(0,0,3) = 2.1427708E+03;
    Omega22(0,1,0) = -8.3493693E-03;  Omega22(0,1,1) = 1.7808911E-01;   Omega22(0,1,2) = -1.4466155E+00;  Omega22(0,1,3) = 1.9324210E+03;
    Omega22(1,0,0) = -8.3493693E-03;  Omega22(1,0,1) = 1.7808911E-01;   Omega22(1,0,2) = -1.4466155E+00;  Omega22(1,0,3) = 1.9324210E+03;
    Omega22(1,1,0) = -7.7439615E-03;  Omega22(1,1,1) = 1.7129007E-01;   Omega22(1,1,2) = -1.4809088E+00;  Omega22(1,1,3) = 2.1284951E+03;

    /*--- Catalytic wall table ---*/
    // Creation/Destruction (+1/-1), Index of monoatomic reactants.
    // Monoatomic species (N,O) recombine into diaatomic (N2, O2)
    CatRecombTable(0,0) =  1; CatRecombTable(0,1) = 1;
    CatRecombTable(1,0) = -1; CatRecombTable(1,1) = 1;

    /*--- Values used in the Sutherland's formula. ---*/
    if (viscous) {
      //F.M. White, Viscous Fluid Flow, 3rd ed., McGraw-Hill, 2006.
      k_ref[0] = 0.0242;
      mu_ref[0] = 1.663E-5;
      Sm_ref[0] = 107.0;
      Sk_ref[0] = 150.0;
    }

  } else if (gas_model == "AIR-5"){

    /*--- Check for errors in the initialization ---*/
    if (nSpecies != 5) {
      SU2_MPI::Error("CONFIG ERROR: nSpecies mismatch between gas model & gas composition",CURRENT_FUNCTION);
    }
    mf = 0.0;
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      mf += MassFrac_Freestream[iSpecies];
    if (mf != 1.0) {
      SU2_MPI::Error("CONFIG ERROR: Intial gas mass fractions do not sum to 1!", CURRENT_FUNCTION);
    }

    /*--- Define parameters of the gas model ---*/
    gamma       = 1.4;
    nReactions  = 17;

    Reactions.resize(nReactions,2,6,0.0);
    ArrheniusCoefficient.resize(nReactions,0.0);
    ArrheniusEta.resize(nReactions,0.0);
    ArrheniusTheta.resize(nReactions,0.0);
    Tcf_a.resize(nReactions,0.0);
    Tcf_b.resize(nReactions,0.0);
    Tcb_a.resize(nReactions,0.0);
    Tcb_b.resize(nReactions,0.0);

    /*--- Assign gas properties ---*/
    // Rotational modes of energy storage
    RotationModes[0] = 2.0;
    RotationModes[1] = 2.0;
    RotationModes[2] = 2.0;
    RotationModes[3] = 0.0;
    RotationModes[4] = 0.0;
    // Molar mass [kg/kmol]
    MolarMass[0] = 2.0*14.0067;
    MolarMass[1] = 2.0*15.9994;
    MolarMass[2] = 14.0067+15.9994;
    MolarMass[3] = 14.0067;
    MolarMass[4] = 15.9994;
    //Characteristic vibrational temperatures
    CharVibTemp[0] = 3395.0;
    CharVibTemp[1] = 2239.0;
    CharVibTemp[2] = 2817.0;
    CharVibTemp[3] = 0.0;
    CharVibTemp[4] = 0.0;
    // Formation enthalpy: (Scalabrin values, J/kg)
    Enthalpy_Formation[0] = 0.0;      //N2
    Enthalpy_Formation[1] = 0.0;      //O2
    Enthalpy_Formation[2] = 3.0E6;    //NO
    Enthalpy_Formation[3] = 3.36E7;   //N
    Enthalpy_Formation[4] = 1.54E7;   //O
    // Reference temperature (JANAF values, [K])
    Ref_Temperature[0] = 0.0;
    Ref_Temperature[1] = 0.0;
    Ref_Temperature[2] = 0.0;
    Ref_Temperature[3] = 0.0;
    Ref_Temperature[4] = 0.0;
    // Blottner viscosity coefficients
    // A                        // B                        // C
    Blottner(0,0) = 2.68E-2;   Blottner(0,1) =  3.18E-1;  Blottner(0,2) = -1.13E1;  // N2
    Blottner(1,0) = 4.49E-2;   Blottner(1,1) = -8.26E-2;  Blottner(1,2) = -9.20E0;  // O2
    Blottner(2,0) = 4.36E-2;   Blottner(2,1) = -3.36E-2;  Blottner(2,2) = -9.58E0;  // NO
    Blottner(3,0) = 1.16E-2;   Blottner(3,1) =  6.03E-1;  Blottner(3,2) = -1.24E1;  // N
    Blottner(4,0) = 2.03E-2;   Blottner(4,1) =  4.29E-1;  Blottner(4,2) = -1.16E1;  // O
    // Number of electron states
    nElStates[0] = 15;                    // N2
    nElStates[1] = 7;                     // O2
    nElStates[2] = 16;                    // NO
    nElStates[3] = 3;                     // N
    nElStates[4] = 5;                     // O
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      maxEl = max(maxEl, nElStates[iSpecies]);
    /*--- Allocate and initialize electron data arrays ---*/
    CharElTemp.resize(nSpecies,maxEl) = su2double(0.0);
    ElDegeneracy.resize(nSpecies,maxEl) = su2double(0.0);

    //N2: 15 states
    CharElTemp(0,0)  = 0.000000000000000E+00;
    CharElTemp(0,1)  = 7.223156514095200E+04;
    CharElTemp(0,2)  = 8.577862640384000E+04;
    CharElTemp(0,3)  = 8.605026716160000E+04;
    CharElTemp(0,4)  = 9.535118627874400E+04;
    CharElTemp(0,5)  = 9.805635702203200E+04;
    CharElTemp(0,6)  = 9.968267656935200E+04;
    CharElTemp(0,7)  = 1.048976467715200E+05;
    CharElTemp(0,8)  = 1.116489555200000E+05;
    CharElTemp(0,9)  = 1.225836470400000E+05;
    CharElTemp(0,10) = 1.248856873600000E+05;
    CharElTemp(0,11) = 1.282476158188320E+05;
    CharElTemp(0,12) = 1.338060936000000E+05;
    CharElTemp(0,13) = 1.404296391107200E+05;
    CharElTemp(0,14) = 1.504958859200000E+05;
    ElDegeneracy(0,0)  = 1;
    ElDegeneracy(0,1)  = 3;
    ElDegeneracy(0,2)  = 6;
    ElDegeneracy(0,3)  = 6;
    ElDegeneracy(0,4)  = 3;
    ElDegeneracy(0,5)  = 1;
    ElDegeneracy(0,6)  = 2;
    ElDegeneracy(0,7)  = 2;
    ElDegeneracy(0,8)  = 5;
    ElDegeneracy(0,9)  = 1;
    ElDegeneracy(0,10) = 6;
    ElDegeneracy(0,11) = 6;
    ElDegeneracy(0,12) = 10;
    ElDegeneracy(0,13) = 6;
    ElDegeneracy(0,14) = 6;
    // O2: 7 states
    CharElTemp(1,0) = 0.000000000000000E+00;
    CharElTemp(1,1) = 1.139156019700800E+04;
    CharElTemp(1,2) = 1.898473947826400E+04;
    CharElTemp(1,3) = 4.755973576639200E+04;
    CharElTemp(1,4) = 4.991242097343200E+04;
    CharElTemp(1,5) = 5.092268575561600E+04;
    CharElTemp(1,6) = 7.189863255967200E+04;
    ElDegeneracy(1,0) = 3;
    ElDegeneracy(1,1) = 2;
    ElDegeneracy(1,2) = 1;
    ElDegeneracy(1,3) = 1;
    ElDegeneracy(1,4) = 6;
    ElDegeneracy(1,5) = 3;
    ElDegeneracy(1,6) = 3;
    // NO: 16 states
    CharElTemp(2,0)  = 0.000000000000000E+00;
    CharElTemp(2,1)  = 5.467345760000000E+04;
    CharElTemp(2,2)  = 6.317139627802400E+04;
    CharElTemp(2,3)  = 6.599450342445600E+04;
    CharElTemp(2,4)  = 6.906120960000000E+04;
    CharElTemp(2,5)  = 7.049998480000000E+04;
    CharElTemp(2,6)  = 7.491055017560000E+04;
    CharElTemp(2,7)  = 7.628875293968000E+04;
    CharElTemp(2,8)  = 8.676188537552000E+04;
    CharElTemp(2,9)  = 8.714431182368000E+04;
    CharElTemp(2,10) = 8.886077063728000E+04;
    CharElTemp(2,11) = 8.981755614528000E+04;
    CharElTemp(2,12) = 8.988445919208000E+04;
    CharElTemp(2,13) = 9.042702132000000E+04;
    CharElTemp(2,14) = 9.064283760000000E+04;
    CharElTemp(2,15) = 9.111763341600000E+04;
    ElDegeneracy(2,0)  = 4;
    ElDegeneracy(2,1)  = 8;
    ElDegeneracy(2,2)  = 2;
    ElDegeneracy(2,3)  = 4;
    ElDegeneracy(2,4)  = 4;
    ElDegeneracy(2,5)  = 4;
    ElDegeneracy(2,6)  = 4;
    ElDegeneracy(2,7)  = 2;
    ElDegeneracy(2,8)  = 4;
    ElDegeneracy(2,9)  = 2;
    ElDegeneracy(2,10) = 4;
    ElDegeneracy(2,11) = 4;
    ElDegeneracy(2,12) = 2;
    ElDegeneracy(2,13) = 2;
    ElDegeneracy(2,14) = 2;
    ElDegeneracy(2,15) = 4;
    // N: 3 states
    CharElTemp(3,0) = 0.000000000000000E+00;
    CharElTemp(3,1) = 2.766469645581980E+04;
    CharElTemp(3,2) = 4.149309313560210E+04;
    ElDegeneracy(3,0)= 4;
    ElDegeneracy(3,1)= 10;
    ElDegeneracy(3,2)= 6;
    // O: 5 states
    CharElTemp(4,0) = 0.000000000000000E+00;
    CharElTemp(4,1) = 2.277077570280000E+02;
    CharElTemp(4,2) = 3.265688785704000E+02;
    CharElTemp(4,3) = 2.283028632262240E+04;
    CharElTemp(4,4) = 4.861993036434160E+04;
    ElDegeneracy(4,0) = 5;
    ElDegeneracy(4,1) = 3;
    ElDegeneracy(4,2) = 1;
    ElDegeneracy(4,3) = 5;
    ElDegeneracy(4,4) = 1;
    /*--- Set reaction maps ---*/
    // N2 dissociation
    Reactions(0,0,0)=0;    Reactions(0,0,1)=0;   Reactions(0,0,2)=nSpecies;    Reactions(0,1,0)=3;   Reactions(0,1,1)=3;   Reactions(0,1,2) =0;
    Reactions(1,0,0)=0;    Reactions(1,0,1)=1;   Reactions(1,0,2)=nSpecies;    Reactions(1,1,0)=3;   Reactions(1,1,1)=3;   Reactions(1,1,2) =1;
    Reactions(2,0,0)=0;    Reactions(2,0,1)=2;   Reactions(2,0,2)=nSpecies;    Reactions(2,1,0)=3;   Reactions(2,1,1)=3;   Reactions(2,1,2) =2;
    Reactions(3,0,0)=0;    Reactions(3,0,1)=3;   Reactions(3,0,2)=nSpecies;    Reactions(3,1,0)=3;   Reactions(3,1,1)=3;   Reactions(3,1,2) =3;
    Reactions(4,0,0)=0;    Reactions(4,0,1)=4;   Reactions(4,0,2)=nSpecies;    Reactions(4,1,0)=3;   Reactions(4,1,1)=3;   Reactions(4,1,2) =4;
    // O2 dissociation
    Reactions(5,0,0)=1;    Reactions(5,0,1)=0;   Reactions(5,0,2)=nSpecies;    Reactions(5,1,0)=4;   Reactions(5,1,1)=4;   Reactions(5,1,2) =0;
    Reactions(6,0,0)=1;    Reactions(6,0,1)=1;   Reactions(6,0,2)=nSpecies;    Reactions(6,1,0)=4;   Reactions(6,1,1)=4;   Reactions(6,1,2) =1;
    Reactions(7,0,0)=1;    Reactions(7,0,1)=2;   Reactions(7,0,2)=nSpecies;    Reactions(7,1,0)=4;   Reactions(7,1,1)=4;   Reactions(7,1,2) =2;
    Reactions(8,0,0)=1;    Reactions(8,0,1)=3;   Reactions(8,0,2)=nSpecies;    Reactions(8,1,0)=4;   Reactions(8,1,1)=4;   Reactions(8,1,2) =3;
    Reactions(9,0,0)=1;    Reactions(9,0,1)=4;   Reactions(9,0,2)=nSpecies;    Reactions(9,1,0)=4;   Reactions(9,1,1)=4;   Reactions(9,1,2) =4;
    // NO dissociation
    Reactions(10,0,0)=2;   Reactions(10,0,1)=0;  Reactions(10,0,2)=nSpecies;   Reactions(10,1,0)=3;  Reactions(10,1,1)=4;    Reactions(10,1,2) =0;
    Reactions(11,0,0)=2;   Reactions(11,0,1)=1;  Reactions(11,0,2)=nSpecies;   Reactions(11,1,0)=3;  Reactions(11,1,1)=4;    Reactions(11,1,2) =1;
    Reactions(12,0,0)=2;   Reactions(12,0,1)=2;  Reactions(12,0,2)=nSpecies;   Reactions(12,1,0)=3;  Reactions(12,1,1)=4;    Reactions(12,1,2) =2;
    Reactions(13,0,0)=2;   Reactions(13,0,1)=3;  Reactions(13,0,2)=nSpecies;   Reactions(13,1,0)=3;  Reactions(13,1,1)=4;    Reactions(13,1,2) =3;
    Reactions(14,0,0)=2;   Reactions(14,0,1)=4;  Reactions(14,0,2)=nSpecies;   Reactions(14,1,0)=3;  Reactions(14,1,1)=4;    Reactions(14,1,2) =4;
    // N2 + O -> NO + N
    Reactions(15,0,0)=0;   Reactions(15,0,1)=4;  Reactions(15,0,2)=nSpecies;   Reactions(15,1,0)=2;  Reactions(15,1,1)=3;    Reactions(15,1,2)= nSpecies;
    // NO + O -> O2 + N
    Reactions(16,0,0)=2;   Reactions(16,0,1)=4;  Reactions(16,0,2)=nSpecies;   Reactions(16,1,0)=1;  Reactions(16,1,1)=3;    Reactions(16,1,2)= nSpecies;
    /*--- Set Arrhenius coefficients for reactions ---*/
    // Pre-exponential factor
    ArrheniusCoefficient[0]  = 7.0E21;
    ArrheniusCoefficient[1]  = 7.0E21;
    ArrheniusCoefficient[2]  = 7.0E21;
    ArrheniusCoefficient[3]  = 3.0E22;
    ArrheniusCoefficient[4]  = 3.0E22;
    ArrheniusCoefficient[5]  = 2.0E21;
    ArrheniusCoefficient[6]  = 2.0E21;
    ArrheniusCoefficient[7]  = 2.0E21;
    ArrheniusCoefficient[8]  = 1.0E22;
    ArrheniusCoefficient[9]  = 1.0E22;
    ArrheniusCoefficient[10] = 5.0E15;
    ArrheniusCoefficient[11] = 5.0E15;
    ArrheniusCoefficient[12] = 5.0E15;
    ArrheniusCoefficient[13] = 1.1E17;
    ArrheniusCoefficient[14] = 1.1E17;
    ArrheniusCoefficient[15] = 6.4E17;
    ArrheniusCoefficient[16] = 8.4E12;
    // Rate-controlling temperature exponent
    ArrheniusEta[0]  = -1.60;
    ArrheniusEta[1]  = -1.60;
    ArrheniusEta[2]  = -1.60;
    ArrheniusEta[3]  = -1.60;
    ArrheniusEta[4]  = -1.60;
    ArrheniusEta[5]  = -1.50;
    ArrheniusEta[6]  = -1.50;
    ArrheniusEta[7]  = -1.50;
    ArrheniusEta[8]  = -1.50;
    ArrheniusEta[9]  = -1.50;
    ArrheniusEta[10] = 0.0;
    ArrheniusEta[11] = 0.0;
    ArrheniusEta[12] = 0.0;
    ArrheniusEta[13] = 0.0;
    ArrheniusEta[14] = 0.0;
    ArrheniusEta[15] = -1.0;
    ArrheniusEta[16] = 0.0;
    // Characteristic temperature
    ArrheniusTheta[0]  = 113200.0;
    ArrheniusTheta[1]  = 113200.0;
    ArrheniusTheta[2]  = 113200.0;
    ArrheniusTheta[3]  = 113200.0;
    ArrheniusTheta[4]  = 113200.0;
    ArrheniusTheta[5]  = 59500.0;
    ArrheniusTheta[6]  = 59500.0;
    ArrheniusTheta[7]  = 59500.0;
    ArrheniusTheta[8]  = 59500.0;
    ArrheniusTheta[9]  = 59500.0;
    ArrheniusTheta[10] = 75500.0;
    ArrheniusTheta[11] = 75500.0;
    ArrheniusTheta[12] = 75500.0;
    ArrheniusTheta[13] = 75500.0;
    ArrheniusTheta[14] = 75500.0;
    ArrheniusTheta[15] = 38400.0;
    ArrheniusTheta[16] = 19450.0;
    /*--- Set rate-controlling temperature exponents ---*/
    //  -----------  Tc = Ttr^a * Tve^b  -----------
    //
    // Forward Reactions
    //   Dissociation:      a = 0.5, b = 0.5  (OR a = 0.7, b =0.3)
    //   Exchange:          a = 1,   b = 0
    //   Impact ionization: a = 0,   b = 1
    //
    // Backward Reactions
    //   Recomb ionization:      a = 0, b = 1
    //   Impact ionization:      a = 0, b = 1
    //   N2 impact dissociation: a = 0, b = 1
    //   Others:                 a = 1, b = 0
    Tcf_a[0]  = 0.5; Tcf_b[0]  = 0.5; Tcb_a[0]  = 1;  Tcb_b[0] = 0;
    Tcf_a[1]  = 0.5; Tcf_b[1]  = 0.5; Tcb_a[1]  = 1;  Tcb_b[1] = 0;
    Tcf_a[2]  = 0.5; Tcf_b[2]  = 0.5; Tcb_a[2]  = 1;  Tcb_b[2] = 0;
    Tcf_a[3]  = 0.5; Tcf_b[3]  = 0.5; Tcb_a[3]  = 1;  Tcb_b[3] = 0;
    Tcf_a[4]  = 0.5; Tcf_b[4]  = 0.5; Tcb_a[4]  = 1;  Tcb_b[4] = 0;
    Tcf_a[5]  = 0.5; Tcf_b[5]  = 0.5; Tcb_a[5]  = 1;  Tcb_b[5] = 0;
    Tcf_a[6]  = 0.5; Tcf_b[6]  = 0.5; Tcb_a[6]  = 1;  Tcb_b[6] = 0;
    Tcf_a[7]  = 0.5; Tcf_b[7]  = 0.5; Tcb_a[7]  = 1;  Tcb_b[7] = 0;
    Tcf_a[8]  = 0.5; Tcf_b[8]  = 0.5; Tcb_a[8]  = 1;  Tcb_b[8] = 0;
    Tcf_a[9]  = 0.5; Tcf_b[9]  = 0.5; Tcb_a[9]  = 1;  Tcb_b[9] = 0;
    Tcf_a[10] = 0.5; Tcf_b[10] = 0.5; Tcb_a[10] = 1;  Tcb_b[10] = 0;
    Tcf_a[11] = 0.5; Tcf_b[11] = 0.5; Tcb_a[11] = 1;  Tcb_b[11] = 0;
    Tcf_a[12] = 0.5; Tcf_b[12] = 0.5; Tcb_a[12] = 1;  Tcb_b[12] = 0;
    Tcf_a[13] = 0.5; Tcf_b[13] = 0.5; Tcb_a[13] = 1;  Tcb_b[13] = 0;
    Tcf_a[14] = 0.5; Tcf_b[14] = 0.5; Tcb_a[14] = 1;  Tcb_b[14] = 0;
    Tcf_a[15] = 1.0; Tcf_b[15] = 0.0; Tcb_a[15] = 1;  Tcb_b[15] = 0;
    Tcf_a[16] = 1.0; Tcf_b[16] = 0.0; Tcb_a[16] = 1;  Tcb_b[16] = 0;
    /*--- Collision integral data ---*/
    // Index 1: collider
    // Index 2: partner
    // Index 3: A1, A2, A3
    // Omega^(1,1) ----------------------
    //N2
    Omega11(0,0,0) = -6.0614558E-03;  Omega11(0,0,1) = 1.2689102E-01;   Omega11(0,0,2) = -1.0616948E+00;  Omega11(0,0,3) = 8.0955466E+02;
    Omega11(0,1,0) = -3.7959091E-03;  Omega11(0,1,1) = 9.5708295E-02;   Omega11(0,1,2) = -1.0070611E+00;  Omega11(0,1,3) = 8.9392313E+02;
    Omega11(0,2,0) = -1.9295666E-03;  Omega11(0,2,1) = 2.7995735E-02;   Omega11(0,2,2) = -3.1588514E-01;  Omega11(0,2,3) = 1.2880734E+02;
    Omega11(0,3,0) = -1.0796249E-02;  Omega11(0,3,1) = 2.2656509E-01;   Omega11(0,3,2) = -1.7910602E+00;  Omega11(0,3,3) = 4.0455218E+03;
    Omega11(0,4,0) = -2.7244269E-03;  Omega11(0,4,1) = 6.9587171E-02;   Omega11(0,4,2) = -7.9538667E-01;  Omega11(0,4,3) = 4.0673730E+02;
    //O2
    Omega11(1,0,0) = -3.7959091E-03;  Omega11(1,0,1) = 9.5708295E-02;   Omega11(1,0,2) = -1.0070611E+00;  Omega11(1,0,3) = 8.9392313E+02;
    Omega11(1,1,0) = -8.0682650E-04;  Omega11(1,1,1) = 1.6602480E-02;   Omega11(1,1,2) = -3.1472774E-01;  Omega11(1,1,3) = 1.4116458E+02;
    Omega11(1,2,0) = -6.4433840E-04;  Omega11(1,2,1) = 8.5378580E-03;   Omega11(1,2,2) = -2.3225102E-01;  Omega11(1,2,3) = 1.1371608E+02;
    Omega11(1,3,0) = -1.1453028E-03;  Omega11(1,3,1) = 1.2654140E-02;   Omega11(1,3,2) = -2.2435218E-01;  Omega11(1,3,3) = 7.7201588E+01;
    Omega11(1,4,0) = -4.8405803E-03;  Omega11(1,4,1) = 1.0297688E-01;   Omega11(1,4,2) = -9.6876576E-01;  Omega11(1,4,3) = 6.1629812E+02;
    //NO
    Omega11(2,0,0) = -1.9295666E-03;  Omega11(2,0,1) = 2.7995735E-02;   Omega11(2,0,2) = -3.1588514E-01;  Omega11(2,0,3) = 1.2880734E+02;
    Omega11(2,1,0) = -6.4433840E-04;  Omega11(2,1,1) = 8.5378580E-03;   Omega11(2,1,2) = -2.3225102E-01;  Omega11(2,1,3) = 1.1371608E+02;
    Omega11(2,2,0) = -0.0000000E+00;  Omega11(2,2,1) = -1.1056066E-02;  Omega11(2,2,2) = -5.9216250E-02;  Omega11(2,2,3) = 7.2542367E+01;
    Omega11(2,3,0) = -1.5770918E-03;  Omega11(2,3,1) = 1.9578381E-02;   Omega11(2,3,2) = -2.7873624E-01;  Omega11(2,3,3) = 9.9547944E+01;
    Omega11(2,4,0) = -1.0885815E-03;  Omega11(2,4,1) = 1.1883688E-02;   Omega11(2,4,2) = -2.1844909E-01;  Omega11(2,4,3) = 7.5512560E+01;
    //N
    Omega11(3,0,0) = -1.0796249E-02;  Omega11(3,0,1) = 2.2656509E-01;   Omega11(3,0,2) = -1.7910602E+00;  Omega11(3,0,3) = 4.0455218E+03;
    Omega11(3,1,0) = -1.1453028E-03;  Omega11(3,1,1) = 1.2654140E-02;   Omega11(3,1,2) = -2.2435218E-01;  Omega11(3,1,3) = 7.7201588E+01;
    Omega11(3,2,0) = -1.5770918E-03;  Omega11(3,2,1) = 1.9578381E-02;   Omega11(3,2,2) = -2.7873624E-01;  Omega11(3,2,3) = 9.9547944E+01;
    Omega11(3,3,0) = -9.6083779E-03;  Omega11(3,3,1) = 2.0938971E-01;   Omega11(3,3,2) = -1.7386904E+00;  Omega11(3,3,3) = 3.3587983E+03;
    Omega11(3,4,0) = -7.8147689E-03;  Omega11(3,4,1) = 1.6792705E-01;   Omega11(3,4,2) = -1.4308628E+00;  Omega11(3,4,3) = 1.6628859E+03;
    //O
    Omega11(4,0,0) = -2.7244269E-03;  Omega11(4,0,1) = 6.9587171E-02;   Omega11(4,0,2) = -7.9538667E-01;  Omega11(4,0,3) = 4.0673730E+02;
    Omega11(4,1,0) = -4.8405803E-03;  Omega11(4,1,1) = 1.0297688E-01;   Omega11(4,1,2) = -9.6876576E-01;  Omega11(4,1,3) = 6.1629812E+02;
    Omega11(4,2,0) = -1.0885815E-03;  Omega11(4,2,1) = 1.1883688E-02;   Omega11(4,2,2) = -2.1844909E-01;  Omega11(4,2,3) = 7.5512560E+01;
    Omega11(4,3,0) = -7.8147689E-03;  Omega11(4,3,1) = 1.6792705E-01;   Omega11(4,3,2) = -1.4308628E+00;  Omega11(4,3,3) = 1.6628859E+03;
    Omega11(4,4,0) = -6.4040535E-03;  Omega11(4,4,1) = 1.4629949E-01;   Omega11(4,4,2) = -1.3892121E+00;  Omega11(4,4,3) = 2.0903441E+03;

    // Omega^(2,2) ----------------------
    //N2
    Omega22(0,0,0) = -7.6303990E-03;  Omega22(0,0,1) = 1.6878089E-01;   Omega22(0,0,2) = -1.4004234E+00;  Omega22(0,0,3) = 2.1427708E+03;
    Omega22(0,1,0) = -8.0457321E-03;  Omega22(0,1,1) = 1.9228905E-01;   Omega22(0,1,2) = -1.7102854E+00;  Omega22(0,1,3) = 5.2213857E+03;
    Omega22(0,2,0) = -6.8237776E-03;  Omega22(0,2,1) = 1.4360616E-01;   Omega22(0,2,2) = -1.1922240E+00;  Omega22(0,2,3) = 1.2433086E+03;
    Omega22(0,3,0) = -8.3493693E-03;  Omega22(0,3,1) = 1.7808911E-01;   Omega22(0,3,2) = -1.4466155E+00;  Omega22(0,3,3) = 1.9324210E+03;
    Omega22(0,4,0) = -8.3110691E-03;  Omega22(0,4,1) = 1.9617877E-01;   Omega22(0,4,2) = -1.7205427E+00;  Omega22(0,4,3) = 4.0812829E+03;
    //O2
    Omega22(1,0,0) = -8.0457321E-03;  Omega22(1,0,1) = 1.9228905E-01;   Omega22(1,0,2) = -1.7102854E+00;  Omega22(1,0,3) = 5.2213857E+03;
    Omega22(1,1,0) = -6.2931612E-03;  Omega22(1,1,1) = 1.4624645E-01;   Omega22(1,1,2) = -1.3006927E+00;  Omega22(1,1,3) = 1.8066892E+03;
    Omega22(1,2,0) = -6.8508672E-03;  Omega22(1,2,1) = 1.5524564E-01;   Omega22(1,2,2) = -1.3479583E+00;  Omega22(1,2,3) = 2.0037890E+03;
    Omega22(1,3,0) = -1.0608832E-03;  Omega22(1,3,1) = 1.1782595E-02;   Omega22(1,3,2) = -2.1246301E-01;  Omega22(1,3,3) = 8.4561598E+01;
    Omega22(1,4,0) = -3.7969686E-03;  Omega22(1,4,1) = 7.6789981E-02;   Omega22(1,4,2) = -7.3056809E-01;  Omega22(1,4,3) = 3.3958171E+02;
    //NO
    Omega22(2,0,0) = -6.8237776E-03;  Omega22(2,0,1) = 1.4360616E-01;   Omega22(2,0,2) = -1.1922240E+00;  Omega22(2,0,3) = 1.2433086E+03;
    Omega22(2,1,0) = -6.8508672E-03;  Omega22(2,1,1) = 1.5524564E-01;   Omega22(2,1,2) = -1.3479583E+00;  Omega22(2,1,3) = 2.0037890E+03;
    Omega22(2,2,0) = -7.4942466E-03;  Omega22(2,2,1) = 1.6626193E-01;   Omega22(2,2,2) = -1.4107027E+00;  Omega22(2,2,3) = 2.3097604E+03;
    Omega22(2,3,0) = -1.4719259E-03;  Omega22(2,3,1) = 1.8446968E-02;   Omega22(2,3,2) = -2.6460411E-01;  Omega22(2,3,3) = 1.0911124E+02;
    Omega22(2,4,0) = -1.0066279E-03;  Omega22(2,4,1) = 1.1029264E-02;   Omega22(2,4,2) = -2.0671266E-01;  Omega22(2,4,3) = 8.2644384E+01;
    //N
    Omega22(3,0,0) = -8.3493693E-03;  Omega22(3,0,1) = 1.7808911E-01;   Omega22(3,0,2) = -1.4466155E+00;  Omega22(3,0,3) = 1.9324210E+03;
    Omega22(3,1,0) = -1.0608832E-03;  Omega22(3,1,1) = 1.1782595E-02;   Omega22(3,1,2) = -2.1246301E-01;  Omega22(3,1,3) = 8.4561598E+01;
    Omega22(3,2,0) = -1.4719259E-03;  Omega22(3,2,1) = 1.8446968E-02;   Omega22(3,2,2) = -2.6460411E-01;  Omega22(3,2,3) = 1.0911124E+02;
    Omega22(3,3,0) = -7.7439615E-03;  Omega22(3,3,1) = 1.7129007E-01;   Omega22(3,3,2) = -1.4809088E+00;  Omega22(3,3,3) = 2.1284951E+03;
    Omega22(3,4,0) = -5.0478143E-03;  Omega22(3,4,1) = 1.0236186E-01;   Omega22(3,4,2) = -9.0058935E-01;  Omega22(3,4,3) = 4.4472565E+02;
    //O
    Omega22(4,0,0) = -8.3110691E-03;  Omega22(4,0,1) = 1.9617877E-01;   Omega22(4,0,2) = -1.7205427E+00;  Omega22(4,0,3) = 4.0812829E+03;
    Omega22(4,1,0) = -3.7969686E-03;  Omega22(4,1,1) = 7.6789981E-02;   Omega22(4,1,2) = -7.3056809E-01;  Omega22(4,1,3) = 3.3958171E+02;
    Omega22(4,2,0) = -1.0066279E-03;  Omega22(4,2,1) = 1.1029264E-02;   Omega22(4,2,2) = -2.0671266E-01;  Omega22(4,2,3) = 8.2644384E+01;
    Omega22(4,3,0) = -5.0478143E-03;  Omega22(4,3,1) = 1.0236186E-01;   Omega22(4,3,2) = -9.0058935E-01;  Omega22(4,3,3) = 4.4472565E+02;
    Omega22(4,4,0) = -4.2451096E-03;  Omega22(4,4,1) = 9.6820337E-02;   Omega22(4,4,2) = -9.9770795E-01;  Omega22(4,4,3) = 8.3320644E+02;

    // Creation/Destruction (+1/-1), Index of monoatomic reactants
    // Monoatomic species (N,O) recombine into diaatomic (N2, O2)
    CatRecombTable(0,0) =  1; CatRecombTable(0,1) = 3;
    CatRecombTable(1,0) =  1; CatRecombTable(1,1) = 4;
    CatRecombTable(2,0) =  0; CatRecombTable(2,1) = 0;
    CatRecombTable(3,0) = -1; CatRecombTable(3,1) = 3;
    CatRecombTable(4,0) = -1; CatRecombTable(4,1) = 4;

    /*--- Values used in the Sutherland's formula. ---*/
    if (viscous) {
      //F.M. White, Viscous Fluid Flow, 3rd ed., McGraw-Hill, 2006.
      k_ref[0] = 0.0241;
      mu_ref[0] = 1.716E-5;
      Sm_ref[0] = 111.0;
      Sk_ref[0] = 194.0;
    }

  } else if (gas_model == "AIR-7"){

    /*--- Check for errors in the initialization ---*/
    if (nSpecies != 7) {
      SU2_MPI::Error("CONFIG ERROR: nSpecies mismatch between gas model & gas composition", CURRENT_FUNCTION);
    }

    mf = 0.0;
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      mf += MassFrac_Freestream[iSpecies];
    if (mf != 1.0) {
      SU2_MPI::Error("CONFIG ERROR: Intial gas mass fractions do not sum to 1!", CURRENT_FUNCTION);
    }

    /*--- Define parameters of the gas model ---*/
    gamma       = 1.4;
    nReactions  = 22;
    ionization  = true;

    Reactions.resize(nReactions,2,6,0.0);
    ArrheniusCoefficient.resize(nReactions,0.0);
    ArrheniusEta.resize(nReactions,0.0);
    ArrheniusTheta.resize(nReactions,0.0);
    Tcf_a.resize(nReactions,0.0);
    Tcf_b.resize(nReactions,0.0);
    Tcb_a.resize(nReactions,0.0);
    Tcb_b.resize(nReactions,0.0);

    /*--- Assign gas properties ---*/
    // Rotational modes of energy storage
    RotationModes[0] = 0.0; // e-
    RotationModes[1] = 2.0; // N2
    RotationModes[2] = 2.0; // O2
    RotationModes[3] = 2.0; // NO
    RotationModes[4] = 0.0; // N
    RotationModes[5] = 0.0; // O
    RotationModes[6] = 2.0; // NO+

    // Molar mass [kg/kmol]
    MolarMass[0] = 5.4858E-04;      // e-
    MolarMass[1] = 2.0*14.0067;     // N2
    MolarMass[2] = 2.0*15.9994;     // O2
    MolarMass[3] = 14.0067+15.9994; // NO
    MolarMass[4] = 14.0067;         // N
    MolarMass[5] = 15.9994;         // O
    MolarMass[6] = 14.0067+15.9994; // NO+

    //Characteristic vibrational temperatures
    CharVibTemp[0] = 0.0;    // e-
    CharVibTemp[1] = 3395.0; // N2
    CharVibTemp[2] = 2239.0; // O2
    CharVibTemp[3] = 2817.0; // NO
    CharVibTemp[4] = 0.0;    // N
    CharVibTemp[5] = 0.0;    // O
    CharVibTemp[6] = 2817.0; // NO+

    // Formation enthalpy: (Scalabrin values, J/kg)
    Enthalpy_Formation[0] = 0.0;    // e-
    Enthalpy_Formation[1] = 0.0;    // N2
    Enthalpy_Formation[2] = 0.0;    // O2
    Enthalpy_Formation[3] = 3.0E6;  // NO
    Enthalpy_Formation[4] = 3.36E7; // N
    Enthalpy_Formation[5] = 1.54E7; // O
    Enthalpy_Formation[6] = 3.28E7; // NO+

    // Reference temperature (JANAF values, [K])
    Ref_Temperature[0] = 0.0;
    Ref_Temperature[1] = 0.0;
    Ref_Temperature[2] = 0.0;
    Ref_Temperature[3] = 0.0;
    Ref_Temperature[4] = 0.0;
    Ref_Temperature[5] = 0.0;
    Ref_Temperature[6] = 0.0;

    // Blottner viscosity coefficients
    // A                        // B                        // C
    Blottner(0,0) = 0.00E+0;   Blottner(0,1) =  0.00E+0;  Blottner(0,2) = -1.20E1;  // e-
    Blottner(1,0) = 2.68E-2;   Blottner(1,1) =  3.18E-1;  Blottner(1,2) = -1.13E1;  // N2
    Blottner(2,0) = 4.49E-2;   Blottner(2,1) = -8.26E-2;  Blottner(2,2) = -9.20E0;  // O2
    Blottner(3,0) = 4.36E-2;   Blottner(3,1) = -3.36E-2;  Blottner(3,2) = -9.58E0;  // NO
    Blottner(4,0) = 1.16E-2;   Blottner(4,1) =  6.03E-1;  Blottner(4,2) = -1.24E1;  // N
    Blottner(5,0) = 2.03E-2;   Blottner(5,1) =  4.29E-1;  Blottner(5,2) = -1.16E1;  // O
    Blottner(6,0) = 3.02E-1;   Blottner(6,1) =  -3.50E0;  Blottner(6,2) = -3.74E0;  // NO+

    // Number of electron states
    nElStates[0] = 1;  // e-
    nElStates[1] = 15; // N2
    nElStates[2] = 7;  // O2
    nElStates[3] = 16; // NO
    nElStates[4] = 3;  // N
    nElStates[5] = 5;  // O
    nElStates[6] = 8;  // NO+

    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      maxEl = max(maxEl, nElStates[iSpecies]);

    /*--- Allocate and initialize electron data arrays ---*/
    CharElTemp.resize(nSpecies,maxEl) = su2double(0.0);
    ElDegeneracy.resize(nSpecies,maxEl) = su2double(0.0);

    // e: 1 state
    CharElTemp(0,0) = 0.000000000000000E+00;
    ElDegeneracy(0,0) = 1;

    //N2: 15 states
    CharElTemp(1,0)  = 0.000000000000000E+00;
    CharElTemp(1,1)  = 7.223156514095200E+04;
    CharElTemp(1,2)  = 8.577862640384000E+04;
    CharElTemp(1,3)  = 8.605026716160000E+04;
    CharElTemp(1,4)  = 9.535118627874400E+04;
    CharElTemp(1,5)  = 9.805635702203200E+04;
    CharElTemp(1,6)  = 9.968267656935200E+04;
    CharElTemp(1,7)  = 1.048976467715200E+05;
    CharElTemp(1,8)  = 1.116489555200000E+05;
    CharElTemp(1,9)  = 1.225836470400000E+05;
    CharElTemp(1,10) = 1.248856873600000E+05;
    CharElTemp(1,11) = 1.282476158188320E+05;
    CharElTemp(1,12) = 1.338060936000000E+05;
    CharElTemp(1,13) = 1.404296391107200E+05;
    CharElTemp(1,14) = 1.504958859200000E+05;
    ElDegeneracy(1,0)  = 1;
    ElDegeneracy(1,1)  = 3;
    ElDegeneracy(1,2)  = 6;
    ElDegeneracy(1,3)  = 6;
    ElDegeneracy(1,4)  = 3;
    ElDegeneracy(1,5)  = 1;
    ElDegeneracy(1,6)  = 2;
    ElDegeneracy(1,7)  = 2;
    ElDegeneracy(1,8)  = 5;
    ElDegeneracy(1,9)  = 1;
    ElDegeneracy(1,10) = 6;
    ElDegeneracy(1,11) = 6;
    ElDegeneracy(1,12) = 10;
    ElDegeneracy(1,13) = 6;
    ElDegeneracy(1,14) = 6;
    // O2: 7 states
    CharElTemp(2,0) = 0.000000000000000E+00;
    CharElTemp(2,1) = 1.139156019700800E+04;
    CharElTemp(2,2) = 1.898473947826400E+04;
    CharElTemp(2,3) = 4.755973576639200E+04;
    CharElTemp(2,4) = 4.991242097343200E+04;
    CharElTemp(2,5) = 5.092268575561600E+04;
    CharElTemp(2,6) = 7.189863255967200E+04;
    ElDegeneracy(2,0) = 3;
    ElDegeneracy(2,1) = 2;
    ElDegeneracy(2,2) = 1;
    ElDegeneracy(2,3) = 1;
    ElDegeneracy(2,4) = 6;
    ElDegeneracy(2,5) = 3;
    ElDegeneracy(2,6) = 3;
    // NO: 16 states
    CharElTemp(3,0)  = 0.000000000000000E+00;
    CharElTemp(3,1)  = 5.467345760000000E+04;
    CharElTemp(3,2)  = 6.317139627802400E+04;
    CharElTemp(3,3)  = 6.599450342445600E+04;
    CharElTemp(3,4)  = 6.906120960000000E+04;
    CharElTemp(3,5)  = 7.049998480000000E+04;
    CharElTemp(3,6)  = 7.491055017560000E+04;
    CharElTemp(3,7)  = 7.628875293968000E+04;
    CharElTemp(3,8)  = 8.676188537552000E+04;
    CharElTemp(3,9)  = 8.714431182368000E+04;
    CharElTemp(3,10) = 8.886077063728000E+04;
    CharElTemp(3,11) = 8.981755614528000E+04;
    CharElTemp(3,12) = 8.988445919208000E+04;
    CharElTemp(3,13) = 9.042702132000000E+04;
    CharElTemp(3,14) = 9.064283760000000E+04;
    CharElTemp(3,15) = 9.111763341600000E+04;
    ElDegeneracy(3,0)  = 4;
    ElDegeneracy(3,1)  = 8;
    ElDegeneracy(3,2)  = 2;
    ElDegeneracy(3,3)  = 4;
    ElDegeneracy(3,4)  = 4;
    ElDegeneracy(3,5)  = 4;
    ElDegeneracy(3,6)  = 4;
    ElDegeneracy(3,7)  = 2;
    ElDegeneracy(3,8)  = 4;
    ElDegeneracy(3,9)  = 2;
    ElDegeneracy(3,10) = 4;
    ElDegeneracy(3,11) = 4;
    ElDegeneracy(3,12) = 2;
    ElDegeneracy(3,13) = 2;
    ElDegeneracy(3,14) = 2;
    ElDegeneracy(3,15) = 4;
    // N: 3 states
    CharElTemp(4,0) = 0.000000000000000E+00;
    CharElTemp(4,1) = 2.766469645581980E+04;
    CharElTemp(4,2) = 4.149309313560210E+04;
    ElDegeneracy(4,0)= 4;
    ElDegeneracy(4,1)= 10;
    ElDegeneracy(4,2)= 6;
    // O: 5 states
    CharElTemp(5,0) = 0.000000000000000E+00;
    CharElTemp(5,1) = 2.277077570280000E+02;
    CharElTemp(5,2) = 3.265688785704000E+02;
    CharElTemp(5,3) = 2.283028632262240E+04;
    CharElTemp(5,4) = 4.861993036434160E+04;
    ElDegeneracy(5,0) = 5;
    ElDegeneracy(5,1) = 3;
    ElDegeneracy(5,2) = 1;
    ElDegeneracy(5,3) = 5;
    ElDegeneracy(5,4) = 1;
    // NO+: 8 states
    CharElTemp(6,0) = 0.000000000000000E+00;
    CharElTemp(6,1) = 7.508967768800000E+04;
    CharElTemp(6,2) = 8.525462447600000E+04;
    CharElTemp(6,3) = 8.903572570160000E+04;
    CharElTemp(6,4) = 9.746982592400000E+04;
    CharElTemp(6,5) = 1.000553049584000E+05;
    CharElTemp(6,6) = 1.028033655904000E+05;
    CharElTemp(6,7) = 1.057138639424800E+05;
    ElDegeneracy(6,0) = 1;
    ElDegeneracy(6,1) = 3;
    ElDegeneracy(6,2) = 6;
    ElDegeneracy(6,3) = 6;
    ElDegeneracy(6,4) = 3;
    ElDegeneracy(6,5) = 1;
    ElDegeneracy(6,6) = 2;
    ElDegeneracy(6,7) = 2;

    /*--- Set reaction maps ---*/
    // N2 dissociation
    Reactions(0,0,0)=1;    Reactions(0,0,1)=1;   Reactions(0,0,2)=nSpecies;    Reactions(0,1,0)=4;   Reactions(0,1,1)=4;   Reactions(0,1,2) =1;
    Reactions(1,0,0)=1;    Reactions(1,0,1)=2;   Reactions(1,0,2)=nSpecies;    Reactions(1,1,0)=4;   Reactions(1,1,1)=4;   Reactions(1,1,2) =2;
    Reactions(2,0,0)=1;    Reactions(2,0,1)=3;   Reactions(2,0,2)=nSpecies;    Reactions(2,1,0)=4;   Reactions(2,1,1)=4;   Reactions(2,1,2) =3;
    Reactions(3,0,0)=1;    Reactions(3,0,1)=4;   Reactions(3,0,2)=nSpecies;    Reactions(3,1,0)=4;   Reactions(3,1,1)=4;   Reactions(3,1,2) =4;
    Reactions(4,0,0)=1;    Reactions(4,0,1)=5;   Reactions(4,0,2)=nSpecies;    Reactions(4,1,0)=4;   Reactions(4,1,1)=4;   Reactions(4,1,2) =5;
    Reactions(5,0,0)=1;    Reactions(5,0,1)=6;   Reactions(5,0,2)=nSpecies;    Reactions(5,1,0)=4;   Reactions(5,1,1)=4;   Reactions(5,1,2) =6;
    // O2 dissociation
    Reactions(6,0,0)=2;    Reactions(6,0,1)=1;   Reactions(6,0,2)=nSpecies;    Reactions(6,1,0)=5;   Reactions(6,1,1)=5;   Reactions(6,1,2) =1;
    Reactions(7,0,0)=2;    Reactions(7,0,1)=2;   Reactions(7,0,2)=nSpecies;    Reactions(7,1,0)=5;   Reactions(7,1,1)=5;   Reactions(7,1,2) =2;
    Reactions(8,0,0)=2;    Reactions(8,0,1)=3;   Reactions(8,0,2)=nSpecies;    Reactions(8,1,0)=5;   Reactions(8,1,1)=5;   Reactions(8,1,2) =3;
    Reactions(9,0,0)=2;    Reactions(9,0,1)=4;   Reactions(9,0,2)=nSpecies;    Reactions(9,1,0)=5;   Reactions(9,1,1)=5;   Reactions(9,1,2) =4;
    Reactions(10,0,0)=2;   Reactions(10,0,1)=5;  Reactions(10,0,2)=nSpecies;   Reactions(10,1,0)=5;  Reactions(10,1,1)=5;  Reactions(10,1,2) =5;
    Reactions(11,0,0)=2;   Reactions(11,0,1)=6;  Reactions(11,0,2)=nSpecies;   Reactions(11,1,0)=5;  Reactions(11,1,1)=5;  Reactions(11,1,2) =6;
    // NO dissociation
    Reactions(12,0,0)=3;   Reactions(12,0,1)=1;  Reactions(12,0,2)=nSpecies;   Reactions(12,1,0)=4;  Reactions(12,1,1)=5;  Reactions(12,1,2) =1;
    Reactions(13,0,0)=3;   Reactions(13,0,1)=2;  Reactions(13,0,2)=nSpecies;   Reactions(13,1,0)=4;  Reactions(13,1,1)=5;  Reactions(13,1,2) =2;
    Reactions(14,0,0)=3;   Reactions(14,0,1)=3;  Reactions(14,0,2)=nSpecies;   Reactions(14,1,0)=4;  Reactions(14,1,1)=5;  Reactions(14,1,2) =3;
    Reactions(15,0,0)=3;   Reactions(15,0,1)=4;  Reactions(15,0,2)=nSpecies;   Reactions(15,1,0)=4;  Reactions(15,1,1)=5;  Reactions(15,1,2) =4;
    Reactions(16,0,0)=3;   Reactions(16,0,1)=5;  Reactions(16,0,2)=nSpecies;   Reactions(16,1,0)=4;  Reactions(16,1,1)=5;  Reactions(16,1,2) =5;
    Reactions(17,0,0)=3;   Reactions(17,0,1)=6;  Reactions(17,0,2)=nSpecies;   Reactions(17,1,0)=4;  Reactions(17,1,1)=5;  Reactions(17,1,2) =6;
    // N2 + O -> NO + N
    Reactions(18,0,0)=1;   Reactions(18,0,1)=5;  Reactions(18,0,2)=nSpecies;   Reactions(18,1,0)=3;  Reactions(18,1,1)=4;  Reactions(18,1,2)= nSpecies;
    // NO + O -> O2 + N
    Reactions(19,0,0)=3;   Reactions(19,0,1)=5;  Reactions(19,0,2)=nSpecies;   Reactions(19,1,0)=2;  Reactions(19,1,1)=4;  Reactions(19,1,2)= nSpecies;
    //N + O -> NO+ + e
    Reactions(20,0,0)=4;   Reactions(20,0,1)=5;  Reactions(20,0,2)=nSpecies;   Reactions(20,1,0)=6;  Reactions(20,1,1)=0;  Reactions(20,1,2)= nSpecies;
    //N2 + e -> N + N + e
    Reactions(21,0,0)=1;   Reactions(21,0,1)=0;  Reactions(21,0,2)=nSpecies;   Reactions(21,1,0)=4;  Reactions(21,1,1)=4;  Reactions(21,1,2)= 0;

    /*--- Set Arrhenius coefficients for reactions ---*/
    // Pre-exponential factor
    ArrheniusCoefficient[0]  = 7.0E21;
    ArrheniusCoefficient[1]  = 7.0E21;
    ArrheniusCoefficient[2]  = 7.0E21;
    ArrheniusCoefficient[3]  = 3.0E22;
    ArrheniusCoefficient[4]  = 3.0E22;
    ArrheniusCoefficient[5]  = 7.0E21;
    ArrheniusCoefficient[6]  = 2.0E21;
    ArrheniusCoefficient[7]  = 2.0E21;
    ArrheniusCoefficient[8]  = 2.0E21;
    ArrheniusCoefficient[9]  = 1.0E22;
    ArrheniusCoefficient[10] = 1.0E22;
    ArrheniusCoefficient[11] = 2.0E21;
    ArrheniusCoefficient[12] = 5.0E15;
    ArrheniusCoefficient[13] = 5.0E15;
    ArrheniusCoefficient[14] = 5.0E15;
    ArrheniusCoefficient[15] = 1.1E17;
    ArrheniusCoefficient[16] = 1.1E17;
    ArrheniusCoefficient[17] = 5.0E15;
    ArrheniusCoefficient[18] = 6.4E17;
    ArrheniusCoefficient[19] = 8.4E12;
    ArrheniusCoefficient[20] = 5.3E12;
    ArrheniusCoefficient[21] = 3.0E24;

    // Rate-controlling temperature exponent
    ArrheniusEta[0]  = -1.60;
    ArrheniusEta[1]  = -1.60;
    ArrheniusEta[2]  = -1.60;
    ArrheniusEta[3]  = -1.60;
    ArrheniusEta[4]  = -1.60;
    ArrheniusEta[5]  = -1.60;
    ArrheniusEta[6]  = -1.50;
    ArrheniusEta[7]  = -1.50;
    ArrheniusEta[8]  = -1.50;
    ArrheniusEta[9]  = -1.50;
    ArrheniusEta[10] = -1.50;
    ArrheniusEta[11] = -1.50;
    ArrheniusEta[12] = 0.0;
    ArrheniusEta[13] = 0.0;
    ArrheniusEta[14] = 0.0;
    ArrheniusEta[15] = 0.0;
    ArrheniusEta[16] = 0.0;
    ArrheniusEta[17] = 0.0;
    ArrheniusEta[18] = -1.0;
    ArrheniusEta[19] = 0.0;
    ArrheniusEta[20] = 0.0;
    ArrheniusEta[21] = -1.60;

    // Characteristic temperature
    ArrheniusTheta[0]  = 113200.0;
    ArrheniusTheta[1]  = 113200.0;
    ArrheniusTheta[2]  = 113200.0;
    ArrheniusTheta[3]  = 113200.0;
    ArrheniusTheta[4]  = 113200.0;
    ArrheniusTheta[5]  = 113200.0;
    ArrheniusTheta[6]  = 59500.0;
    ArrheniusTheta[7]  = 59500.0;
    ArrheniusTheta[8]  = 59500.0;
    ArrheniusTheta[9]  = 59500.0;
    ArrheniusTheta[10]  = 59500.0;
    ArrheniusTheta[11]  = 59500.0;
    ArrheniusTheta[12] = 75500.0;
    ArrheniusTheta[13] = 75500.0;
    ArrheniusTheta[14] = 75500.0;
    ArrheniusTheta[15] = 75500.0;
    ArrheniusTheta[16] = 75500.0;
    ArrheniusTheta[17] = 75500.0;
    ArrheniusTheta[18] = 38400.0;
    ArrheniusTheta[19] = 19450.0;
    ArrheniusTheta[20] = 31900.0;
    ArrheniusTheta[21] = 113200.0;

    /*--- Set rate-controlling temperature exponents ---*/
    //  -----------  Tc = Ttr^a * Tve^b  -----------
    //
    // Forward Reactions
    //   Dissociation:         a = 0.5, b = 0.5  (OR a = 0.7, b =0.3)
    //   Exchange:             a = 1,   b = 0
    //   Associative ion...    a = 1,   b = 0  ???
    //   E Impact dissociation a = 0,   b = 1
    //   E Impact ionization:  a = 0,   b = 1
    //
    // Backward Reactions
    //   Dissociation:           a = 1,   b = 0
    //   Exchange:               a = 1,   b = 0
    //   Associative  ion...     a = 0.5, b = 0.5
    //   E Impact ionization:    a = 0,   b = 1
    //   E Impact dissocitation: a = 0.5, b = 0.5 ???
    //   N2 impact dissociation: a = 0,   b = 1
    //   Others:                 a = 1,   b = 0
    Tcf_a[0]  = 0.5; Tcf_b[0]  = 0.5; Tcb_a[0]  = 1;   Tcb_b[0] = 0;
    Tcf_a[1]  = 0.5; Tcf_b[1]  = 0.5; Tcb_a[1]  = 1;   Tcb_b[1] = 0;
    Tcf_a[2]  = 0.5; Tcf_b[2]  = 0.5; Tcb_a[2]  = 1;   Tcb_b[2] = 0;
    Tcf_a[3]  = 0.5; Tcf_b[3]  = 0.5; Tcb_a[3]  = 1;   Tcb_b[3] = 0;
    Tcf_a[4]  = 0.5; Tcf_b[4]  = 0.5; Tcb_a[4]  = 1;   Tcb_b[4] = 0;
    Tcf_a[5]  = 0.5; Tcf_b[5]  = 0.5; Tcb_a[5]  = 1;   Tcb_b[5] = 0;
    Tcf_a[6]  = 0.5; Tcf_b[6]  = 0.5; Tcb_a[6]  = 1;   Tcb_b[6] = 0;
    Tcf_a[7]  = 0.5; Tcf_b[7]  = 0.5; Tcb_a[7]  = 1;   Tcb_b[7] = 0;
    Tcf_a[8]  = 0.5; Tcf_b[8]  = 0.5; Tcb_a[8]  = 1;   Tcb_b[8] = 0;
    Tcf_a[9]  = 0.5; Tcf_b[9]  = 0.5; Tcb_a[9]  = 1;   Tcb_b[9] = 0;
    Tcf_a[10] = 0.5; Tcf_b[10] = 0.5; Tcb_a[10] = 1;   Tcb_b[10] = 0;
    Tcf_a[11] = 0.5; Tcf_b[11] = 0.5; Tcb_a[11] = 1;   Tcb_b[11] = 0;
    Tcf_a[12] = 0.5; Tcf_b[12] = 0.5; Tcb_a[12] = 1;   Tcb_b[12] = 0;
    Tcf_a[13] = 0.5; Tcf_b[13] = 0.5; Tcb_a[13] = 1;   Tcb_b[13] = 0;
    Tcf_a[14] = 0.5; Tcf_b[14] = 0.5; Tcb_a[14] = 1;   Tcb_b[14] = 0;
    Tcf_a[15] = 0.5; Tcf_b[15] = 0.5; Tcb_a[15] = 1;   Tcb_b[15] = 0;
    Tcf_a[16] = 0.5; Tcf_b[16] = 0.5; Tcb_a[16] = 1;   Tcb_b[16] = 0;
    Tcf_a[17] = 0.5; Tcf_b[17] = 0.5; Tcb_a[17] = 1;   Tcb_b[17] = 0;
    Tcf_a[18] = 1.0; Tcf_b[18] = 0.0; Tcb_a[18] = 1;   Tcb_b[18] = 0;
    Tcf_a[19] = 1.0; Tcf_b[19] = 0.0; Tcb_a[19] = 1;   Tcb_b[19] = 0;
    Tcf_a[20] = 1.0; Tcf_b[20] = 0.0; Tcb_a[20] = 0.5; Tcb_b[20] = 0.5;
    Tcf_a[21] = 0.0; Tcf_b[21] = 1.0; Tcb_a[21] = 0;   Tcb_b[21] = 1;

    /*--- Collision integral data ---*/
    // Index 1: collider
    // Index 2: partner
    // Index 3: A1, A2, A3

    // Omega^(1,1) ----------------------
    Omega11(0,0,0) = -1.000000E+00;  Omega11(0,0,1) = -1.000000E+00;  Omega11(0,0,2) = -1.000000E+00;  Omega11(0,0,3) = -1.000000E+00;
    Omega11(0,1,0) = -1.0525124E-02; Omega11(0,1,1) = 1.3498950E-01;  Omega11(0,1,2) = 1.2524805E-01;  Omega11(0,1,3) = 1.5066506E-01;
    Omega11(0,2,0) = 2.3527001E-02;  Omega11(0,2,1) = -6.9632323E-01; Omega11(0,2,2) = 6.8035475E+00;  Omega11(0,2,3) = 1.8335509E-09;
    Omega11(0,3,0) = 1.0414818E-01;  Omega11(0,3,1) = -2.8369126E+00; Omega11(0,3,2) = 2.5323135E+01;  Omega11(0,3,3) = 7.7138358E-32;
    Omega11(0,4,0) = 0.0000000E+00;  Omega11(0,4,1) = 1.6554247E-01;  Omega11(0,4,2) = -3.4986344E+00; Omega11(0,4,3) = 5.9268038E+08;
    Omega11(0,5,0) = 9.9865506E-03;  Omega11(0,5,1) = -2.7407431E-01; Omega11(0,5,2) = 2.6561032E+00;  Omega11(0,5,3) = 4.3080676E-04;
    Omega11(0,6,0) = 1.0000000E+00;  Omega11(0,6,1) = 1.0000000E+00;  Omega11(0,6,2) = 1.0000000E+00;  Omega11(0,6,3) = 1.0000000E+00;
    //N2
    Omega11(1,0,0) = -1.0525124E-02; Omega11(1,0,1) = 1.3498950E-01;  Omega11(1,0,2) = 1.2524805E-01;  Omega11(1,0,3) = 1.5066506E-01;
    Omega11(1,1,0) = -6.0614558E-03; Omega11(1,1,1) = 1.2689102E-01;  Omega11(1,1,2) = -1.0616948E+00; Omega11(1,1,3) = 8.0955466E+02;
    Omega11(1,2,0) = -3.7959091E-03; Omega11(1,2,1) = 9.5708295E-02;  Omega11(1,2,2) = -1.0070611E+00; Omega11(1,2,3) = 8.9392313E+02;
    Omega11(1,3,0) = -1.9295666E-03; Omega11(1,3,1) = 2.7995735E-02;  Omega11(1,3,2) = -3.1588514E-01; Omega11(1,3,3) = 1.2880734E+02;
    Omega11(1,4,0) = -1.0796249E-02; Omega11(1,4,1) = 2.2656509E-01;  Omega11(1,4,2) = -1.7910602E+00; Omega11(1,4,3) = 4.0455218E+03;
    Omega11(1,5,0) = -2.7244269E-03; Omega11(1,5,1) = 6.9587171E-02;  Omega11(1,5,2) = -7.9538667E-01; Omega11(1,5,3) = 4.0673730E+02;
    Omega11(1,6,0) = 0.0000000E+00;  Omega11(1,6,1) = 9.1205839E-02;  Omega11(1,6,2) = -1.8728231E+00; Omega11(1,6,3) = 2.4432020E+05;
    //O2
    Omega11(2,0,0) = 2.3527001E-02;  Omega11(2,0,1) = -6.9632323E-01; Omega11(2,0,2) = 6.8035475E+00;  Omega11(2,0,3) = 1.8335509E-09;
    Omega11(2,1,0) = -3.7959091E-03; Omega11(2,1,1) = 9.5708295E-02;  Omega11(2,1,2) = -1.0070611E+00; Omega11(2,1,3) = 8.9392313E+02;
    Omega11(2,2,0) = -8.0682650E-04; Omega11(2,2,1) = 1.6602480E-02;  Omega11(2,2,2) = -3.1472774E-01; Omega11(2,2,3) = 1.4116458E+02;
    Omega11(2,3,0) = -6.4433840E-04; Omega11(2,3,1) = 8.5378580E-03;  Omega11(2,3,2) = -2.3225102E-01; Omega11(2,3,3) = 1.1371608E+02;
    Omega11(2,4,0) = -1.1453028E-03; Omega11(2,4,1) = 1.2654140E-02;  Omega11(2,4,2) = -2.2435218E-01; Omega11(2,4,3) = 7.7201588E+01;
    Omega11(2,5,0) = -4.8405803E-03; Omega11(2,5,1) = 1.0297688E-01;  Omega11(2,5,2) = -9.6876576E-01; Omega11(2,5,3) = 6.1629812E+02;
    Omega11(2,6,0) = -3.7822765E-03; Omega11(2,6,1) = 1.7967016E-01;  Omega11(2,6,2) = -2.5409098E+00; Omega11(2,6,3) = 1.1840435E+06;
    //NO
    Omega11(3,0,0) = 1.0414818E-01;  Omega11(3,0,1) = -2.8369126E+00; Omega11(3,0,2) = 2.5323135E+01;  Omega11(3,0,3) = 7.7138358E-32;
    Omega11(3,1,0) = -1.9295666E-03; Omega11(3,1,1) = 2.7995735E-02;  Omega11(3,1,2) = -3.1588514E-01; Omega11(3,1,3) = 1.2880734E+02;
    Omega11(3,2,0) = -6.4433840E-04; Omega11(3,2,1) = 8.5378580E-03;  Omega11(3,2,2) = -2.3225102E-01; Omega11(3,2,3) = 1.1371608E+02;
    Omega11(3,3,0) = -0.0000000E+00; Omega11(3,3,1) = -1.1056066E-02; Omega11(3,3,2) = -5.9216250E-02; Omega11(3,3,3) = 7.2542367E+01;
    Omega11(3,4,0) = -1.5770918E-03; Omega11(3,4,1) = 1.9578381E-02;  Omega11(3,4,2) = -2.7873624E-01; Omega11(3,4,3) = 9.9547944E+01;
    Omega11(3,5,0) = -1.0885815E-03; Omega11(3,5,1) = 1.1883688E-02;  Omega11(3,5,2) = -2.1844909E-01; Omega11(3,5,3) = 7.5512560E+01;
    Omega11(3,6,0) = -8.1158474E-03; Omega11(3,6,1) = 2.1474280E-01;  Omega11(3,6,2) = -2.0148450E+00; Omega11(3,6,3) = 6.2986385E+04;
    //N
    Omega11(4,0,0) = 0.0000000E+00;  Omega11(4,0,1) = 1.6554247E-01;  Omega11(4,0,2) = -3.4986344E+00; Omega11(4,0,3) = 5.9268038E+08;
    Omega11(4,1,0) = -1.0796249E-02; Omega11(4,1,1) = 2.2656509E-01;  Omega11(4,1,2) = -1.7910602E+00; Omega11(4,1,3) = 4.0455218E+03;
    Omega11(4,2,0) = -1.1453028E-03; Omega11(4,2,1) = 1.2654140E-02;  Omega11(4,2,2) = -2.2435218E-01; Omega11(4,2,3) = 7.7201588E+01;
    Omega11(4,3,0) = -1.5770918E-03; Omega11(4,3,1) = 1.9578381E-02;  Omega11(4,3,2) = -2.7873624E-01; Omega11(4,3,3) = 9.9547944E+01;
    Omega11(4,4,0) = -9.6083779E-03; Omega11(4,4,1) = 2.0938971E-01;  Omega11(4,4,2) = -1.7386904E+00; Omega11(4,4,3) = 3.3587983E+03;
    Omega11(4,5,0) = -7.8147689E-03; Omega11(4,5,1) = 1.6792705E-01;  Omega11(4,5,2) = -1.4308628E+00; Omega11(4,5,3) = 1.6628859E+03;
    Omega11(4,6,0) = -1.9605234E-02; Omega11(4,6,1) = 5.5570872E-01;  Omega11(4,6,2) = -5.4285702E+00; Omega11(4,6,3) = 1.3574446E+09;
    //O
    Omega11(5,0,0) = 9.9865506E-03;  Omega11(5,0,1) = -2.7407431E-01; Omega11(5,0,2) = 2.6561032E+00;  Omega11(5,0,3) = 4.3080676E-04;
    Omega11(5,1,0) = -2.7244269E-03; Omega11(5,1,1) = 6.9587171E-02;  Omega11(5,1,2) = -7.9538667E-01; Omega11(5,1,3) = 4.0673730E+02;
    Omega11(5,2,0) = -4.8405803E-03; Omega11(5,2,1) = 1.0297688E-01;  Omega11(5,2,2) = -9.6876576E-01; Omega11(5,2,3) = 6.1629812E+02;
    Omega11(5,3,0) = -1.0885815E-03; Omega11(5,3,1) = 1.1883688E-02;  Omega11(5,3,2) = -2.1844909E-01; Omega11(5,3,3) = 7.5512560E+01;
    Omega11(5,4,0) = -7.8147689E-03; Omega11(5,4,1) = 1.6792705E-01;  Omega11(5,4,2) = -1.4308628E+00; Omega11(5,4,3) = 1.6628859E+03;
    Omega11(5,5,0) = -6.4040535E-03; Omega11(5,5,1) = 1.4629949E-01;  Omega11(5,5,2) = -1.3892121E+00; Omega11(5,5,3) = 2.0903441E+03;
    Omega11(5,6,0) = -1.6409054E-02; Omega11(5,6,1) = 4.6352852E-01;  Omega11(5,6,2) = -4.5479735E+00; Omega11(5,6,3) = 7.4250671E+07;
    //NO+
    Omega11(6,0,0) = 1.0000000E+00;  Omega11(6,0,1) = 1.0000000E+00;  Omega11(6,0,2) = 1.0000000E+00;  Omega11(6,0,3) = 1.0000000E+00;
    Omega11(6,1,0) = 0.0000000E+00;  Omega11(6,1,1) = 9.1205839E-02;  Omega11(6,1,2) = -1.8728231E+00; Omega11(6,1,3) = 2.4432020E+05;
    Omega11(6,2,0) = -3.7822765E-03; Omega11(6,2,1) = 1.7967016E-01;  Omega11(6,2,2) = -2.5409098E+00; Omega11(6,2,3) = 1.1840435E+06;
    Omega11(6,3,0) = -8.1158474E-03; Omega11(6,3,1) = 2.1474280E-01;  Omega11(6,3,2) = -2.0148450E+00; Omega11(6,3,3) = 6.2986385E+04;
    Omega11(6,4,0) = -1.9605234E-02; Omega11(6,4,1) = 5.5570872E-01;  Omega11(6,4,2) = -5.4285702E+00; Omega11(6,4,3) = 1.3574446E+09;
    Omega11(6,5,0) = -1.6409054E-02; Omega11(6,5,1) = 4.6352852E-01;  Omega11(6,5,2) = -4.5479735E+00; Omega11(6,5,3) = 7.4250671E+07;
    Omega11(6,6,0) = -1.000000E+00;  Omega11(6,6,1) = -1.000000E+00;  Omega11(6,6,2) = -1.000000E+00;  Omega11(6,6,3) = -1.000000E+00;

    // Omega^(2,2) ----------------------
    Omega22(0,0,0) = -1.000000E+00;  Omega22(0,0,1) = -1.000000E+00;  Omega22(0,0,2) = -1.000000E+00;  Omega22(0,0,3) = -1.000000E+00;
    Omega22(0,1,0) = -4.2254948E-03; Omega22(0,1,1) = -5.2965163E-02; Omega22(0,1,2) = 1.9157708E+00;  Omega22(0,1,3) = 6.3263309E-04;
    Omega22(0,2,0) = 9.6744867E-03;  Omega22(0,2,1) = -3.3759583E-01; Omega22(0,2,2) = 3.7952121E+00;  Omega22(0,2,3) = 6.8468036E-06;
    Omega22(0,3,0) = 0.0000000E+00;  Omega22(0,3,1) = 5.4444485E-02;  Omega22(0,3,2) = -1.2854128E+00; Omega22(0,3,3) = 1.3857556E+04;
    Omega22(0,4,0) = -1.0903638E-01; Omega22(0,4,1) = 2.8678381E+00;  Omega22(0,4,2) = -2.5297550E+01; Omega22(0,4,3) = 3.4838798E+33;
    Omega22(0,5,0) = -1.7924100E-02; Omega22(0,5,1) = 4.0402656E-01;  Omega22(0,5,2) = -2.6712374E+00; Omega22(0,5,3) = 4.1447669E+02;
    Omega22(0,6,0) = 1.0000000E+00;  Omega22(0,6,1) = 1.0000000E+00;  Omega22(0,6,2) = 1.0000000E+00;  Omega22(0,6,3) = 1.0000000E+00;
    //N2
    Omega22(1,0,0) = -4.2254948E-03; Omega22(1,0,1) = -5.2965163E-02; Omega22(1,0,2) = 1.9157708E+00;  Omega22(1,0,3) = 6.3263309E-04;
    Omega22(1,1,0) = -7.6303990E-03; Omega22(1,1,1) = 1.6878089E-01;  Omega22(1,1,2) = -1.4004234E+00; Omega22(1,1,3) = 2.1427708E+03;
    Omega22(1,2,0) = -8.0457321E-03; Omega22(1,2,1) = 1.9228905E-01;  Omega22(1,2,2) = -1.7102854E+00; Omega22(1,2,3) = 5.2213857E+03;
    Omega22(1,3,0) = -6.8237776E-03; Omega22(1,3,1) = 1.4360616E-01;  Omega22(1,3,2) = -1.1922240E+00; Omega22(1,3,3) = 1.2433086E+03;
    Omega22(1,4,0) = -8.3493693E-03; Omega22(1,4,1) = 1.7808911E-01;  Omega22(1,4,2) = -1.4466155E+00; Omega22(1,4,3) = 1.9324210E+03;
    Omega22(1,5,0) = -8.3110691E-03; Omega22(1,5,1) = 1.9617877E-01;  Omega22(1,5,2) = -1.7205427E+00; Omega22(1,5,3) = 4.0812829E+03;
    Omega22(1,6,0) = 0.0000000E+00;  Omega22(1,6,1) = 8.5112236E-02;  Omega22(1,6,2) = -1.7460044E+00; Omega22(1,6,3) = 1.4498969E+05;
    //O2
    Omega22(2,0,0) = 9.6744867E-03;  Omega22(2,0,1) = -3.3759583E-01; Omega22(2,0,2) = 3.7952121E+00;  Omega22(2,0,3) = 6.8468036E-06;
    Omega22(2,1,0) = -8.0457321E-03; Omega22(2,1,1) = 1.9228905E-01;  Omega22(2,1,2) = -1.7102854E+00; Omega22(2,1,3) = 5.2213857E+03;
    Omega22(2,2,0) = -6.2931612E-03; Omega22(2,2,1) = 1.4624645E-01;  Omega22(2,2,2) = -1.3006927E+00; Omega22(2,2,3) = 1.8066892E+03;
    Omega22(2,3,0) = -6.8508672E-03; Omega22(2,3,1) = 1.5524564E-01;  Omega22(2,3,2) = -1.3479583E+00; Omega22(2,3,3) = 2.0037890E+03;
    Omega22(2,4,0) = -1.0608832E-03; Omega22(2,4,1) = 1.1782595E-02;  Omega22(2,4,2) = -2.1246301E-01; Omega22(2,4,3) = 8.4561598E+01;
    Omega22(2,5,0) = -3.7969686E-03; Omega22(2,5,1) = 7.6789981E-02;  Omega22(2,5,2) = -7.3056809E-01; Omega22(2,5,3) = 3.3958171E+02;
    Omega22(2,6,0) = 0.0000000E+00;  Omega22(2,6,1) = 8.4737359E-02;  Omega22(2,6,2) = -1.7290488E+00; Omega22(2,6,3) = 1.2485194E+05;
    //NO
    Omega22(3,0,0) = 0.0000000E+00;  Omega22(0,3,1) = 5.4444485E-02;  Omega22(0,3,2) = -1.2854128E+00; Omega22(0,3,3) = 1.3857556E+04;
    Omega22(3,1,0) = -6.8237776E-03; Omega22(3,1,1) = 1.4360616E-01;  Omega22(3,1,2) = -1.1922240E+00; Omega22(3,1,3) = 1.2433086E+03;
    Omega22(3,2,0) = -6.8508672E-03; Omega22(3,2,1) = 1.5524564E-01;  Omega22(3,2,2) = -1.3479583E+00; Omega22(3,2,3) = 2.0037890E+03;
    Omega22(3,3,0) = -7.4942466E-03; Omega22(3,3,1) = 1.6626193E-01;  Omega22(3,3,2) = -1.4107027E+00; Omega22(3,3,3) = 2.3097604E+03;
    Omega22(3,4,0) = -1.4719259E-03; Omega22(3,4,1) = 1.8446968E-02;  Omega22(3,4,2) = -2.6460411E-01; Omega22(3,4,3) = 1.0911124E+02;
    Omega22(3,5,0) = -1.0066279E-03; Omega22(3,5,1) = 1.1029264E-02;  Omega22(3,5,2) = -2.0671266E-01; Omega22(3,5,3) = 8.2644384E+01;
    Omega22(3,6,0) = 1.1055777E-02;  Omega22(3,6,1) = -1.6621846E-01; Omega22(3,6,2) = 1.4372166E-01;  Omega22(3,6,3) = 1.3182061E+03;
    //N
    Omega22(4,0,0) = -1.0903638E-01; Omega22(4,0,1) = 2.8678381E+00;  Omega22(4,0,2) = -2.5297550E+01; Omega22(4,0,3) = 3.4838798E+33;
    Omega22(4,1,0) = -8.3493693E-03; Omega22(4,1,1) = 1.7808911E-01;  Omega22(4,1,2) = -1.4466155E+00; Omega22(4,1,3) = 1.9324210E+03;
    Omega22(4,2,0) = -1.0608832E-03; Omega22(4,2,1) = 1.1782595E-02;  Omega22(4,2,2) = -2.1246301E-01; Omega22(4,2,3) = 8.4561598E+01;
    Omega22(4,3,0) = -1.4719259E-03; Omega22(4,3,1) = 1.8446968E-02;  Omega22(4,3,2) = -2.6460411E-01; Omega22(4,3,3) = 1.0911124E+02;
    Omega22(4,4,0) = -7.7439615E-03; Omega22(4,4,1) = 1.7129007E-01;  Omega22(4,4,2) = -1.4809088E+00; Omega22(4,4,3) = 2.1284951E+03;
    Omega22(4,5,0) = -5.0478143E-03; Omega22(4,5,1) = 1.0236186E-01;  Omega22(4,5,2) = -9.0058935E-01; Omega22(4,5,3) = 4.4472565E+02;
    Omega22(4,6,0) = -2.1009546E-02; Omega22(4,6,1) = 5.8910426E-01;  Omega22(4,6,2) = -5.6681361E+00; Omega22(4,6,3) = 2.4486594E+09;
    //O
    Omega22(5,0,0) = -1.7924100E-02; Omega22(5,0,1) = 4.0402656E-01;  Omega22(5,0,2) = -2.6712374E+00; Omega22(5,0,3) = 4.1447669E+02;
    Omega22(5,1,0) = -8.3110691E-03; Omega22(5,1,1) = 1.9617877E-01;  Omega22(5,1,2) = -1.7205427E+00; Omega22(5,1,3) = 4.0812829E+03;
    Omega22(5,2,0) = -3.7969686E-03; Omega22(5,2,1) = 7.6789981E-02;  Omega22(5,2,2) = -7.3056809E-01; Omega22(5,2,3) = 3.3958171E+02;
    Omega22(5,3,0) = -1.0066279E-03; Omega22(5,3,1) = 1.1029264E-02;  Omega22(5,3,2) = -2.0671266E-01; Omega22(5,3,3) = 8.2644384E+01;
    Omega22(5,4,0) = -5.0478143E-03; Omega22(5,4,1) = 1.0236186E-01;  Omega22(5,4,2) = -9.0058935E-01; Omega22(5,4,3) = 4.4472565E+02;
    Omega22(5,5,0) = -4.2451096E-03; Omega22(5,5,1) = 9.6820337E-02;  Omega22(5,5,2) = -9.9770795E-01; Omega22(5,5,3) = 8.3320644E+02;
    Omega22(5,6,0) = -1.5315132E-02; Omega22(5,6,1) = 4.3541627E-01;  Omega22(5,6,2) = -4.2864279E+00; Omega22(5,6,3) = 3.5125207E+07;
    //NO+
    Omega22(6,0,0) = 1.0000000E+00;  Omega22(6,0,1) = 1.0000000E+00;  Omega22(6,0,2) = 1.0000000E+00;  Omega22(6,0,3) = 1.0000000E+00;
    Omega22(6,1,0) = 0.0000000E+00;  Omega22(6,1,1) = 8.5112236E-02;  Omega22(6,1,2) = -1.7460044E+00; Omega22(6,1,3) = 1.4498969E+05;
    Omega22(6,2,0) = 0.0000000E+00;  Omega22(6,2,1) = 8.4737359E-02;  Omega22(6,2,2) = -1.7290488E+00; Omega22(6,2,3) = 1.2485194E+05;
    Omega22(6,3,0) = 1.1055777E-02;  Omega22(6,3,1) = -1.6621846E-01; Omega22(6,3,2) = 1.4372166E-01;  Omega22(6,3,3) = 1.3182061E+03;
    Omega22(6,4,0) = -2.1009546E-02; Omega22(6,4,1) = 5.8910426E-01;  Omega22(6,4,2) = -5.6681361E+00; Omega22(6,4,3) = 2.4486594E+09;
    Omega22(6,5,0) = -1.5315132E-02; Omega22(6,5,1) = 4.3541627E-01;  Omega22(6,5,2) = -4.2864279E+00; Omega22(6,5,3) = 3.5125207E+07;
    Omega22(6,6,0) = -1.000000E+00;  Omega22(6,6,1) = -1.000000E+00;  Omega22(6,6,2) = -1.000000E+00;  Omega22(6,6,3) = -1.000000E+00;

    // Creation/Destruction (+1/-1), Index of monoatomic reactants
    // Monoatomic species (N,O) recombine into diaatomic (N2, O2)
    CatRecombTable(0,0) =  0; CatRecombTable(0,1) = 1;
    CatRecombTable(1,0) =  1; CatRecombTable(1,1) = 4;
    CatRecombTable(2,0) =  1; CatRecombTable(2,1) = 5;
    CatRecombTable(3,0) =  0; CatRecombTable(3,1) = 1;
    CatRecombTable(4,0) = -1; CatRecombTable(4,1) = 4;
    CatRecombTable(5,0) = -1; CatRecombTable(5,1) = 5;
    CatRecombTable(6,0) =  0; CatRecombTable(6,1) = 1;

    /*--- Values for Sutherland's formula. ---*/
    if (viscous) {
      //F.M. White, Viscous Fluid Flow, 3rd ed., McGraw-Hill, 2006.
      k_ref[0] = 0.0241;
      mu_ref[0] = 1.716E-5;
      Sm_ref[0] = 111.0;
      Sk_ref[0] = 194.0;
    }
  } else if (gas_model == "AIR-11"){

    /*--- Check for errors in the initialization ---*/
    if (nSpecies != 11) {
      SU2_MPI::Error("CONFIG ERROR: nSpecies mismatch between gas model & gas composition", CURRENT_FUNCTION);
    }

    mf = 0.0;
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      mf += MassFrac_Freestream[iSpecies];
    if (mf != 1.0) {
      SU2_MPI::Error("CONFIG ERROR: Intial gas mass fractions do not sum to 1!", CURRENT_FUNCTION);
    }

    /*--- Define parameters of the gas model ---*/
    gamma       = 1.4;
    nReactions  = 47;
    ionization  = true;

    Reactions.resize(nReactions,2,6,0.0);
    ArrheniusCoefficient.resize(nReactions,0.0);
    ArrheniusEta.resize(nReactions,0.0);
    ArrheniusTheta.resize(nReactions,0.0);
    Tcf_a.resize(nReactions,0.0);
    Tcf_b.resize(nReactions,0.0);
    Tcb_a.resize(nReactions,0.0);
    Tcb_b.resize(nReactions,0.0);

    /*--- Assign gas properties ---*/
    // Rotational modes of energy storage
    RotationModes[0] = 0.0; // e-
    RotationModes[1] = 2.0; // N2
    RotationModes[2] = 2.0; // O2
    RotationModes[3] = 2.0; // NO
    RotationModes[4] = 0.0; // N
    RotationModes[5] = 0.0; // O
    RotationModes[6] = 2.0; // NO+
    RotationModes[7] = 2.0; // N2+
    RotationModes[8] = 2.0; // O2+
    RotationModes[9] = 0.0; // N+
    RotationModes[10]= 0.0; // O+

    // Molar mass [kg/kmol]
    MolarMass[0] = 5.4858E-04;      // e-
    MolarMass[1] = 2.0*14.0067;     // N2
    MolarMass[2] = 2.0*15.9994;     // O2
    MolarMass[3] = 14.0067+15.9994; // NO
    MolarMass[4] = 14.0067;         // N
    MolarMass[5] = 15.9994;         // O
    MolarMass[6] = 14.0067+15.9994 - 5.4858E-04; // NO+
    MolarMass[7] = 2.0*14.0067     - 5.4858E-04; // N2+
    MolarMass[8] = 2.0*15.9994     - 5.4858E-04; // O2+
    MolarMass[9] = 14.0067         - 5.4858E-04; // N+
    MolarMass[10]= 15.9994         - 5.4858E-04; // O+

    //Characteristic vibrational temperatures
    CharVibTemp[0] = 0.0;    // e-
    CharVibTemp[1] = 3408.464; // N2
    CharVibTemp[2] = 2276.979; // O2
    CharVibTemp[3] = 2759.293; // NO
    CharVibTemp[4] = 0.0;    // N
    CharVibTemp[5] = 0.0;    // O
    CharVibTemp[6] = 3473.491; // NO+
    CharVibTemp[7] = 3253.157; // N2+
    CharVibTemp[8] = 2887.139; // O2+
    CharVibTemp[9] = 0.0;    // N+
    CharVibTemp[10]= 0.0;    // O+

    // Formation enthalpy: (Scalabrin values, J/kg)
    Enthalpy_Formation[0] = 0.0;    // e-
    Enthalpy_Formation[1] = 0.0;    // N2
    Enthalpy_Formation[2] = 0.0;    // O2
    Enthalpy_Formation[3] = 3.0357E6;  // NO
    Enthalpy_Formation[4] = 3.373E7; // N
    Enthalpy_Formation[5] = 1.5577E7; // O
    Enthalpy_Formation[6] = 3.3016E7; // NO+
    Enthalpy_Formation[7] = 5.3886E7;  // N2+
    Enthalpy_Formation[8] = 3.6609E7;  // O2+
    Enthalpy_Formation[9] = 1.3436E8;  // N+
    Enthalpy_Formation[10]= 9.8059E7;  // O+

    // Reference temperature (JANAF values, [K])
    Ref_Temperature[0] = 0.0;
    Ref_Temperature[1] = 0.0;
    Ref_Temperature[2] = 0.0;
    Ref_Temperature[3] = 0.0;
    Ref_Temperature[4] = 0.0;
    Ref_Temperature[5] = 0.0;
    Ref_Temperature[6] = 0.0;
    Ref_Temperature[7] = 0.0;
    Ref_Temperature[8] = 0.0;
    Ref_Temperature[9] = 0.0;
    Ref_Temperature[10]= 0.0;

    // Blottner viscosity coefficients
    // A                        // B                        // C
    Blottner(0,0) = 0.00E+0;   Blottner(0,1) =  0.00E+0;  Blottner(0,2) = -1.20E1;  // e-
    Blottner(1,0) = 2.68E-2;   Blottner(1,1) =  3.18E-1;  Blottner(1,2) = -1.13E1;  // N2
    Blottner(2,0) = 4.49E-2;   Blottner(2,1) = -8.26E-2;  Blottner(2,2) = -9.20E0;  // O2
    Blottner(3,0) = 4.36E-2;   Blottner(3,1) = -3.36E-2;  Blottner(3,2) = -9.58E0;  // NO
    Blottner(4,0) = 1.16E-2;   Blottner(4,1) =  6.03E-1;  Blottner(4,2) = -1.24E1;  // N
    Blottner(5,0) = 2.03E-2;   Blottner(5,1) =  4.29E-1;  Blottner(5,2) = -1.16E1;  // O
    Blottner(6,0) = 3.02E-1;   Blottner(6,1) =  -3.50E0;  Blottner(6,2) = -3.74E0;  // NO+
    //Check the following Blottner values are just copied, not used in the code
    Blottner(7,0)  = 2.68E-2;  Blottner(7,1)  =  3.18E-1;  Blottner(7,2)  = -1.13E1; // N2+
    Blottner(8,0)  = 4.49E-2;  Blottner(8,1)  = -8.26E-2;  Blottner(8,2)  = -9.20E0; // O2+
    Blottner(9,0)  = 1.16E-2;  Blottner(9,1)  =  6.03E-1;  Blottner(9,2)  = -1.24E1; // N+
    Blottner(10,0) = 2.03E-2;  Blottner(10,1) =  4.29E-1;  Blottner(10,2) = -1.16E1; // O+

    // Number of electron states
    nElStates[0] = 1;  // e-
    nElStates[1] = 15; // N2
    nElStates[2] = 7;  // O2
    nElStates[3] = 16; // NO
    nElStates[4] = 3;  // N
    nElStates[5] = 5;  // O
    nElStates[6] = 8;  // NO+
    nElStates[7] = 9;  // N2+
    nElStates[8] = 5;  // O2+
    nElStates[9] = 9;  // N+
    nElStates[10]= 5;  // O+

    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      maxEl = max(maxEl, nElStates[iSpecies]);

    /*--- Allocate and initialize electron data arrays ---*/
    CharElTemp.resize(nSpecies,maxEl) = su2double(0.0);
    ElDegeneracy.resize(nSpecies,maxEl) = su2double(0.0);

    // e: 1 state
    CharElTemp(0,0) = 0.000000000000000E+00;
    ElDegeneracy(0,0) = 1;

    //N2: 15 states
    CharElTemp(1,0)  = 0.000000000000000E+00;
    CharElTemp(1,1)  = 7.223156514095200E+04;
    CharElTemp(1,2)  = 8.577862640384000E+04;
    CharElTemp(1,3)  = 8.605026716160000E+04;
    CharElTemp(1,4)  = 9.535118627874400E+04;
    CharElTemp(1,5)  = 9.805635702203200E+04;
    CharElTemp(1,6)  = 9.968267656935200E+04;
    CharElTemp(1,7)  = 1.048976467715200E+05;
    CharElTemp(1,8)  = 1.116489555200000E+05;
    CharElTemp(1,9)  = 1.225836470400000E+05;
    CharElTemp(1,10) = 1.248856873600000E+05;
    CharElTemp(1,11) = 1.282476158188320E+05;
    CharElTemp(1,12) = 1.338060936000000E+05;
    CharElTemp(1,13) = 1.404296391107200E+05;
    CharElTemp(1,14) = 1.504958859200000E+05;
    ElDegeneracy(1,0)  = 1;
    ElDegeneracy(1,1)  = 3;
    ElDegeneracy(1,2)  = 6;
    ElDegeneracy(1,3)  = 6;
    ElDegeneracy(1,4)  = 3;
    ElDegeneracy(1,5)  = 1;
    ElDegeneracy(1,6)  = 2;
    ElDegeneracy(1,7)  = 2;
    ElDegeneracy(1,8)  = 5;
    ElDegeneracy(1,9)  = 1;
    ElDegeneracy(1,10) = 6;
    ElDegeneracy(1,11) = 6;
    ElDegeneracy(1,12) = 10;
    ElDegeneracy(1,13) = 6;
    ElDegeneracy(1,14) = 6;
    // O2: 7 states
    CharElTemp(2,0) = 0.000000000000000E+00;
    CharElTemp(2,1) = 1.139156019700800E+04;
    CharElTemp(2,2) = 1.898473947826400E+04;
    CharElTemp(2,3) = 4.755973576639200E+04;
    CharElTemp(2,4) = 4.991242097343200E+04;
    CharElTemp(2,5) = 5.092268575561600E+04;
    CharElTemp(2,6) = 7.189863255967200E+04;
    ElDegeneracy(2,0) = 3;
    ElDegeneracy(2,1) = 2;
    ElDegeneracy(2,2) = 1;
    ElDegeneracy(2,3) = 1;
    ElDegeneracy(2,4) = 6;
    ElDegeneracy(2,5) = 3;
    ElDegeneracy(2,6) = 3;
    // NO: 16 states
    CharElTemp(3,0)  = 0.000000000000000E+00;
    CharElTemp(3,1)  = 5.467345760000000E+04;
    CharElTemp(3,2)  = 6.317139627802400E+04;
    CharElTemp(3,3)  = 6.599450342445600E+04;
    CharElTemp(3,4)  = 6.906120960000000E+04;
    CharElTemp(3,5)  = 7.049998480000000E+04;
    CharElTemp(3,6)  = 7.491055017560000E+04;
    CharElTemp(3,7)  = 7.628875293968000E+04;
    CharElTemp(3,8)  = 8.676188537552000E+04;
    CharElTemp(3,9)  = 8.714431182368000E+04;
    CharElTemp(3,10) = 8.886077063728000E+04;
    CharElTemp(3,11) = 8.981755614528000E+04;
    CharElTemp(3,12) = 8.988445919208000E+04;
    CharElTemp(3,13) = 9.042702132000000E+04;
    CharElTemp(3,14) = 9.064283760000000E+04;
    CharElTemp(3,15) = 9.111763341600000E+04;
    ElDegeneracy(3,0)  = 4;
    ElDegeneracy(3,1)  = 8;
    ElDegeneracy(3,2)  = 2;
    ElDegeneracy(3,3)  = 4;
    ElDegeneracy(3,4)  = 4;
    ElDegeneracy(3,5)  = 4;
    ElDegeneracy(3,6)  = 4;
    ElDegeneracy(3,7)  = 2;
    ElDegeneracy(3,8)  = 4;
    ElDegeneracy(3,9)  = 2;
    ElDegeneracy(3,10) = 4;
    ElDegeneracy(3,11) = 4;
    ElDegeneracy(3,12) = 2;
    ElDegeneracy(3,13) = 2;
    ElDegeneracy(3,14) = 2;
    ElDegeneracy(3,15) = 4;
    // N: 3 states
    CharElTemp(4,0) = 0.000000000000000E+00;
    CharElTemp(4,1) = 2.766469645581980E+04;
    CharElTemp(4,2) = 4.149309313560210E+04;
    ElDegeneracy(4,0)= 4;
    ElDegeneracy(4,1)= 10;
    ElDegeneracy(4,2)= 6;
    // O: 5 states
    CharElTemp(5,0) = 0.000000000000000E+00;
    CharElTemp(5,1) = 2.277077570280000E+02;
    CharElTemp(5,2) = 3.265688785704000E+02;
    CharElTemp(5,3) = 2.283028632262240E+04;
    CharElTemp(5,4) = 4.861993036434160E+04;
    ElDegeneracy(5,0) = 5;
    ElDegeneracy(5,1) = 3;
    ElDegeneracy(5,2) = 1;
    ElDegeneracy(5,3) = 5;
    ElDegeneracy(5,4) = 1;
    // NO+: 8 states
    CharElTemp(6,0) = 0.000000000000000E+00;
    CharElTemp(6,1) = 7.508967768800000E+04;
    CharElTemp(6,2) = 8.525462447600000E+04;
    CharElTemp(6,3) = 8.903572570160000E+04;
    CharElTemp(6,4) = 9.746982592400000E+04;
    CharElTemp(6,5) = 1.000553049584000E+05;
    CharElTemp(6,6) = 1.028033655904000E+05;
    CharElTemp(6,7) = 1.057138639424800E+05;
    ElDegeneracy(6,0) = 1;
    ElDegeneracy(6,1) = 3;
    ElDegeneracy(6,2) = 6;
    ElDegeneracy(6,3) = 6;
    ElDegeneracy(6,4) = 3;
    ElDegeneracy(6,5) = 1;
    ElDegeneracy(6,6) = 2;
    ElDegeneracy(6,7) = 2;

    // N2+: 9 states
    CharElTemp(7,0) = 0.000000000000000E+00;
    CharElTemp(7,1) = 1.318926784230000E+04;
    CharElTemp(7,2) = 3.663269865090000E+04;
    CharElTemp(7,3) = 3.668881095000000E+04;
    CharElTemp(7,4) = 5.985311904000000E+04;
    CharElTemp(7,5) = 6.618373740000000E+04;
    CharElTemp(7,6) = 7.598900197350000E+04;
    CharElTemp(7,7) = 7.625517570000000E+04;
    CharElTemp(7,8) = 8.201028330000000E+04;
    ElDegeneracy(7,0) = 2;
    ElDegeneracy(7,1) = 4;
    ElDegeneracy(7,2) = 2;
    ElDegeneracy(7,3) = 4;
    ElDegeneracy(7,4) = 8;
    ElDegeneracy(7,5) = 8;
    ElDegeneracy(7,6) = 4;
    ElDegeneracy(7,7) = 4;
    ElDegeneracy(7,8) = 4;

    // O2+: 5 states
    CharElTemp(8,0) = 0.000000000000000E+00;
    CharElTemp(8,1) = 4.735446410970000E+04;
    CharElTemp(8,2) = 5.837405638680000E+04;
    CharElTemp(8,3) = 5.841434214000000E+04;
    CharElTemp(8,4) = 6.229903977000000E+04;
    ElDegeneracy(8,0) = 4;
    ElDegeneracy(8,1) = 8;
    ElDegeneracy(8,2) = 4;
    ElDegeneracy(8,3) = 4;
    ElDegeneracy(8,4) = 6;

    // N+: 9 states
    CharElTemp(9,0) = 0.000000000000000E+00;
    CharElTemp(9,1) = 7.006843503000000E+01;
    CharElTemp(9,2) = 1.881920185200000E+02;
    CharElTemp(9,3) = 2.203659475578000E+04;
    CharElTemp(9,4) = 4.703189032872000E+04;
    CharElTemp(9,5) = 6.731260175574000E+04;
    CharElTemp(9,6) = 1.327087526806800E+05;
    CharElTemp(9,7) = 1.327276006580700E+05;
    CharElTemp(9,8) = 1.327297588234200E+05;
    ElDegeneracy(9,0) = 1;
    ElDegeneracy(9,1) = 3;
    ElDegeneracy(9,2) = 5;
    ElDegeneracy(9,3) = 5;
    ElDegeneracy(9,4) = 1;
    ElDegeneracy(9,5) = 5;
    ElDegeneracy(9,6) = 7;
    ElDegeneracy(9,7) = 5;
    ElDegeneracy(9,8) = 3;

    // O+: 5 states
    CharElTemp(10,0) = 0.000000000000000E+00;
    CharElTemp(10,1) = 3.857130664596000E+04;
    CharElTemp(10,2) = 3.860152096086000E+04;
    CharElTemp(10,3) = 5.822284093461000E+04;
    CharElTemp(10,4) = 5.822499909996000E+04;
    ElDegeneracy(10,0) = 4;
    ElDegeneracy(10,1) = 6;
    ElDegeneracy(10,2) = 4;
    ElDegeneracy(10,3) = 4;
    ElDegeneracy(10,4) = 2;

    /*--- Set reaction maps ---*/
    // N2 dissociation    // N2 + M -> 2N + M  (M = N2, O2, NO, N, O, NO+, N2+, O2+, N+, O+)
    Reactions(0,0,0)=1;   Reactions(0,0,1)=1;   Reactions(0,0,2)=nSpecies;   Reactions(0,1,0)=4;   Reactions(0,1,1)=4;   Reactions(0,1,2)=1;
    Reactions(1,0,0)=1;   Reactions(1,0,1)=2;   Reactions(1,0,2)=nSpecies;   Reactions(1,1,0)=4;   Reactions(1,1,1)=4;   Reactions(1,1,2)=2;
    Reactions(2,0,0)=1;   Reactions(2,0,1)=3;   Reactions(2,0,2)=nSpecies;   Reactions(2,1,0)=4;   Reactions(2,1,1)=4;   Reactions(2,1,2)=3;
    Reactions(3,0,0)=1;   Reactions(3,0,1)=4;   Reactions(3,0,2)=nSpecies;   Reactions(3,1,0)=4;   Reactions(3,1,1)=4;   Reactions(3,1,2)=4;
    Reactions(4,0,0)=1;   Reactions(4,0,1)=5;   Reactions(4,0,2)=nSpecies;   Reactions(4,1,0)=4;   Reactions(4,1,1)=4;   Reactions(4,1,2)=5;
    Reactions(5,0,0)=1;   Reactions(5,0,1)=6;   Reactions(5,0,2)=nSpecies;   Reactions(5,1,0)=4;   Reactions(5,1,1)=4;   Reactions(5,1,2)=6;
    Reactions(6,0,0)=1;   Reactions(6,0,1)=7;   Reactions(6,0,2)=nSpecies;   Reactions(6,1,0)=4;   Reactions(6,1,1)=4;   Reactions(6,1,2)=7;
    Reactions(7,0,0)=1;   Reactions(7,0,1)=8;   Reactions(7,0,2)=nSpecies;   Reactions(7,1,0)=4;   Reactions(7,1,1)=4;   Reactions(7,1,2)=8;
    Reactions(8,0,0)=1;   Reactions(8,0,1)=9;   Reactions(8,0,2)=nSpecies;   Reactions(8,1,0)=4;   Reactions(8,1,1)=4;   Reactions(8,1,2)=9;
    Reactions(9,0,0)=1;   Reactions(9,0,1)=10;  Reactions(9,0,2)=nSpecies;   Reactions(9,1,0)=4;   Reactions(9,1,1)=4;   Reactions(9,1,2)=10;

    // O2 dissociation    // O2 + M -> 2O + M  (M = N2, O2, NO, N, O, NO+, N2+, O2+, N+, O+)
    Reactions(10,0,0)=2;   Reactions(10,0,1)=1;   Reactions(10,0,2)=nSpecies;   Reactions(10,1,0)=5;   Reactions(10,1,1)=5;   Reactions(10,1,2)=1;
    Reactions(11,0,0)=2;   Reactions(11,0,1)=2;   Reactions(11,0,2)=nSpecies;   Reactions(11,1,0)=5;   Reactions(11,1,1)=5;   Reactions(11,1,2)=2;
    Reactions(12,0,0)=2;   Reactions(12,0,1)=3;   Reactions(12,0,2)=nSpecies;   Reactions(12,1,0)=5;   Reactions(12,1,1)=5;   Reactions(12,1,2)=3;
    Reactions(13,0,0)=2;   Reactions(13,0,1)=4;   Reactions(13,0,2)=nSpecies;   Reactions(13,1,0)=5;   Reactions(13,1,1)=5;   Reactions(13,1,2)=4;
    Reactions(14,0,0)=2;   Reactions(14,0,1)=5;   Reactions(14,0,2)=nSpecies;   Reactions(14,1,0)=5;   Reactions(14,1,1)=5;   Reactions(14,1,2)=5;
    Reactions(15,0,0)=2;   Reactions(15,0,1)=6;   Reactions(15,0,2)=nSpecies;   Reactions(15,1,0)=5;   Reactions(15,1,1)=5;   Reactions(15,1,2)=6;
    Reactions(16,0,0)=2;   Reactions(16,0,1)=7;   Reactions(16,0,2)=nSpecies;   Reactions(16,1,0)=5;   Reactions(16,1,1)=5;   Reactions(16,1,2)=7;
    Reactions(17,0,0)=2;   Reactions(17,0,1)=8;   Reactions(17,0,2)=nSpecies;   Reactions(17,1,0)=5;   Reactions(17,1,1)=5;   Reactions(17,1,2)=8;
    Reactions(18,0,0)=2;   Reactions(18,0,1)=9;   Reactions(18,0,2)=nSpecies;   Reactions(18,1,0)=5;   Reactions(18,1,1)=5;   Reactions(18,1,2)=9;
    Reactions(19,0,0)=2;   Reactions(19,0,1)=10;  Reactions(19,0,2)=nSpecies;   Reactions(19,1,0)=5;   Reactions(19,1,1)=5;   Reactions(19,1,2)=10;

    // NO dissociation    // NO + M -> N + O + M  (M = N2, O2, NO, N, O, NO+, N2+, O2+, N+, O+)
    Reactions(20,0,0)=3;   Reactions(20,0,1)=1;   Reactions(20,0,2)=nSpecies;   Reactions(20,1,0)=4;   Reactions(20,1,1)=5;   Reactions(20,1,2)=1;
    Reactions(21,0,0)=3;   Reactions(21,0,1)=2;   Reactions(21,0,2)=nSpecies;   Reactions(21,1,0)=4;   Reactions(21,1,1)=5;   Reactions(21,1,2)=2;
    Reactions(22,0,0)=3;   Reactions(22,0,1)=3;   Reactions(22,0,2)=nSpecies;   Reactions(22,1,0)=4;   Reactions(22,1,1)=5;   Reactions(22,1,2)=3;
    Reactions(23,0,0)=3;   Reactions(23,0,1)=4;   Reactions(23,0,2)=nSpecies;   Reactions(23,1,0)=4;   Reactions(23,1,1)=5;   Reactions(23,1,2)=4;
    Reactions(24,0,0)=3;   Reactions(24,0,1)=5;   Reactions(24,0,2)=nSpecies;   Reactions(24,1,0)=4;   Reactions(24,1,1)=5;   Reactions(24,1,2)=5;
    Reactions(25,0,0)=3;   Reactions(25,0,1)=6;   Reactions(25,0,2)=nSpecies;   Reactions(25,1,0)=4;   Reactions(25,1,1)=5;   Reactions(25,1,2)=6;
    Reactions(26,0,0)=3;   Reactions(26,0,1)=7;   Reactions(26,0,2)=nSpecies;   Reactions(26,1,0)=4;   Reactions(26,1,1)=5;   Reactions(26,1,2)=7;
    Reactions(27,0,0)=3;   Reactions(27,0,1)=8;   Reactions(27,0,2)=nSpecies;   Reactions(27,1,0)=4;   Reactions(27,1,1)=5;   Reactions(27,1,2)=8;
    Reactions(28,0,0)=3;   Reactions(28,0,1)=9;   Reactions(28,0,2)=nSpecies;   Reactions(28,1,0)=4;   Reactions(28,1,1)=5;   Reactions(28,1,2)=9;
    Reactions(29,0,0)=3;   Reactions(29,0,1)=10;  Reactions(29,0,2)=nSpecies;   Reactions(29,1,0)=4;   Reactions(29,1,1)=5;   Reactions(29,1,2)=10;

    // N2 + O -> NO + N
    Reactions(30,0,0)=1;   Reactions(30,0,1)=5;   Reactions(30,0,2)=nSpecies;   Reactions(30,1,0)=3;   Reactions(30,1,1)=4;   Reactions(30,1,2)=nSpecies;
    // NO + O -> O2 + N
    Reactions(31,0,0)=3;   Reactions(31,0,1)=5;   Reactions(31,0,2)=nSpecies;   Reactions(31,1,0)=2;   Reactions(31,1,1)=4;   Reactions(31,1,2)=nSpecies;
    //N + O -> NO+ + e
    Reactions(32,0,0)=4;   Reactions(32,0,1)=5;   Reactions(32,0,2)=nSpecies;   Reactions(32,1,0)=6;   Reactions(32,1,1)=0;   Reactions(32,1,2)=nSpecies;
    //N2 + e -> N + N + e
    Reactions(33,0,0)=1;   Reactions(33,0,1)=0;   Reactions(33,0,2)=nSpecies;   Reactions(33,1,0)=4;   Reactions(33,1,1)=4;   Reactions(33,1,2)=0;

    // O + O -> O2+ + e-
    Reactions(34,0,0)=5;   Reactions(34,0,1)=5;   Reactions(34,0,2)=nSpecies;   Reactions(34,1,0)=8;   Reactions(34,1,1)=0;   Reactions(34,1,2)=nSpecies;

    // N + N -> N2+ + e-
    Reactions(35,0,0)=4;   Reactions(35,0,1)=4;   Reactions(35,0,2)=nSpecies;   Reactions(35,1,0)=7;   Reactions(35,1,1)=0;   Reactions(35,1,2)=nSpecies;

    // NO+ + O -> N+ + O2
    Reactions(36,0,0)=6;   Reactions(36,0,1)=5;   Reactions(36,0,2)=nSpecies;   Reactions(36,1,0)=9;   Reactions(36,1,1)=2;   Reactions(36,1,2)=nSpecies;

    // N+ + N2 -> N2+ + N
    //Reactions(37,0,0)=9;   Reactions(37,0,1)=1;   Reactions(37,0,2)=nSpecies;   Reactions(37,1,0)=7;   Reactions(37,1,1)=4;   Reactions(37,1,2)=nSpecies;

    // O2+ + N -> N+ + O2
    Reactions(37,0,0)=8;   Reactions(37,0,1)=4;   Reactions(37,0,2)=nSpecies;   Reactions(37,1,0)=9;   Reactions(37,1,1)=2;   Reactions(37,1,2)=nSpecies;

    // O+ + NO -> N+ + O2
    Reactions(38,0,0)=10;  Reactions(38,0,1)=3;   Reactions(38,0,2)=nSpecies;   Reactions(38,1,0)=9;   Reactions(38,1,1)=2;   Reactions(38,1,2)=nSpecies;

    // O2+ + N2 -> N2+ + O2
    Reactions(39,0,0)=8;   Reactions(39,0,1)=1;   Reactions(39,0,2)=nSpecies;   Reactions(39,1,0)=7;   Reactions(39,1,1)=2;   Reactions(39,1,2)=nSpecies;

    // O2+ + O -> O+ + O2
    //Reactions(41,0,0)=8;   Reactions(41,0,1)=5;   Reactions(41,0,2)=nSpecies;   Reactions(41,1,0)=10;  Reactions(41,1,1)=2;   Reactions(41,1,2)=nSpecies;

    // NO+ + N -> O+ + N2
    Reactions(40,0,0)=6;   Reactions(40,0,1)=4;   Reactions(40,0,2)=nSpecies;   Reactions(40,1,0)=10;  Reactions(40,1,1)=1;   Reactions(40,1,2)=nSpecies;

    // NO+ + O2 -> O2+ + NO
    Reactions(41,0,0)=6;   Reactions(41,0,1)=2;   Reactions(41,0,2)=nSpecies;   Reactions(41,1,0)=8;   Reactions(41,1,1)=3;   Reactions(41,1,2)=nSpecies;

    // NO+ + O -> O2+ + N
    Reactions(42,0,0)=6;   Reactions(42,0,1)=5;   Reactions(42,0,2)=nSpecies;   Reactions(42,1,0)=8;   Reactions(42,1,1)=4;   Reactions(42,1,2)=nSpecies;

    // O+ + N2 -> N2+ + O
    Reactions(43,0,0)=10;  Reactions(43,0,1)=1;   Reactions(43,0,2)=nSpecies;   Reactions(43,1,0)=7;   Reactions(43,1,1)=5;   Reactions(43,1,2)=nSpecies;

    // NO+ + N -> N2+ + O
    Reactions(44,0,0)=6;   Reactions(44,0,1)=4;   Reactions(44,0,2)=nSpecies;   Reactions(44,1,0)=7;   Reactions(44,1,1)=5;   Reactions(44,1,2)=nSpecies;

    // O + e- -> O+ + e- + e-
    Reactions(45,0,0)=5;   Reactions(45,0,1)=0;   Reactions(45,0,2)=nSpecies;   Reactions(45,1,0)=10;  Reactions(45,1,1)=0;   Reactions(45,1,2)=0;

    // N + e- -> N+ + e- + e-
    Reactions(46,0,0)=4;   Reactions(46,0,1)=0;   Reactions(46,0,2)=nSpecies;   Reactions(46,1,0)=9;   Reactions(46,1,1)=0;   Reactions(46,1,2)=0;

    /*--- Set Arrhenius coefficients for reactions ---*/
    // Pre-exponential factor
    ArrheniusCoefficient[0]  = 7.0E+21;   // M = N2,  eff = 1.0
    ArrheniusCoefficient[1]  = 7.0E+21;   // M = O2,  eff = 1.0
    ArrheniusCoefficient[2]  = 7.0E+21;   // M = NO,  eff = 1.0
    ArrheniusCoefficient[3]  = 3.0E+22;   // M = N,   eff = 4.2857
    ArrheniusCoefficient[4]  = 3.0E+22;   // M = O,   eff = 4.2857
    ArrheniusCoefficient[5]  = 7.0E+21;   // M = NO+, eff = 1.0
    ArrheniusCoefficient[6]  = 7.0E+21;   // M = N2+, eff = 1.0
    ArrheniusCoefficient[7]  = 7.0E+21;   // M = O2+, eff = 1.0
    ArrheniusCoefficient[8]  = 7.0E+21;   // M = N+,  eff = 1.0
    ArrheniusCoefficient[9]  = 7.0E+21;   // M = O+,  eff = 1.0

    ArrheniusCoefficient[10] = 2.0E+21;   // M = N2,  eff = 1.0
    ArrheniusCoefficient[11] = 2.0E+21;   // M = O2,  eff = 1.0
    ArrheniusCoefficient[12] = 2.0E+21;   // M = NO,  eff = 1.0
    ArrheniusCoefficient[13] = 1.0E+22;   // M = N,   eff = 5.0
    ArrheniusCoefficient[14] = 1.0E+22;   // M = O,   eff = 5.0
    ArrheniusCoefficient[15] = 2.0E+21;   // M = NO+, eff = 1.0
    ArrheniusCoefficient[16] = 2.0E+21;   // M = N2+, eff = 1.0
    ArrheniusCoefficient[17] = 2.0E+21;   // M = O2+, eff = 1.0
    ArrheniusCoefficient[18] = 2.0E+21;   // M = N+,  eff = 1.0
    ArrheniusCoefficient[19] = 2.0E+21;   // M = O+,  eff = 1.0

    ArrheniusCoefficient[20] = 5.0E+15;   // M = N2,  eff = 1.0
    ArrheniusCoefficient[21] = 5.0E+15;   // M = O2,  eff = 1.0
    ArrheniusCoefficient[22] = 1.1E+17;   // M = NO,  eff = 22.0
    ArrheniusCoefficient[23] = 1.1E+17;   // M = N,   eff = 22.0
    ArrheniusCoefficient[24] = 1.1E+17;   // M = O,   eff = 22.0
    ArrheniusCoefficient[25] = 5.0E+15;   // M = NO+, eff = 1.0
    ArrheniusCoefficient[26] = 5.0E+15;   // M = N2+, eff = 1.0
    ArrheniusCoefficient[27] = 5.0E+15;   // M = O2+, eff = 1.0
    ArrheniusCoefficient[28] = 1.1E+17;   // M = N+,  eff = 22.0
    ArrheniusCoefficient[29] = 1.1E+17;   // M = O+,  eff = 22.0

    // Rate-controlling temperature exponent
    ArrheniusEta[0]  = -1.6;   ArrheniusTheta[0]  = 113200.0;
    ArrheniusEta[1]  = -1.6;   ArrheniusTheta[1]  = 113200.0;
    ArrheniusEta[2]  = -1.6;   ArrheniusTheta[2]  = 113200.0;
    ArrheniusEta[3]  = -1.6;   ArrheniusTheta[3]  = 113200.0;
    ArrheniusEta[4]  = -1.6;   ArrheniusTheta[4]  = 113200.0;
    ArrheniusEta[5]  = -1.6;   ArrheniusTheta[5]  = 113200.0;
    ArrheniusEta[6]  = -1.6;   ArrheniusTheta[6]  = 113200.0;
    ArrheniusEta[7]  = -1.6;   ArrheniusTheta[7]  = 113200.0;
    ArrheniusEta[8]  = -1.6;   ArrheniusTheta[8]  = 113200.0;
    ArrheniusEta[9]  = -1.6;   ArrheniusTheta[9]  = 113200.0;

    ArrheniusEta[10] = -1.5;   ArrheniusTheta[10] = 59360.0;
    ArrheniusEta[11] = -1.5;   ArrheniusTheta[11] = 59360.0;
    ArrheniusEta[12] = -1.5;   ArrheniusTheta[12] = 59360.0;
    ArrheniusEta[13] = -1.5;   ArrheniusTheta[13] = 59360.0;
    ArrheniusEta[14] = -1.5;   ArrheniusTheta[14] = 59360.0;
    ArrheniusEta[15] = -1.5;   ArrheniusTheta[15] = 59360.0;
    ArrheniusEta[16] = -1.5;   ArrheniusTheta[16] = 59360.0;
    ArrheniusEta[17] = -1.5;   ArrheniusTheta[17] = 59360.0;
    ArrheniusEta[18] = -1.5;   ArrheniusTheta[18] = 59360.0;
    ArrheniusEta[19] = -1.5;   ArrheniusTheta[19] = 59360.0;

    ArrheniusEta[20] = 0.0;   ArrheniusTheta[20] = 75500.0;
    ArrheniusEta[21] = 0.0;   ArrheniusTheta[21] = 75500.0;
    ArrheniusEta[22] = 0.0;   ArrheniusTheta[22] = 75500.0;
    ArrheniusEta[23] = 0.0;   ArrheniusTheta[23] = 75500.0;
    ArrheniusEta[24] = 0.0;   ArrheniusTheta[24] = 75500.0;
    ArrheniusEta[25] = 0.0;   ArrheniusTheta[25] = 75500.0;
    ArrheniusEta[26] = 0.0;   ArrheniusTheta[26] = 75500.0;
    ArrheniusEta[27] = 0.0;   ArrheniusTheta[27] = 75500.0;
    ArrheniusEta[28] = 0.0;   ArrheniusTheta[28] = 75500.0;
    ArrheniusEta[29] = 0.0;   ArrheniusTheta[29] = 75500.0;

    // Reaction 30: N2 + O -> NO + N
    ArrheniusCoefficient[30] = 5.7E+12;
    ArrheniusEta[30]         = 0.42;
    ArrheniusTheta[30]       = 42938.0;

    // Reaction 31: NO + O -> O2 + N
    ArrheniusCoefficient[31] = 8.4E+12;
    ArrheniusEta[31]         = 0.0;
    ArrheniusTheta[31]       = 19400.0;

    // Reaction 32: N + O -> NO+ + e-
    ArrheniusCoefficient[32] = 5.3E+12;
    ArrheniusEta[32]         = 0.0;
    ArrheniusTheta[32]       = 31900.0;

    // Reaction 33: N2 + e- -> 2N + e-
    ArrheniusCoefficient[33] = 3.0E+24;
    ArrheniusEta[33]         = -1.6;
    ArrheniusTheta[33]       = 113200.0;

    // Reaction 34: O + O -> O2+ + e-
    ArrheniusCoefficient[34] = 7.1E+02;
    ArrheniusEta[34]         = 2.7;
    ArrheniusTheta[34]       = 80600.0;

    // Reaction 35: N + N -> N2+ + e-
    ArrheniusCoefficient[35] = 4.4E+07;
    ArrheniusEta[35]         = 1.5;
    ArrheniusTheta[35]       = 67500.0;

    // Reaction 36: NO+ + O -> N+ + O2
    ArrheniusCoefficient[36] = 1.0E+12;
    ArrheniusEta[36]         = 0.5;
    ArrheniusTheta[36]       = 77200.0;

    // Reaction 37: N+ + N2 -> N2+ + N
    //ArrheniusCoefficient[37] = 1.0E+12;
    //ArrheniusEta[37]         = 0.5;
    //ArrheniusTheta[37]       = 12200.0;

    // Reaction 37: O2+ + N -> N+ + O2
    ArrheniusCoefficient[37] = 8.7E+13;
    ArrheniusEta[37]         = 0.14;
    ArrheniusTheta[37]       = 28600.0;

    // Reaction 38: O+ + NO -> N+ + O2
    ArrheniusCoefficient[38] = 1.4E+05;
    ArrheniusEta[38]         = 1.9;
    ArrheniusTheta[38]       = 26600.0;

    // Reaction 39: O2+ + N2 -> N2+ + O2
    ArrheniusCoefficient[39] = 9.9E+12;
    ArrheniusEta[39]         = 0.0;
    ArrheniusTheta[39]       = 40700.0;

    // Reaction 41: O2+ + O -> O+ + O2
    //ArrheniusCoefficient[41] = 4.0E+12;
    //ArrheniusEta[41]         = 0.09;
    //ArrheniusTheta[41]       = 18000.0;

    // Reaction 40: NO+ + N -> O+ + N2
    ArrheniusCoefficient[40] = 3.4E+13;
    ArrheniusEta[40]         = -1.08;
    ArrheniusTheta[40]       = 12800.0;

    // Reaction 41: NO+ + O2 -> O2+ + NO
    ArrheniusCoefficient[41] = 2.4E+13;
    ArrheniusEta[41]         = 0.41;
    ArrheniusTheta[41]       = 32600.0;

    // Reaction 42: NO+ + O -> O2+ + N
    ArrheniusCoefficient[42] = 7.2E+12;
    ArrheniusEta[42]         = 0.29;
    ArrheniusTheta[42]       = 48600.0;

    // Reaction 43: O+ + N2 -> N2+ + O
    ArrheniusCoefficient[43] = 9.1E+11;
    ArrheniusEta[43]         = 0.36;
    ArrheniusTheta[43]       = 22800.0;

    // Reaction 44: NO+ + N -> N2+ + O
    ArrheniusCoefficient[44] = 7.2E+13;
    ArrheniusEta[44]         = 0.0;
    ArrheniusTheta[44]       = 35500.0;

    // Reaction 45: O + e- -> O+ + e- + e-
    ArrheniusCoefficient[45] = 3.9E+33;
    ArrheniusEta[45]         = -3.78;
    ArrheniusTheta[45]       = 158500.0;

    // Reaction 46: N + e- -> N+ + e- + e-
    ArrheniusCoefficient[46] = 2.5E+34;
    ArrheniusEta[46]         = -3.82;
    ArrheniusTheta[46]       = 168200.0;
    /*--- Set rate-controlling temperature exponents ---*/
    //  -----------  Tc = Ttr^a * Tve^b  -----------
    //
    // Forward Reactions
    //   Dissociation:         a = 0.5, b = 0.5  (OR a = 0.7, b =0.3)
    //   Exchange:             a = 1,   b = 0
    //   Associative ion...    a = 1,   b = 0  ???
    //   E Impact dissociation a = 0,   b = 1
    //   E Impact ionization:  a = 0,   b = 1
    //
    // Backward Reactions
    //   Dissociation:           a = 1,   b = 0
    //   Exchange:               a = 1,   b = 0
    //   Associative  ion...     a = 0.5, b = 0.5
    //   E Impact ionization:    a = 0,   b = 1
    //   E Impact dissocitation: a = 0.5, b = 0.5 ???
    //   N2 impact dissociation: a = 0,   b = 1
    //   Others:                 a = 1,   b = 0

    // Reactions 0-9: N2+M dissociation
     Tcf_a[0]  = 0.5; Tcf_b[0]  = 0.5; Tcb_a[0]  = 1.0; Tcb_b[0]  = 0.0;
     Tcf_a[1]  = 0.5; Tcf_b[1]  = 0.5; Tcb_a[1]  = 1.0; Tcb_b[1]  = 0.0;
     Tcf_a[2]  = 0.5; Tcf_b[2]  = 0.5; Tcb_a[2]  = 1.0; Tcb_b[2]  = 0.0;
     Tcf_a[3]  = 0.5; Tcf_b[3]  = 0.5; Tcb_a[3]  = 1.0; Tcb_b[3]  = 0.0;
     Tcf_a[4]  = 0.5; Tcf_b[4]  = 0.5; Tcb_a[4]  = 1.0; Tcb_b[4]  = 0.0;
     Tcf_a[5]  = 0.5; Tcf_b[5]  = 0.5; Tcb_a[5]  = 1.0; Tcb_b[5]  = 0.0;
     Tcf_a[6]  = 0.5; Tcf_b[6]  = 0.5; Tcb_a[6]  = 1.0; Tcb_b[6]  = 0.0;
     Tcf_a[7]  = 0.5; Tcf_b[7]  = 0.5; Tcb_a[7]  = 1.0; Tcb_b[7]  = 0.0;
     Tcf_a[8]  = 0.5; Tcf_b[8]  = 0.5; Tcb_a[8]  = 1.0; Tcb_b[8]  = 0.0;
     Tcf_a[9]  = 0.5; Tcf_b[9]  = 0.5; Tcb_a[9]  = 1.0; Tcb_b[9]  = 0.0;

     // Reactions 10-19: O2+M dissociation
     Tcf_a[10] = 0.5; Tcf_b[10] = 0.5; Tcb_a[10] = 1.0; Tcb_b[10] = 0.0;
     Tcf_a[11] = 0.5; Tcf_b[11] = 0.5; Tcb_a[11] = 1.0; Tcb_b[11] = 0.0;
     Tcf_a[12] = 0.5; Tcf_b[12] = 0.5; Tcb_a[12] = 1.0; Tcb_b[12] = 0.0;
     Tcf_a[13] = 0.5; Tcf_b[13] = 0.5; Tcb_a[13] = 1.0; Tcb_b[13] = 0.0;
     Tcf_a[14] = 0.5; Tcf_b[14] = 0.5; Tcb_a[14] = 1.0; Tcb_b[14] = 0.0;
     Tcf_a[15] = 0.5; Tcf_b[15] = 0.5; Tcb_a[15] = 1.0; Tcb_b[15] = 0.0;
     Tcf_a[16] = 0.5; Tcf_b[16] = 0.5; Tcb_a[16] = 1.0; Tcb_b[16] = 0.0;
     Tcf_a[17] = 0.5; Tcf_b[17] = 0.5; Tcb_a[17] = 1.0; Tcb_b[17] = 0.0;
     Tcf_a[18] = 0.5; Tcf_b[18] = 0.5; Tcb_a[18] = 1.0; Tcb_b[18] = 0.0;
     Tcf_a[19] = 0.5; Tcf_b[19] = 0.5; Tcb_a[19] = 1.0; Tcb_b[19] = 0.0;

     // Reactions 20-29: NO+M dissociation
     Tcf_a[20] = 0.5; Tcf_b[20] = 0.5; Tcb_a[20] = 1.0; Tcb_b[20] = 0.0;
     Tcf_a[21] = 0.5; Tcf_b[21] = 0.5; Tcb_a[21] = 1.0; Tcb_b[21] = 0.0;
     Tcf_a[22] = 0.5; Tcf_b[22] = 0.5; Tcb_a[22] = 1.0; Tcb_b[22] = 0.0;
     Tcf_a[23] = 0.5; Tcf_b[23] = 0.5; Tcb_a[23] = 1.0; Tcb_b[23] = 0.0;
     Tcf_a[24] = 0.5; Tcf_b[24] = 0.5; Tcb_a[24] = 1.0; Tcb_b[24] = 0.0;
     Tcf_a[25] = 0.5; Tcf_b[25] = 0.5; Tcb_a[25] = 1.0; Tcb_b[25] = 0.0;
     Tcf_a[26] = 0.5; Tcf_b[26] = 0.5; Tcb_a[26] = 1.0; Tcb_b[26] = 0.0;
     Tcf_a[27] = 0.5; Tcf_b[27] = 0.5; Tcb_a[27] = 1.0; Tcb_b[27] = 0.0;
     Tcf_a[28] = 0.5; Tcf_b[28] = 0.5; Tcb_a[28] = 1.0; Tcb_b[28] = 0.0;
     Tcf_a[29] = 0.5; Tcf_b[29] = 0.5; Tcb_a[29] = 1.0; Tcb_b[29] = 0.0;

    // Reaction 30: N2+O -> NO+N  (exchange)
    Tcf_a[30] = 1.0; Tcf_b[30] = 0.0; Tcb_a[30] = 1.0; Tcb_b[30] = 0.0;

    // Reaction 31: NO+O -> O2+N  (exchange)
    Tcf_a[31] = 1.0; Tcf_b[31] = 0.0; Tcb_a[31] = 1.0; Tcb_b[31] = 0.0;

    // Reaction 32: N+O -> NO++e-  (associative ionization)
    Tcf_a[32] = 1.0; Tcf_b[32] = 0.0; Tcb_a[32] = 0.5; Tcb_b[32] = 0.5;

    // Reaction 33: N2+e- -> 2N+e-  (electron impact dissociation)
    Tcf_a[33] = 0.0; Tcf_b[33] = 1.0; Tcb_a[33] = 0.0; Tcb_b[33] = 1.0;

    // Reaction 34: O+O -> O2++e-  (associative ionization)
    Tcf_a[34] = 1.0; Tcf_b[34] = 0.0; Tcb_a[34] = 0.5; Tcb_b[34] = 0.5;

    // Reaction 35: N+N -> N2++e-  (associative ionization)
    Tcf_a[35] = 1.0; Tcf_b[35] = 0.0; Tcb_a[35] = 0.5; Tcb_b[35] = 0.5;

    // Reactions 36-44: charge exchange
    Tcf_a[36] = 1.0; Tcf_b[36] = 0.0; Tcb_a[36] = 1.0; Tcb_b[36] = 0.0;
    //Tcf_a[37] = 1.0; Tcf_b[37] = 0.0; Tcb_a[37] = 1.0; Tcb_b[37] = 0.0;
    Tcf_a[37] = 1.0; Tcf_b[37] = 0.0; Tcb_a[37] = 1.0; Tcb_b[37] = 0.0;
    Tcf_a[38] = 1.0; Tcf_b[38] = 0.0; Tcb_a[38] = 1.0; Tcb_b[38] = 0.0;
    Tcf_a[39] = 1.0; Tcf_b[39] = 0.0; Tcb_a[39] = 1.0; Tcb_b[39] = 0.0;
    //Tcf_a[41] = 1.0; Tcf_b[41] = 0.0; Tcb_a[41] = 1.0; Tcb_b[41] = 0.0;
    Tcf_a[40] = 1.0; Tcf_b[40] = 0.0; Tcb_a[40] = 1.0; Tcb_b[40] = 0.0;
    Tcf_a[41] = 1.0; Tcf_b[41] = 0.0; Tcb_a[41] = 1.0; Tcb_b[41] = 0.0;
    Tcf_a[42] = 1.0; Tcf_b[42] = 0.0; Tcb_a[42] = 1.0; Tcb_b[42] = 0.0;
    Tcf_a[43] = 1.0; Tcf_b[43] = 0.0; Tcb_a[43] = 1.0; Tcb_b[43] = 0.0;
    Tcf_a[44] = 1.0; Tcf_b[44] = 0.0; Tcb_a[44] = 1.0; Tcb_b[44] = 0.0;

    // Reaction 45: O+e- -> O++e-+e-  (electron impact ionization)
    Tcf_a[45] = 0.0; Tcf_b[45] = 1.0; Tcb_a[45] = 0.0; Tcb_b[45] = 1.0;

    // Reaction 46: N+e- -> N++e-+e-  (electron impact ionization)
    Tcf_a[46] = 0.0; Tcf_b[46] = 1.0; Tcb_a[46] = 0.0; Tcb_b[46] = 1.0;
    /*--- Collision integral data ---*/
    // Index 1: collider
    // Index 2: partner
    // Index 3: A1, A2, A3

    // Omega^(1,1) ----------------------
    // e-(0)
    Omega11(0,0,0) = -1.0000000E+00;  Omega11(0,0,1) = -1.0000000E+00;  Omega11(0,0,2) = -1.0000000E+00;  Omega11(0,0,3) = -1.0000000E+00;
    Omega11(0,1,0) = -1.0525124E-02;  Omega11(0,1,1) =  1.3498950E-01;  Omega11(0,1,2) =  1.2524805E-01;  Omega11(0,1,3) =  1.5066506E-01;
    Omega11(0,2,0) =  2.3527001E-02;  Omega11(0,2,1) = -6.9632323E-01;  Omega11(0,2,2) =  6.8035475E+00;  Omega11(0,2,3) =  1.8335509E-09;
    Omega11(0,3,0) =  1.0414818E-01;  Omega11(0,3,1) = -2.8369126E+00;  Omega11(0,3,2) =  2.5323135E+01;  Omega11(0,3,3) =  7.7138358E-32;
    Omega11(0,4,0) =  0.0000000E+00;  Omega11(0,4,1) =  1.6554247E-01;  Omega11(0,4,2) = -3.4986344E+00;  Omega11(0,4,3) =  5.9268038E+08;
    Omega11(0,5,0) =  9.9865506E-03;  Omega11(0,5,1) = -2.7407431E-01;  Omega11(0,5,2) =  2.6561032E+00;  Omega11(0,5,3) =  4.3080676E-04;
    Omega11(0,6,0) =  1.0000000E+00;  Omega11(0,6,1) =  1.0000000E+00;  Omega11(0,6,2) =  1.0000000E+00;  Omega11(0,6,3) =  1.0000000E+00;
    Omega11(0,7,0) =  1.0000000E+00;  Omega11(0,7,1) =  1.0000000E+00;  Omega11(0,7,2) =  1.0000000E+00;  Omega11(0,7,3) =  1.0000000E+00;
    Omega11(0,8,0) =  1.0000000E+00;  Omega11(0,8,1) =  1.0000000E+00;  Omega11(0,8,2) =  1.0000000E+00;  Omega11(0,8,3) =  1.0000000E+00;
    Omega11(0,9,0) =  1.0000000E+00;  Omega11(0,9,1) =  1.0000000E+00;  Omega11(0,9,2) =  1.0000000E+00;  Omega11(0,9,3) =  1.0000000E+00;
    Omega11(0,10,0) = 1.0000000E+00;  Omega11(0,10,1) = 1.0000000E+00;  Omega11(0,10,2) = 1.0000000E+00;  Omega11(0,10,3) = 1.0000000E+00;
    // N2(1)
    Omega11(1,0,0) = -1.0525124E-02;  Omega11(1,0,1) =  1.3498950E-01;  Omega11(1,0,2) =  1.2524805E-01;  Omega11(1,0,3) =  1.5066506E-01;
    Omega11(1,1,0) = -6.0614558E-03;  Omega11(1,1,1) =  1.2689102E-01;  Omega11(1,1,2) = -1.0616948E+00;  Omega11(1,1,3) =  8.0955466E+02;
    Omega11(1,2,0) = -3.7959091E-03;  Omega11(1,2,1) =  9.5708295E-02;  Omega11(1,2,2) = -1.0070611E+00;  Omega11(1,2,3) =  8.9392313E+02;
    Omega11(1,3,0) = -1.9295666E-03;  Omega11(1,3,1) =  2.7995735E-02;  Omega11(1,3,2) = -3.1588514E-01;  Omega11(1,3,3) =  1.2880734E+02;
    Omega11(1,4,0) = -1.0796249E-02;  Omega11(1,4,1) =  2.2656509E-01;  Omega11(1,4,2) = -1.7910602E+00;  Omega11(1,4,3) =  4.0455218E+03;
    Omega11(1,5,0) = -2.7244269E-03;  Omega11(1,5,1) =  6.9587171E-02;  Omega11(1,5,2) = -7.9538667E-01;  Omega11(1,5,3) =  4.0673730E+02;
    Omega11(1,6,0) =  0.0000000E+00;  Omega11(1,6,1) =  9.1205839E-02;  Omega11(1,6,2) = -1.8728231E+00;  Omega11(1,6,3) =  2.4432020E+05;
    Omega11(1,7,0) = -2.9123716E-03;  Omega11(1,7,1) =  9.6850678E-02;  Omega11(1,7,2) = -1.1416540E+00;  Omega11(1,7,3) =  7.9252169E+03;
    Omega11(1,8,0) =  1.2405624E-02;  Omega11(1,8,1) = -2.0452111E-01;  Omega11(1,8,2) =  3.5478475E-01;  Omega11(1,8,3) =  1.0778357E+03;
    Omega11(1,9,0) = -1.0687805E-02;  Omega11(1,9,1) =  2.4479697E-01;  Omega11(1,9,2) = -2.3192863E+00;  Omega11(1,9,3) =  1.0689229E+05;
    Omega11(1,10,0) = 1.0352091E-02;  Omega11(1,10,1) = -1.5733723E-01;  Omega11(1,10,2) =  2.9326150E-02;  Omega11(1,10,3) =  2.1003616E+03;
    // O2(2)
    Omega11(2,0,0) =  2.3527001E-02;  Omega11(2,0,1) = -6.9632323E-01;  Omega11(2,0,2) =  6.8035475E+00;  Omega11(2,0,3) =  1.8335509E-09;
    Omega11(2,1,0) = -3.7959091E-03;  Omega11(2,1,1) =  9.5708295E-02;  Omega11(2,1,2) = -1.0070611E+00;  Omega11(2,1,3) =  8.9392313E+02;
    Omega11(2,2,0) = -8.0682650E-04;  Omega11(2,2,1) =  1.6602480E-02;  Omega11(2,2,2) = -3.1472774E-01;  Omega11(2,2,3) =  1.4116458E+02;
    Omega11(2,3,0) = -6.4433840E-04;  Omega11(2,3,1) =  8.5378580E-03;  Omega11(2,3,2) = -2.3225102E-01;  Omega11(2,3,3) =  1.1371608E+02;
    Omega11(2,4,0) = -1.1453028E-03;  Omega11(2,4,1) =  1.2654140E-02;  Omega11(2,4,2) = -2.2435218E-01;  Omega11(2,4,3) =  7.7201588E+01;
    Omega11(2,5,0) = -4.8405803E-03;  Omega11(2,5,1) =  1.0297688E-01;  Omega11(2,5,2) = -9.6876576E-01;  Omega11(2,5,3) =  6.1629812E+02;
    Omega11(2,6,0) = -3.7822765E-03;  Omega11(2,6,1) =  1.7967016E-01;  Omega11(2,6,2) = -2.5409098E+00;  Omega11(2,6,3) =  1.1840435E+06;
    Omega11(2,7,0) = -4.0893007E-03;  Omega11(2,7,1) =  1.7795266E-01;  Omega11(2,7,2) = -2.3800543E+00;  Omega11(2,7,3) =  5.1949298E+05;
    Omega11(2,8,0) = -8.9520932E-03;  Omega11(2,8,1) =  2.2749642E-01;  Omega11(2,8,2) = -2.0758341E+00;  Omega11(2,8,3) =  6.7674419E+04;
    Omega11(2,9,0) =  0.0000000E+00;  Omega11(2,9,1) =  8.7745537E-02;  Omega11(2,9,2) = -1.8347158E+00;  Omega11(2,9,3) =  1.9830120E+05;
    Omega11(2,10,0) = 0.0000000E+00;  Omega11(2,10,1) =  9.3559978E-02;  Omega11(2,10,2) = -1.9842999E+00;  Omega11(2,10,3) =  4.3097490E+05;
    // NO(3)
    Omega11(3,0,0) =  1.0414818E-01;  Omega11(3,0,1) = -2.8369126E+00;  Omega11(3,0,2) =  2.5323135E+01;  Omega11(3,0,3) =  7.7138358E-32;
    Omega11(3,1,0) = -1.9295666E-03;  Omega11(3,1,1) =  2.7995735E-02;  Omega11(3,1,2) = -3.1588514E-01;  Omega11(3,1,3) =  1.2880734E+02;
    Omega11(3,2,0) = -6.4433840E-04;  Omega11(3,2,1) =  8.5378580E-03;  Omega11(3,2,2) = -2.3225102E-01;  Omega11(3,2,3) =  1.1371608E+02;
    Omega11(3,3,0) =  0.0000000E+00;  Omega11(3,3,1) = -1.1056066E-02;  Omega11(3,3,2) = -5.9216250E-02;  Omega11(3,3,3) =  7.2542367E+01;
    Omega11(3,4,0) = -1.5770918E-03;  Omega11(3,4,1) =  1.9578381E-02;  Omega11(3,4,2) = -2.7873624E-01;  Omega11(3,4,3) =  9.9547944E+01;
    Omega11(3,5,0) = -1.0885815E-03;  Omega11(3,5,1) =  1.1883688E-02;  Omega11(3,5,2) = -2.1844909E-01;  Omega11(3,5,3) =  7.5512560E+01;
    Omega11(3,6,0) = -8.1158474E-03;  Omega11(3,6,1) =  2.1474280E-01;  Omega11(3,6,2) = -2.0148450E+00;  Omega11(3,6,3) =  6.2986385E+04;
    Omega11(3,7,0) = -9.2292933E-03;  Omega11(3,7,1) =  2.9813226E-01;  Omega11(3,7,2) = -3.2899475E+00;  Omega11(3,7,3) =  4.9147046E+06;
    Omega11(3,8,0) =  1.3731123E-02;  Omega11(3,8,1) = -2.3920299E-01;  Omega11(3,8,2) =  6.7093226E-01;  Omega11(3,8,3) =  4.0068731E+02;
    Omega11(3,9,0) =  0.0000000E+00;  Omega11(3,9,1) =  8.6531037E-02;  Omega11(3,9,2) = -1.8117931E+00;  Omega11(3,9,3) =  1.8621272E+05;
    Omega11(3,10,0) = 8.3856973E-03;  Omega11(3,10,1) = -1.0972656E-01;  Omega11(3,10,2) = -3.3896281E-01;  Omega11(3,10,3) =  5.2004690E+03;
    // N(4)
    Omega11(4,0,0) =  0.0000000E+00;  Omega11(4,0,1) =  1.6554247E-01;  Omega11(4,0,2) = -3.4986344E+00;  Omega11(4,0,3) =  5.9268038E+08;
    Omega11(4,1,0) = -1.0796249E-02;  Omega11(4,1,1) =  2.2656509E-01;  Omega11(4,1,2) = -1.7910602E+00;  Omega11(4,1,3) =  4.0455218E+03;
    Omega11(4,2,0) = -1.1453028E-03;  Omega11(4,2,1) =  1.2654140E-02;  Omega11(4,2,2) = -2.2435218E-01;  Omega11(4,2,3) =  7.7201588E+01;
    Omega11(4,3,0) = -1.5770918E-03;  Omega11(4,3,1) =  1.9578381E-02;  Omega11(4,3,2) = -2.7873624E-01;  Omega11(4,3,3) =  9.9547944E+01;
    Omega11(4,4,0) = -9.6083779E-03;  Omega11(4,4,1) =  2.0938971E-01;  Omega11(4,4,2) = -1.7386904E+00;  Omega11(4,4,3) =  3.3587983E+03;
    Omega11(4,5,0) = -7.8147689E-03;  Omega11(4,5,1) =  1.6792705E-01;  Omega11(4,5,2) = -1.4308628E+00;  Omega11(4,5,3) =  1.6628859E+03;
    Omega11(4,6,0) = -1.9605234E-02;  Omega11(4,6,1) =  5.5570872E-01;  Omega11(4,6,2) = -5.4285702E+00;  Omega11(4,6,3) =  1.3574446E+09;
    Omega11(4,7,0) = -1.4501284E-02;  Omega11(4,7,1) =  4.1085338E-01;  Omega11(4,7,2) = -4.0115094E+00;  Omega11(4,7,3) =  1.5735451E+07;
    Omega11(4,8,0) =  0.0000000E+00;  Omega11(4,8,1) =  8.3065769E-02;  Omega11(4,8,2) = -1.7501512E+00;  Omega11(4,8,3) =  1.0846799E+05;
    Omega11(4,9,0) = -4.0078980E-03;  Omega11(4,9,1) =  1.0327487E-01;  Omega11(4,9,2) = -9.9473323E-01;  Omega11(4,9,3) =  2.8178290E+03;
    Omega11(4,10,0) = -2.4288224E-02;  Omega11(4,10,1) =  5.6305072E-01;  Omega11(4,10,2) = -4.6849679E+00;  Omega11(4,10,3) =  2.7303024E+07;
    // O(5)
    Omega11(5,0,0) =  9.9865506E-03;  Omega11(5,0,1) = -2.7407431E-01;  Omega11(5,0,2) =  2.6561032E+00;  Omega11(5,0,3) =  4.3080676E-04;
    Omega11(5,1,0) = -2.7244269E-03;  Omega11(5,1,1) =  6.9587171E-02;  Omega11(5,1,2) = -7.9538667E-01;  Omega11(5,1,3) =  4.0673730E+02;
    Omega11(5,2,0) = -4.8405803E-03;  Omega11(5,2,1) =  1.0297688E-01;  Omega11(5,2,2) = -9.6876576E-01;  Omega11(5,2,3) =  6.1629812E+02;
    Omega11(5,3,0) = -1.0885815E-03;  Omega11(5,3,1) =  1.1883688E-02;  Omega11(5,3,2) = -2.1844909E-01;  Omega11(5,3,3) =  7.5512560E+01;
    Omega11(5,4,0) = -7.8147689E-03;  Omega11(5,4,1) =  1.6792705E-01;  Omega11(5,4,2) = -1.4308628E+00;  Omega11(5,4,3) =  1.6628859E+03;
    Omega11(5,5,0) = -6.4040535E-03;  Omega11(5,5,1) =  1.4629949E-01;  Omega11(5,5,2) = -1.3892121E+00;  Omega11(5,5,3) =  2.0903441E+03;
    Omega11(5,6,0) = -1.6409054E-02;  Omega11(5,6,1) =  4.6352852E-01;  Omega11(5,6,2) = -4.5479735E+00;  Omega11(5,6,3) =  7.4250671E+07;
    Omega11(5,7,0) = -1.6923472E-02;  Omega11(5,7,1) =  4.6067692E-01;  Omega11(5,7,2) = -4.3294966E+00;  Omega11(5,7,3) =  2.5538927E+07;
    Omega11(5,8,0) = -2.9417970E-03;  Omega11(5,8,1) =  1.5129273E-01;  Omega11(5,8,2) = -2.2497964E+00;  Omega11(5,8,3) =  2.9325215E+05;
    Omega11(5,9,0) = -1.5767256E-02;  Omega11(5,9,1) =  3.5405830E-01;  Omega11(5,9,2) = -3.0686783E+00;  Omega11(5,9,3) =  4.6336779E+05;
    Omega11(5,10,0) = -3.8347988E-03;  Omega11(5,10,1) =  9.9930498E-02;  Omega11(5,10,2) = -9.6288891E-01;  Omega11(5,10,3) =  1.9669897E+03;
    // NO+(6)
    Omega11(6,0,0) =  1.0000000E+00;  Omega11(6,0,1) =  1.0000000E+00;  Omega11(6,0,2) =  1.0000000E+00;  Omega11(6,0,3) =  1.0000000E+00;
    Omega11(6,1,0) =  0.0000000E+00;  Omega11(6,1,1) =  9.1205839E-02;  Omega11(6,1,2) = -1.8728231E+00;  Omega11(6,1,3) =  2.4432020E+05;
    Omega11(6,2,0) = -3.7822765E-03;  Omega11(6,2,1) =  1.7967016E-01;  Omega11(6,2,2) = -2.5409098E+00;  Omega11(6,2,3) =  1.1840435E+06;
    Omega11(6,3,0) = -8.1158474E-03;  Omega11(6,3,1) =  2.1474280E-01;  Omega11(6,3,2) = -2.0148450E+00;  Omega11(6,3,3) =  6.2986385E+04;
    Omega11(6,4,0) = -1.9605234E-02;  Omega11(6,4,1) =  5.5570872E-01;  Omega11(6,4,2) = -5.4285702E+00;  Omega11(6,4,3) =  1.3574446E+09;
    Omega11(6,5,0) = -1.6409054E-02;  Omega11(6,5,1) =  4.6352852E-01;  Omega11(6,5,2) = -4.5479735E+00;  Omega11(6,5,3) =  7.4250671E+07;
    Omega11(6,6,0) = -1.0000000E+00;  Omega11(6,6,1) = -1.0000000E+00;  Omega11(6,6,2) = -1.0000000E+00;  Omega11(6,6,3) = -1.0000000E+00;
    Omega11(6,7,0) = -1.0000000E+00;  Omega11(6,7,1) = -1.0000000E+00;  Omega11(6,7,2) = -1.0000000E+00;  Omega11(6,7,3) = -1.0000000E+00;
    Omega11(6,8,0) = -1.0000000E+00;  Omega11(6,8,1) = -1.0000000E+00;  Omega11(6,8,2) = -1.0000000E+00;  Omega11(6,8,3) = -1.0000000E+00;
    Omega11(6,9,0) = -1.0000000E+00;  Omega11(6,9,1) = -1.0000000E+00;  Omega11(6,9,2) = -1.0000000E+00;  Omega11(6,9,3) = -1.0000000E+00;
    Omega11(6,10,0) = -1.0000000E+00; Omega11(6,10,1) = -1.0000000E+00; Omega11(6,10,2) = -1.0000000E+00; Omega11(6,10,3) = -1.0000000E+00;
    // N2+(7)
    Omega11(7,0,0) =  1.0000000E+00;  Omega11(7,0,1) =  1.0000000E+00;  Omega11(7,0,2) =  1.0000000E+00;  Omega11(7,0,3) =  1.0000000E+00;
    Omega11(7,1,0) = -2.9123716E-03;  Omega11(7,1,1) =  9.6850678E-02;  Omega11(7,1,2) = -1.1416540E+00;  Omega11(7,1,3) =  7.9252169E+03;
    Omega11(7,2,0) = -4.0893007E-03;  Omega11(7,2,1) =  1.7795266E-01;  Omega11(7,2,2) = -2.3800543E+00;  Omega11(7,2,3) =  5.1949298E+05;
    Omega11(7,3,0) = -9.2292933E-03;  Omega11(7,3,1) =  2.9813226E-01;  Omega11(7,3,2) = -3.2899475E+00;  Omega11(7,3,3) =  4.9147046E+06;
    Omega11(7,4,0) = -1.4501284E-02;  Omega11(7,4,1) =  4.1085338E-01;  Omega11(7,4,2) = -4.0115094E+00;  Omega11(7,4,3) =  1.5735451E+07;
    Omega11(7,5,0) = -1.6923472E-02;  Omega11(7,5,1) =  4.6067692E-01;  Omega11(7,5,2) = -4.3294966E+00;  Omega11(7,5,3) =  2.5538927E+07;
    Omega11(7,6,0) = -1.0000000E+00;  Omega11(7,6,1) = -1.0000000E+00;  Omega11(7,6,2) = -1.0000000E+00;  Omega11(7,6,3) = -1.0000000E+00;
    Omega11(7,7,0) = -1.0000000E+00;  Omega11(7,7,1) = -1.0000000E+00;  Omega11(7,7,2) = -1.0000000E+00;  Omega11(7,7,3) = -1.0000000E+00;
    Omega11(7,8,0) = -1.0000000E+00;  Omega11(7,8,1) = -1.0000000E+00;  Omega11(7,8,2) = -1.0000000E+00;  Omega11(7,8,3) = -1.0000000E+00;
    Omega11(7,9,0) = -1.0000000E+00;  Omega11(7,9,1) = -1.0000000E+00;  Omega11(7,9,2) = -1.0000000E+00;  Omega11(7,9,3) = -1.0000000E+00;
    Omega11(7,10,0) = -1.0000000E+00; Omega11(7,10,1) = -1.0000000E+00; Omega11(7,10,2) = -1.0000000E+00; Omega11(7,10,3) = -1.0000000E+00;
    // O2+(8)
    Omega11(8,0,0) =  1.0000000E+00;  Omega11(8,0,1) =  1.0000000E+00;  Omega11(8,0,2) =  1.0000000E+00;  Omega11(8,0,3) =  1.0000000E+00;
    Omega11(8,1,0) =  1.2405624E-02;  Omega11(8,1,1) = -2.0452111E-01;  Omega11(8,1,2) =  3.5478475E-01;  Omega11(8,1,3) =  1.0778357E+03;
    Omega11(8,2,0) = -8.9520932E-03;  Omega11(8,2,1) =  2.2749642E-01;  Omega11(8,2,2) = -2.0758341E+00;  Omega11(8,2,3) =  6.7674419E+04;
    Omega11(8,3,0) =  1.3731123E-02;  Omega11(8,3,1) = -2.3920299E-01;  Omega11(8,3,2) =  6.7093226E-01;  Omega11(8,3,3) =  4.0068731E+02;
    Omega11(8,4,0) =  0.0000000E+00;  Omega11(8,4,1) =  8.3065769E-02;  Omega11(8,4,2) = -1.7501512E+00;  Omega11(8,4,3) =  1.0846799E+05;
    Omega11(8,5,0) = -2.9417970E-03;  Omega11(8,5,1) =  1.5129273E-01;  Omega11(8,5,2) = -2.2497964E+00;  Omega11(8,5,3) =  2.9325215E+05;
    Omega11(8,6,0) = -1.0000000E+00;  Omega11(8,6,1) = -1.0000000E+00;  Omega11(8,6,2) = -1.0000000E+00;  Omega11(8,6,3) = -1.0000000E+00;
    Omega11(8,7,0) = -1.0000000E+00;  Omega11(8,7,1) = -1.0000000E+00;  Omega11(8,7,2) = -1.0000000E+00;  Omega11(8,7,3) = -1.0000000E+00;
    Omega11(8,8,0) = -1.0000000E+00;  Omega11(8,8,1) = -1.0000000E+00;  Omega11(8,8,2) = -1.0000000E+00;  Omega11(8,8,3) = -1.0000000E+00;
    Omega11(8,9,0) = -1.0000000E+00;  Omega11(8,9,1) = -1.0000000E+00;  Omega11(8,9,2) = -1.0000000E+00;  Omega11(8,9,3) = -1.0000000E+00;
    Omega11(8,10,0) = -1.0000000E+00; Omega11(8,10,1) = -1.0000000E+00; Omega11(8,10,2) = -1.0000000E+00; Omega11(8,10,3) = -1.0000000E+00;
    // N+(9)
    Omega11(9,0,0) =  1.0000000E+00;  Omega11(9,0,1) =  1.0000000E+00;  Omega11(9,0,2) =  1.0000000E+00;  Omega11(9,0,3) =  1.0000000E+00;
    Omega11(9,1,0) = -1.0687805E-02;  Omega11(9,1,1) =  2.4479697E-01;  Omega11(9,1,2) = -2.3192863E+00;  Omega11(9,1,3) =  1.0689229E+05;
    Omega11(9,2,0) =  0.0000000E+00;  Omega11(9,2,1) =  8.7745537E-02;  Omega11(9,2,2) = -1.8347158E+00;  Omega11(9,2,3) =  1.9830120E+05;
    Omega11(9,3,0) =  0.0000000E+00;  Omega11(9,3,1) =  8.6531037E-02;  Omega11(9,3,2) = -1.8117931E+00;  Omega11(9,3,3) =  1.8621272E+05;
    Omega11(9,4,0) = -4.0078980E-03;  Omega11(9,4,1) =  1.0327487E-01;  Omega11(9,4,2) = -9.9473323E-01;  Omega11(9,4,3) =  2.8178290E+03;
    Omega11(9,5,0) = -1.5767256E-02;  Omega11(9,5,1) =  3.5405830E-01;  Omega11(9,5,2) = -3.0686783E+00;  Omega11(9,5,3) =  4.6336779E+05;
    Omega11(9,6,0) = -1.0000000E+00;  Omega11(9,6,1) = -1.0000000E+00;  Omega11(9,6,2) = -1.0000000E+00;  Omega11(9,6,3) = -1.0000000E+00;
    Omega11(9,7,0) = -1.0000000E+00;  Omega11(9,7,1) = -1.0000000E+00;  Omega11(9,7,2) = -1.0000000E+00;  Omega11(9,7,3) = -1.0000000E+00;
    Omega11(9,8,0) = -1.0000000E+00;  Omega11(9,8,1) = -1.0000000E+00;  Omega11(9,8,2) = -1.0000000E+00;  Omega11(9,8,3) = -1.0000000E+00;
    Omega11(9,9,0) = -1.0000000E+00;  Omega11(9,9,1) = -1.0000000E+00;  Omega11(9,9,2) = -1.0000000E+00;  Omega11(9,9,3) = -1.0000000E+00;
    Omega11(9,10,0) = -1.0000000E+00; Omega11(9,10,1) = -1.0000000E+00; Omega11(9,10,2) = -1.0000000E+00; Omega11(9,10,3) = -1.0000000E+00;
    // O+(10)
    Omega11(10,0,0) =  1.0000000E+00;  Omega11(10,0,1) =  1.0000000E+00;  Omega11(10,0,2) =  1.0000000E+00;  Omega11(10,0,3) =  1.0000000E+00;
    Omega11(10,1,0) =  1.0352091E-02;  Omega11(10,1,1) = -1.5733723E-01;  Omega11(10,1,2) =  2.9326150E-02;  Omega11(10,1,3) =  2.1003616E+03;
    Omega11(10,2,0) =  0.0000000E+00;  Omega11(10,2,1) =  9.3559978E-02;  Omega11(10,2,2) = -1.9842999E+00;  Omega11(10,2,3) =  4.3097490E+05;
    Omega11(10,3,0) =  8.3856973E-03;  Omega11(10,3,1) = -1.0972656E-01;  Omega11(10,3,2) = -3.3896281E-01;  Omega11(10,3,3) =  5.2004690E+03;
    Omega11(10,4,0) = -2.4288224E-02;  Omega11(10,4,1) =  5.6305072E-01;  Omega11(10,4,2) = -4.6849679E+00;  Omega11(10,4,3) =  2.7303024E+07;
    Omega11(10,5,0) = -3.8347988E-03;  Omega11(10,5,1) =  9.9930498E-02;  Omega11(10,5,2) = -9.6288891E-01;  Omega11(10,5,3) =  1.9669897E+03;
    Omega11(10,6,0) = -1.0000000E+00;  Omega11(10,6,1) = -1.0000000E+00;  Omega11(10,6,2) = -1.0000000E+00;  Omega11(10,6,3) = -1.0000000E+00;
    Omega11(10,7,0) = -1.0000000E+00;  Omega11(10,7,1) = -1.0000000E+00;  Omega11(10,7,2) = -1.0000000E+00;  Omega11(10,7,3) = -1.0000000E+00;
    Omega11(10,8,0) = -1.0000000E+00;  Omega11(10,8,1) = -1.0000000E+00;  Omega11(10,8,2) = -1.0000000E+00;  Omega11(10,8,3) = -1.0000000E+00;
    Omega11(10,9,0) = -1.0000000E+00;  Omega11(10,9,1) = -1.0000000E+00;  Omega11(10,9,2) = -1.0000000E+00;  Omega11(10,9,3) = -1.0000000E+00;
    Omega11(10,10,0) = -1.0000000E+00; Omega11(10,10,1) = -1.0000000E+00; Omega11(10,10,2) = -1.0000000E+00; Omega11(10,10,3) = -1.0000000E+00;


    // Omega^(2,2) ----------------------
    // e-(0)
    Omega22(0,0,0) = -1.0000000E+00;  Omega22(0,0,1) = -1.0000000E+00;  Omega22(0,0,2) = -1.0000000E+00;  Omega22(0,0,3) = -1.0000000E+00;
    Omega22(0,1,0) = -4.2254948E-03;  Omega22(0,1,1) = -5.2965163E-02;  Omega22(0,1,2) =  1.9157708E+00;  Omega22(0,1,3) =  6.3263309E-04;
    Omega22(0,2,0) =  9.6744867E-03;  Omega22(0,2,1) = -3.3759583E-01;  Omega22(0,2,2) =  3.7952121E+00;  Omega22(0,2,3) =  6.8468036E-06;
    Omega22(0,3,0) =  0.0000000E+00;  Omega22(0,3,1) =  5.4444485E-02;  Omega22(0,3,2) = -1.2854128E+00;  Omega22(0,3,3) =  1.3857556E+04;
    Omega22(0,4,0) = -1.0903638E-01;  Omega22(0,4,1) =  2.8678381E+00;  Omega22(0,4,2) = -2.5297550E+01;  Omega22(0,4,3) =  3.4838798E+33;
    Omega22(0,5,0) = -1.7924100E-02;  Omega22(0,5,1) =  4.0402656E-01;  Omega22(0,5,2) = -2.6712374E+00;  Omega22(0,5,3) =  4.1447669E+02;
    Omega22(0,6,0) =  1.0000000E+00;  Omega22(0,6,1) =  1.0000000E+00;  Omega22(0,6,2) =  1.0000000E+00;  Omega22(0,6,3) =  1.0000000E+00;
    Omega22(0,7,0) =  1.0000000E+00;  Omega22(0,7,1) =  1.0000000E+00;  Omega22(0,7,2) =  1.0000000E+00;  Omega22(0,7,3) =  1.0000000E+00;
    Omega22(0,8,0) =  1.0000000E+00;  Omega22(0,8,1) =  1.0000000E+00;  Omega22(0,8,2) =  1.0000000E+00;  Omega22(0,8,3) =  1.0000000E+00;
    Omega22(0,9,0) =  1.0000000E+00;  Omega22(0,9,1) =  1.0000000E+00;  Omega22(0,9,2) =  1.0000000E+00;  Omega22(0,9,3) =  1.0000000E+00;
    Omega22(0,10,0) = 1.0000000E+00;  Omega22(0,10,1) = 1.0000000E+00;  Omega22(0,10,2) = 1.0000000E+00;  Omega22(0,10,3) = 1.0000000E+00;
    // N2(1)
    Omega22(1,0,0) = -4.2254948E-03;  Omega22(1,0,1) = -5.2965163E-02;  Omega22(1,0,2) =  1.9157708E+00;  Omega22(1,0,3) =  6.3263309E-04;
    Omega22(1,1,0) = -7.6303990E-03;  Omega22(1,1,1) =  1.6878089E-01;  Omega22(1,1,2) = -1.4004234E+00;  Omega22(1,1,3) =  2.1427708E+03;
    Omega22(1,2,0) = -8.0457321E-03;  Omega22(1,2,1) =  1.9228905E-01;  Omega22(1,2,2) = -1.7102854E+00;  Omega22(1,2,3) =  5.2213857E+03;
    Omega22(1,3,0) = -6.8237776E-03;  Omega22(1,3,1) =  1.4360616E-01;  Omega22(1,3,2) = -1.1922240E+00;  Omega22(1,3,3) =  1.2433086E+03;
    Omega22(1,4,0) = -8.3493693E-03;  Omega22(1,4,1) =  1.7808911E-01;  Omega22(1,4,2) = -1.4466155E+00;  Omega22(1,4,3) =  1.9324210E+03;
    Omega22(1,5,0) = -8.3110691E-03;  Omega22(1,5,1) =  1.9617877E-01;  Omega22(1,5,2) = -1.7205427E+00;  Omega22(1,5,3) =  4.0812829E+03;
    Omega22(1,6,0) =  0.0000000E+00;  Omega22(1,6,1) =  8.5112236E-02;  Omega22(1,6,2) = -1.7460044E+00;  Omega22(1,6,3) =  1.4498969E+05;
    Omega22(1,7,0) = -1.6447237E-02;  Omega22(1,7,1) =  4.7759522E-01;  Omega22(1,7,2) = -4.7641986E+00;  Omega22(1,7,3) =  2.9127542E+08;
    Omega22(1,8,0) =  2.2455421E-02;  Omega22(1,8,1) = -4.5106797E-01;  Omega22(1,8,2) =  2.3763420E+00;  Omega22(1,8,3) =  4.7754696E+00;
    Omega22(1,9,0) = -7.0776069E-03;  Omega22(1,9,1) =  1.7917938E-01;  Omega22(1,9,2) = -1.9102410E+00;  Omega22(1,9,3) =  4.6736263E+04;
    Omega22(1,10,0) = 1.8733000E-02;  Omega22(1,10,1)= -3.6163781E-01;  Omega22(1,10,2) = 1.6947101E+00;  Omega22(1,10,3) =  2.5244859E+01;
    // O2(2)
    Omega22(2,0,0) =  9.6744867E-03;  Omega22(2,0,1) = -3.3759583E-01;  Omega22(2,0,2) =  3.7952121E+00;  Omega22(2,0,3) =  6.8468036E-06;
    Omega22(2,1,0) = -8.0457321E-03;  Omega22(2,1,1) =  1.9228905E-01;  Omega22(2,1,2) = -1.7102854E+00;  Omega22(2,1,3) =  5.2213857E+03;
    Omega22(2,2,0) = -6.2931612E-03;  Omega22(2,2,1) =  1.4624645E-01;  Omega22(2,2,2) = -1.3006927E+00;  Omega22(2,2,3) =  1.8066892E+03;
    Omega22(2,3,0) = -6.8508672E-03;  Omega22(2,3,1) =  1.5524564E-01;  Omega22(2,3,2) = -1.3479583E+00;  Omega22(2,3,3) =  2.0037890E+03;
    Omega22(2,4,0) = -1.0608832E-03;  Omega22(2,4,1) =  1.1782595E-02;  Omega22(2,4,2) = -2.1246301E-01;  Omega22(2,4,3) =  8.4561598E+01;
    Omega22(2,5,0) = -3.7969686E-03;  Omega22(2,5,1) =  7.6789981E-02;  Omega22(2,5,2) = -7.3056809E-01;  Omega22(2,5,3) =  3.3958171E+02;
    Omega22(2,6,0) =  0.0000000E+00;  Omega22(2,6,1) =  8.4737359E-02;  Omega22(2,6,2) = -1.7290488E+00;  Omega22(2,6,3) =  1.2485194E+05;
    Omega22(2,7,0) = -4.9176811E-03;  Omega22(2,7,1) =  1.9694738E-01;  Omega22(2,7,2) = -2.5025540E+00;  Omega22(2,7,3) =  6.8213629E+05;
    Omega22(2,8,0) =  2.8664463E-02;  Omega22(2,8,1) = -5.8087240E-01;  Omega22(2,8,2) =  3.2564558E+00;  Omega22(2,8,3) =  6.6890428E-01;
    Omega22(2,9,0) =  9.8019578E-03;  Omega22(2,9,1) = -1.4699425E-01;  Omega22(2,9,2) =  3.9382460E-02;  Omega22(2,9,3) =  1.5112165E+03;
    Omega22(2,10,0) = 1.4207970E-02;  Omega22(2,10,1) = -2.4736726E-01;  Omega22(2,10,2) =  7.4561859E-01;  Omega22(2,10,3) =  3.2519188E+02;
    // NO(3)
    Omega22(3,0,0) =  0.0000000E+00;  Omega22(3,0,1) =  5.4444485E-02;  Omega22(3,0,2) = -1.2854128E+00;  Omega22(3,0,3) =  1.3857556E+04;
    Omega22(3,1,0) = -6.8237776E-03;  Omega22(3,1,1) =  1.4360616E-01;  Omega22(3,1,2) = -1.1922240E+00;  Omega22(3,1,3) =  1.2433086E+03;
    Omega22(3,2,0) = -6.8508672E-03;  Omega22(3,2,1) =  1.5524564E-01;  Omega22(3,2,2) = -1.3479583E+00;  Omega22(3,2,3) =  2.0037890E+03;
    Omega22(3,3,0) = -7.4942466E-03;  Omega22(3,3,1) =  1.6626193E-01;  Omega22(3,3,2) = -1.4107027E+00;  Omega22(3,3,3) =  2.3097604E+03;
    Omega22(3,4,0) = -1.4719259E-03;  Omega22(3,4,1) =  1.8446968E-02;  Omega22(3,4,2) = -2.6460411E-01;  Omega22(3,4,3) =  1.0911124E+02;
    Omega22(3,5,0) = -1.0066279E-03;  Omega22(3,5,1) =  1.1029264E-02;  Omega22(3,5,2) = -2.0671266E-01;  Omega22(3,5,3) =  8.2644384E+01;
    Omega22(3,6,0) =  1.1055777E-02;  Omega22(3,6,1) = -1.6621846E-01;  Omega22(3,6,2) =  1.4372166E-01;  Omega22(3,6,3) =  1.3182061E+03;
    Omega22(3,7,0) = -4.0133981E-03;  Omega22(3,7,1) =  1.7290664E-01;  Omega22(3,7,2) = -2.2855449E+00;  Omega22(3,7,3) =  3.6429320E+05;
    Omega22(3,8,0) =  2.2679271E-02;  Omega22(3,8,1) = -4.5710920E-01;  Omega22(3,8,2) =  2.4427275E+00;  Omega22(3,8,3) =  3.6733514E+00;
    Omega22(3,9,0) =  1.1716366E-02;  Omega22(3,9,1) = -1.9289789E-01;  Omega22(3,9,2) =  4.0269474E-01;  Omega22(3,9,3) =  6.0891590E+02;
    Omega22(3,10,0) = 1.8015337E-02;  Omega22(3,10,1) = -3.4415293E-01;  Omega22(3,10,2) = 1.5658151E+00 ;  Omega22(3,10,3) =  3.3303758E+01;
    // N(4)
    Omega22(4,0,0) = -1.0903638E-01;  Omega22(4,0,1) =  2.8678381E+00;  Omega22(4,0,2) = -2.5297550E+01;  Omega22(4,0,3) =  3.4838798E+33;
    Omega22(4,1,0) = -8.3493693E-03;  Omega22(4,1,1) =  1.7808911E-01;  Omega22(4,1,2) = -1.4466155E+00;  Omega22(4,1,3) =  1.9324210E+03;
    Omega22(4,2,0) = -1.0608832E-03;  Omega22(4,2,1) =  1.1782595E-02;  Omega22(4,2,2) = -2.1246301E-01;  Omega22(4,2,3) =  8.4561598E+01;
    Omega22(4,3,0) = -1.4719259E-03;  Omega22(4,3,1) =  1.8446968E-02;  Omega22(4,3,2) = -2.6460411E-01;  Omega22(4,3,3) =  1.0911124E+02;
    Omega22(4,4,0) = -7.7439615E-03;  Omega22(4,4,1) =  1.7129007E-01;  Omega22(4,4,2) = -1.4809088E+00;  Omega22(4,4,3) =  2.1284951E+03;
    Omega22(4,5,0) = -5.0478143E-03;  Omega22(4,5,1) =  1.0236186E-01;  Omega22(4,5,2) = -9.0058935E-01;  Omega22(4,5,3) =  4.4472565E+02;
    Omega22(4,6,0) = -2.1009546E-02;  Omega22(4,6,1) =  5.8910426E-01;  Omega22(4,6,2) = -5.6681361E+00;  Omega22(4,6,3) =  2.4486594E+09;
    Omega22(4,7,0) = -1.2882395E-02;  Omega22(4,7,1) =  3.7306469E-01;  Omega22(4,7,2) = -3.7106760E+00;  Omega22(4,7,3) =  7.5444981E+06;
    Omega22(4,8,0) =  1.1205000E-02;  Omega22(4,8,1) = -1.8182149E-01;  Omega22(4,8,2) =  3.2624972E-01;  Omega22(4,8,3) =  5.5186183E+02;
    Omega22(4,9,0) = -1.4271306E-02;  Omega22(4,9,1) =  3.0401993E-01;  Omega22(4,9,2) = -2.4573879E+00;  Omega22(4,9,3) =  5.3694705E+04;
    Omega22(4,10,0)= -2.1681211E-02;  Omega22(4,10,1) = 5.2300453E-01;  Omega22(4,10,2)= -4.5118623E+00;  Omega22(4,10,3) =  2.3467766E+07;
    // O(5)
    Omega22(5,0,0) = -1.7924100E-02;  Omega22(5,0,1) =  4.0402656E-01;  Omega22(5,0,2) = -2.6712374E+00;  Omega22(5,0,3) =  4.1447669E+02;
    Omega22(5,1,0) = -8.3110691E-03;  Omega22(5,1,1) =  1.9617877E-01;  Omega22(5,1,2) = -1.7205427E+00;  Omega22(5,1,3) =  4.0812829E+03;
    Omega22(5,2,0) = -3.7969686E-03;  Omega22(5,2,1) =  7.6789981E-02;  Omega22(5,2,2) = -7.3056809E-01;  Omega22(5,2,3) =  3.3958171E+02;
    Omega22(5,3,0) = -1.0066279E-03;  Omega22(5,3,1) =  1.1029264E-02;  Omega22(5,3,2) = -2.0671266E-01;  Omega22(5,3,3) =  8.2644384E+01;
    Omega22(5,4,0) = -5.0478143E-03;  Omega22(5,4,1) =  1.0236186E-01;  Omega22(5,4,2) = -9.0058935E-01;  Omega22(5,4,3) =  4.4472565E+02;
    Omega22(5,5,0) = -4.2451096E-03;  Omega22(5,5,1) =  9.6820337E-02;  Omega22(5,5,2) = -9.9770795E-01;  Omega22(5,5,3) =  8.3320644E+02;
    Omega22(5,6,0) = -1.5315132E-02;  Omega22(5,6,1) =  4.3541627E-01;  Omega22(5,6,2) = -4.2864279E+00;  Omega22(5,6,3) =  3.5125207E+07;
    Omega22(5,7,0) = -1.7420606E-02;  Omega22(5,7,1) =  4.7126950E-01;  Omega22(5,7,2) = -4.3841087E+00;  Omega22(5,7,3) =  2.8275095E+07;
    Omega22(5,8,0) =  0.0000000E+00;  Omega22(5,8,1) =  8.3446262E-02;  Omega22(5,8,2) = -1.7191179E+00;  Omega22(5,8,3) =  8.0539928E+04;
    Omega22(5,9,0) = -1.7907392E-02;  Omega22(5,9,1) =  4.1207892E-01;  Omega22(5,9,2) = -3.5343610E+00;  Omega22(5,9,3) =  1.4987678E+06;
    Omega22(5,10,0)= -1.6032919E-02;  Omega22(5,10,1) = 3.7114396E-01;  Omega22(5,10,2)= -3.2050078E+00;  Omega22(5,10,3) =  5.8099314E+05;
    // NO+(6)
    Omega22(6,0,0) =  1.0000000E+00;  Omega22(6,0,1) =  1.0000000E+00;  Omega22(6,0,2) =  1.0000000E+00;  Omega22(6,0,3) =  1.0000000E+00;
    Omega22(6,1,0) =  0.0000000E+00;  Omega22(6,1,1) =  8.5112236E-02;  Omega22(6,1,2) = -1.7460044E+00;  Omega22(6,1,3) =  1.4498969E+05;
    Omega22(6,2,0) =  0.0000000E+00;  Omega22(6,2,1) =  8.4737359E-02;  Omega22(6,2,2) = -1.7290488E+00;  Omega22(6,2,3) =  1.2485194E+05;
    Omega22(6,3,0) =  1.1055777E-02;  Omega22(6,3,1) = -1.6621846E-01;  Omega22(6,3,2) =  1.4372166E-01;  Omega22(6,3,3) =  1.3182061E+03;
    Omega22(6,4,0) = -2.1009546E-02;  Omega22(6,4,1) =  5.8910426E-01;  Omega22(6,4,2) = -5.6681361E+00;  Omega22(6,4,3) =  2.4486594E+09;
    Omega22(6,5,0) = -1.5315132E-02;  Omega22(6,5,1) =  4.3541627E-01;  Omega22(6,5,2) = -4.2864279E+00;  Omega22(6,5,3) =  3.5125207E+07;
    Omega22(6,6,0) = -1.0000000E+00;  Omega22(6,6,1) = -1.0000000E+00;  Omega22(6,6,2) = -1.0000000E+00;  Omega22(6,6,3) = -1.0000000E+00;
    Omega22(6,7,0) = -1.0000000E+00;  Omega22(6,7,1) = -1.0000000E+00;  Omega22(6,7,2) = -1.0000000E+00;  Omega22(6,7,3) = -1.0000000E+00;
    Omega22(6,8,0) = -1.0000000E+00;  Omega22(6,8,1) = -1.0000000E+00;  Omega22(6,8,2) = -1.0000000E+00;  Omega22(6,8,3) = -1.0000000E+00;
    Omega22(6,9,0) = -1.0000000E+00;  Omega22(6,9,1) = -1.0000000E+00;  Omega22(6,9,2) = -1.0000000E+00;  Omega22(6,9,3) = -1.0000000E+00;
    Omega22(6,10,0) = -1.0000000E+00; Omega22(6,10,1) = -1.0000000E+00; Omega22(6,10,2) = -1.0000000E+00; Omega22(6,10,3) = -1.0000000E+00;
    // N2+(7)
    Omega22(7,0,0) =  1.0000000E+00;  Omega22(7,0,1) =  1.0000000E+00;  Omega22(7,0,2) =  1.0000000E+00;  Omega22(7,0,3) =  1.0000000E+00;
    Omega22(7,1,0) = -1.6447237E-02;  Omega22(7,1,1) =  4.7759522E-01;  Omega22(7,1,2) = -4.7641986E+00;  Omega22(7,1,3) =  2.9127542E+08;
    Omega22(7,2,0) = -4.9176811E-03;  Omega22(7,2,1) =  1.9694738E-01;  Omega22(7,2,2) = -2.5025540E+00;  Omega22(7,2,3) =  6.8213629E+05;
    Omega22(7,3,0) = -4.0133981E-03;  Omega22(7,3,1) =  1.7290664E-01;  Omega22(7,3,2) = -2.2855449E+00;  Omega22(7,3,3) =  3.6429320E+05;
    Omega22(7,4,0) = -1.2882395E-02;  Omega22(7,4,1) =  3.7306469E-01;  Omega22(7,4,2) = -3.7106760E+00;  Omega22(7,4,3) =  7.5444981E+06;
    Omega22(7,5,0) = -1.7420606E-02;  Omega22(7,5,1) =  4.7126950E-01;  Omega22(7,5,2) = -4.3841087E+00;  Omega22(7,5,3) =  2.8275095E+07;
    Omega22(7,6,0) = -1.0000000E+00;  Omega22(7,6,1) = -1.0000000E+00;  Omega22(7,6,2) = -1.0000000E+00;  Omega22(7,6,3) = -1.0000000E+00;
    Omega22(7,7,0) = -1.0000000E+00;  Omega22(7,7,1) = -1.0000000E+00;  Omega22(7,7,2) = -1.0000000E+00;  Omega22(7,7,3) = -1.0000000E+00;
    Omega22(7,8,0) = -1.0000000E+00;  Omega22(7,8,1) = -1.0000000E+00;  Omega22(7,8,2) = -1.0000000E+00;  Omega22(7,8,3) = -1.0000000E+00;
    Omega22(7,9,0) = -1.0000000E+00;  Omega22(7,9,1) = -1.0000000E+00;  Omega22(7,9,2) = -1.0000000E+00;  Omega22(7,9,3) = -1.0000000E+00;
    Omega22(7,10,0) = -1.0000000E+00; Omega22(7,10,1) = -1.0000000E+00; Omega22(7,10,2) = -1.0000000E+00; Omega22(7,10,3) = -1.0000000E+00;
    // O2+(8)
    Omega22(8,0,0) =  1.0000000E+00;  Omega22(8,0,1) =  1.0000000E+00;  Omega22(8,0,2) =  1.0000000E+00;  Omega22(8,0,3) =  1.0000000E+00;
    Omega22(8,1,0) =  2.2455421E-02;  Omega22(8,1,1) = -4.5106797E-01;  Omega22(8,1,2) =  2.3763420E+00;  Omega22(8,1,3) =  4.7754696E+00;
    Omega22(8,2,0) =  2.8664463E-02;  Omega22(8,2,1) = -5.8087240E-01;  Omega22(8,2,2) =  3.2564558E+00;  Omega22(8,2,3) =  6.6890428E-01;
    Omega22(8,3,0) =  2.2679271E-02;  Omega22(8,3,1) = -4.5710920E-01;  Omega22(8,3,2) =  2.4427275E+00;  Omega22(8,3,3) =  3.6733514E+00;
    Omega22(8,4,0) =  1.1205000E-02;  Omega22(8,4,1) = -1.8182149E-01;  Omega22(8,4,2) =  3.2624972E-01;  Omega22(8,4,3) =  5.5186183E+02;
    Omega22(8,5,0) =  0.0000000E+00;  Omega22(8,5,1) =  8.3446262E-02;  Omega22(8,5,2) = -1.7191179E+00;  Omega22(8,5,3) =  8.0539928E+04;
    Omega22(8,6,0) = -1.0000000E+00;  Omega22(8,6,1) = -1.0000000E+00;  Omega22(8,6,2) = -1.0000000E+00;  Omega22(8,6,3) = -1.0000000E+00;
    Omega22(8,7,0) = -1.0000000E+00;  Omega22(8,7,1) = -1.0000000E+00;  Omega22(8,7,2) = -1.0000000E+00;  Omega22(8,7,3) = -1.0000000E+00;
    Omega22(8,8,0) = -1.0000000E+00;  Omega22(8,8,1) = -1.0000000E+00;  Omega22(8,8,2) = -1.0000000E+00;  Omega22(8,8,3) = -1.0000000E+00;
    Omega22(8,9,0) = -1.0000000E+00;  Omega22(8,9,1) = -1.0000000E+00;  Omega22(8,9,2) = -1.0000000E+00;  Omega22(8,9,3) = -1.0000000E+00;
    Omega22(8,10,0) = -1.0000000E+00; Omega22(8,10,1) = -1.0000000E+00; Omega22(8,10,2) = -1.0000000E+00; Omega22(8,10,3) = -1.0000000E+00;
    // N+(9)
    Omega22(9,0,0) =  1.0000000E+00;  Omega22(9,0,1) =  1.0000000E+00;  Omega22(9,0,2) =  1.0000000E+00;  Omega22(9,0,3) =  1.0000000E+00;
    Omega22(9,1,0) = -7.0776069E-03;  Omega22(9,1,1) =  1.7917938E-01;  Omega22(9,1,2) = -1.9102410E+00;  Omega22(9,1,3) =  4.6736263E+04;
    Omega22(9,2,0) =  9.8019578E-03;  Omega22(9,2,1) = -1.4699425E-01;  Omega22(9,2,2) =  3.9382460E-02;  Omega22(9,2,3) =  1.5112165E+03;
    Omega22(9,3,0) =  1.1716366E-02;  Omega22(9,3,1) = -1.9289789E-01;  Omega22(9,3,2) =  4.0269474E-01;  Omega22(9,3,3) =  6.0891590E+02;
    Omega22(9,4,0) = -1.4271306E-02;  Omega22(9,4,1) =  3.0401993E-01;  Omega22(9,4,2) = -2.4573879E+00;  Omega22(9,4,3) =  5.3694705E+04;
    Omega22(9,5,0) = -1.7907392E-02;  Omega22(9,5,1) =  4.1207892E-01;  Omega22(9,5,2) = -3.5343610E+00;  Omega22(9,5,3) =  1.4987678E+06;
    Omega22(9,6,0) = -1.0000000E+00;  Omega22(9,6,1) = -1.0000000E+00;  Omega22(9,6,2) = -1.0000000E+00;  Omega22(9,6,3) = -1.0000000E+00;
    Omega22(9,7,0) = -1.0000000E+00;  Omega22(9,7,1) = -1.0000000E+00;  Omega22(9,7,2) = -1.0000000E+00;  Omega22(9,7,3) = -1.0000000E+00;
    Omega22(9,8,0) = -1.0000000E+00;  Omega22(9,8,1) = -1.0000000E+00;  Omega22(9,8,2) = -1.0000000E+00;  Omega22(9,8,3) = -1.0000000E+00;
    Omega22(9,9,0) = -1.0000000E+00;  Omega22(9,9,1) = -1.0000000E+00;  Omega22(9,9,2) = -1.0000000E+00;  Omega22(9,9,3) = -1.0000000E+00;
    Omega22(9,10,0) = -1.0000000E+00; Omega22(9,10,1) = -1.0000000E+00; Omega22(9,10,2) = -1.0000000E+00; Omega22(9,10,3) = -1.0000000E+00;
    // O+(10)
    Omega22(10,0,0) =  1.0000000E+00;  Omega22(10,0,1) =  1.0000000E+00;  Omega22(10,0,2) =  1.0000000E+00;  Omega22(10,0,3) =  1.0000000E+00;
    Omega22(10,1,0) =  1.8733000E-02;  Omega22(10,1,1) = -3.6163781E-01;  Omega22(10,1,2) =  1.6947101E+00;  Omega22(10,1,3) =  2.5244859E+01;
    Omega22(10,2,0) =  1.4207970E-02;  Omega22(10,2,1) = -2.4736726E-01;  Omega22(10,2,2) =  7.4561859E-01;  Omega22(10,2,3) =  3.2519188E+02;
    Omega22(10,3,0) =  1.8015337E-02;  Omega22(10,3,1) = -3.4415293E-01;  Omega22(10,3,2) =  1.5658151E+00;  Omega22(10,3,3) =  3.3303758E+01;
    Omega22(10,4,0) = -2.1681211E-02;  Omega22(10,4,1) =  5.2300453E-01;  Omega22(10,4,2) = -4.5118623E+00;  Omega22(10,4,3) =  2.3467766E+07;
    Omega22(10,5,0) = -1.6032919E-02;  Omega22(10,5,1) =  3.7114396E-01;  Omega22(10,5,2) = -3.2050078E+00;  Omega22(10,5,3) =  5.8099314E+05;
    Omega22(10,6,0) = -1.0000000E+00;  Omega22(10,6,1) = -1.0000000E+00;  Omega22(10,6,2) = -1.0000000E+00;  Omega22(10,6,3) = -1.0000000E+00;
    Omega22(10,7,0) = -1.0000000E+00;  Omega22(10,7,1) = -1.0000000E+00;  Omega22(10,7,2) = -1.0000000E+00;  Omega22(10,7,3) = -1.0000000E+00;
    Omega22(10,8,0) = -1.0000000E+00;  Omega22(10,8,1) = -1.0000000E+00;  Omega22(10,8,2) = -1.0000000E+00;  Omega22(10,8,3) = -1.0000000E+00;
    Omega22(10,9,0) = -1.0000000E+00;  Omega22(10,9,1) = -1.0000000E+00;  Omega22(10,9,2) = -1.0000000E+00;  Omega22(10,9,3) = -1.0000000E+00;
    Omega22(10,10,0) = -1.0000000E+00; Omega22(10,10,1) = -1.0000000E+00; Omega22(10,10,2) = -1.0000000E+00; Omega22(10,10,3) = -1.0000000E+00;




    // Creation/Destruction (+1/-1), Index of monoatomic reactants
    // Monoatomic species (N,O) recombine into diaatomic (N2, O2)
    CatRecombTable(0,0) =  0; CatRecombTable(0,1) = 1;
    CatRecombTable(1,0) =  1; CatRecombTable(1,1) = 4;
    CatRecombTable(2,0) =  1; CatRecombTable(2,1) = 5;
    CatRecombTable(3,0) =  0; CatRecombTable(3,1) = 1;
    CatRecombTable(4,0) = -1; CatRecombTable(4,1) = 4;
    CatRecombTable(5,0) = -1; CatRecombTable(5,1) = 5;
    CatRecombTable(6,0) =  0; CatRecombTable(6,1) = 1;
    CatRecombTable(7,0) =  0; CatRecombTable(7,1) = 1;
    CatRecombTable(8,0) =  0; CatRecombTable(8,1) = 1;
    CatRecombTable(9,0) =  0; CatRecombTable(9,1) = 1;
    CatRecombTable(10,0) =  0; CatRecombTable(10,1) = 1;


    /*--- Values for Sutherland's formula. ---*/
    if (viscous) {
      //F.M. White, Viscous Fluid Flow, 3rd ed., McGraw-Hill, 2006.
      k_ref[0] = 0.0241;
      mu_ref[0] = 1.716E-5;
      Sm_ref[0] = 111.0;
      Sk_ref[0] = 194.0;
    }
  }

  if (ionization) { nHeavy = nSpecies-1; nEl = 1; }
  else            { nHeavy = nSpecies;   nEl = 0; }
}

CSU2TCLib::~CSU2TCLib()= default;

void CSU2TCLib::SetTDStateRhosTTv(vector<su2double>& val_rhos, su2double val_temperature, su2double val_temperature_ve){

  rhos = val_rhos;
  T    = val_temperature;
  Tve  = val_temperature_ve;

  Density = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    Density += rhos[iSpecies];

  Pressure = ComputePressure();

}

vector<su2double>& CSU2TCLib::GetSpeciesCvTraRot(){

  //if(ionization) Cvtrs[0] = 0.0;

  //for (iSpecies = nEl; iSpecies < nHeavy; iSpecies++)
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    Cvtrs[iSpecies] = (3.0/2.0 + RotationModes[iSpecies]/2.0) * Ru/MolarMass[iSpecies];

  return Cvtrs;
}

vector<su2double>& CSU2TCLib::ComputeSpeciesCvVibEle(su2double val_T){

  su2double thoTve, exptv, num, num2, num3, denom, Cvvs, Cves;
  unsigned short iElectron = 0;

  /*--- Loop through species ---*/
  for(iSpecies = 0; iSpecies < nSpecies; iSpecies++){

    /*--- If requesting electron specific heat ---*/
    if (ionization && iSpecies == iElectron) {
      Cvvs = 0.0;
      Cves = 3.0/2.0 * Ru/MolarMass[iSpecies];
    }

    /*--- Heavy particle specific heat ---*/
    else {

      /*--- Vibrational energy ---*/
      if (CharVibTemp[iSpecies] != 0.0) {
        thoTve = CharVibTemp[iSpecies]/val_T;
        exptv = exp(CharVibTemp[iSpecies]/val_T);
        Cvvs  = Ru/MolarMass[iSpecies] * thoTve*thoTve * exptv / ((exptv-1.0)*(exptv-1.0));
      } else {
        Cvvs = 0.0;
      }

      /*--- Electronic energy ---*/
      if (nElStates[iSpecies] != 0) {
        num = 0.0; num2 = 0.0;
        denom = ElDegeneracy[iSpecies][0] * exp(-CharElTemp[iSpecies][0]/val_T);
        num3  = ElDegeneracy[iSpecies][0] * (CharElTemp[iSpecies][0]/(val_T*val_T))*exp(-CharElTemp[iSpecies][0]/val_T);
        for (iEl = 1; iEl < nElStates[iSpecies]; iEl++) {
          thoTve = CharElTemp[iSpecies][iEl]/val_T;
          exptv = exp(-CharElTemp[iSpecies][iEl]/val_T);

          num   += ElDegeneracy[iSpecies][iEl] * CharElTemp[iSpecies][iEl] * exptv;
          denom += ElDegeneracy[iSpecies][iEl] * exptv;
          num2  += ElDegeneracy[iSpecies][iEl] * (thoTve*thoTve) * exptv;
          num3  += ElDegeneracy[iSpecies][iEl] * thoTve/val_T * exptv;
        }
        Cves = Ru/MolarMass[iSpecies] * (num2/denom - num*num3/(denom*denom));
      } else {
        Cves = 0.0;
      }
    }

    Cvves[iSpecies] = Cvvs + Cves;
  }

  return Cvves;

}

vector<su2double>& CSU2TCLib::ComputeMixtureEnergies(){

  su2double Ev, Ee, Ef, num;

  su2double rhoEmix = 0.0;
  su2double rhoEve  = 0.0;
  su2double denom   = 0.0;

  // Electrons
  for (iSpecies = 0; iSpecies < nEl; iSpecies++) {

    // Species formation energy
    Ef = Enthalpy_Formation[iSpecies] - Ru/MolarMass[iSpecies] * Ref_Temperature[iSpecies];

    // Electron t-r mode contributes to mixture vib-el energy
    rhoEve += rhos[iSpecies]*((3.0/2.0) * Ru/MolarMass[iSpecies] * (Tve - Ref_Temperature[iSpecies]));
  }

  for (iSpecies = nEl; iSpecies < nSpecies; iSpecies++){

    // Species formation energy
    Ef = Enthalpy_Formation[iSpecies] - Ru/MolarMass[iSpecies]*Ref_Temperature[iSpecies];

    // Species vibrational energy
    if (CharVibTemp[iSpecies] != 0.0)
      Ev = Ru/MolarMass[iSpecies] * CharVibTemp[iSpecies] / (exp(CharVibTemp[iSpecies]/Tve)-1.0);
    else
      Ev = 0.0;

    // Species electronic energy
    num = 0.0;
    denom = ElDegeneracy(iSpecies,0) * exp(-CharElTemp(iSpecies,0)/Tve);
    for (iEl = 1; iEl < nElStates[iSpecies]; iEl++) {
      num   += ElDegeneracy(iSpecies,iEl) * CharElTemp(iSpecies,iEl) * exp(-CharElTemp(iSpecies,iEl)/Tve);
      denom += ElDegeneracy(iSpecies,iEl) * exp(-CharElTemp(iSpecies,iEl)/Tve);
    }
    Ee = Ru/MolarMass[iSpecies] * (num/denom);

    // Mixture total energy
    rhoEmix += rhos[iSpecies] * ((3.0/2.0+RotationModes[iSpecies]/2.0) * Ru/MolarMass[iSpecies] * (T-Ref_Temperature[iSpecies]) + Ev + Ee + Ef);

    // Mixture vibrational-electronic energy
    rhoEve += rhos[iSpecies] * (Ev + Ee);

  }

  energies[0] = rhoEmix/Density;
  energies[1] = rhoEve/Density;

  return energies;

}

vector<su2double>& CSU2TCLib::ComputeSpeciesEve(su2double val_T, bool vibe_only){

  su2double Ev, Eel, Ef, num, denom;
  unsigned short iElectron = 0;

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++){

    /*--- Electron species energy ---*/
    if ( ionization && (iSpecies == iElectron)) {

      /*--- Calculate formation energy ---*/
      Ef = Enthalpy_Formation[iSpecies] - Ru/MolarMass[iSpecies] * Ref_Temperature[iSpecies];

      /*--- Electron t-r mode contributes to mixture vib-el energy ---*/
      Eel = (3.0/2.0) * Ru/MolarMass[iSpecies] * (val_T - Ref_Temperature[iSpecies]) + Ef;
      Ev  = 0.0;
    }
    /*--- Heavy particle energy ---*/
    else {

      /*--- Calculate vibrational energy (harmonic-oscillator model) ---*/
      if (CharVibTemp[iSpecies] != 0.0)
        Ev = Ru/MolarMass[iSpecies] * CharVibTemp[iSpecies] / (exp(CharVibTemp[iSpecies]/val_T)-1.0);
      else
        Ev = 0.0;

      /*--- Calculate electronic energy ---*/
      num = 0.0;
      denom = ElDegeneracy[iSpecies][0] * exp(-CharElTemp[iSpecies][0]/val_T);
      for (iEl = 1; iEl < nElStates[iSpecies]; iEl++) {
        num   += ElDegeneracy[iSpecies][iEl] * CharElTemp[iSpecies][iEl] * exp(-CharElTemp[iSpecies][iEl]/val_T);
        denom += ElDegeneracy[iSpecies][iEl] * exp(-CharElTemp[iSpecies][iEl]/val_T);
      }
      Eel = Ru/MolarMass[iSpecies] * (num/denom);
    }
    if (vibe_only) {eves[iSpecies] = Ev;}
    else {eves[iSpecies] = Ev + Eel;}
  }

  return eves;
}

vector<su2double>& CSU2TCLib::ComputeNetProductionRates(bool implicit, const su2double *V, const su2double* eve,
                                                        const su2double* cvve, const su2double* dTdU, const su2double* dTvedU,
                                                        su2double **val_jacobian){

  /*---                          ---*/
  /*--- Nonequilibrium chemistry ---*/
  /*---                          ---*/

  /*--- Initialize variables ---*/
  unsigned short ii, iReaction;
  ws.resize(nSpecies,0.0);
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies ++)
    ws[iSpecies] = 0.0;

  /*--- Define artificial chemistry parameters ---*/
  // Note: These parameters artificially increase the rate-controlling reaction
  //       temperature.  This relaxes some of the stiffness in the chemistry
  //       source term.
  const su2double T_min   = 800.0;
  const su2double epsilon = 80;

  /*--- Define preferential dissociation coefficient ---*/
  //alpha = 0.3; //TODO: make this a config option?

  /*--- Loop over all reactions ---*/
  for (iReaction = 0; iReaction < nReactions; iReaction++) {

    /*--- Determine the rate-controlling temperature ---*/
    af = Tcf_a[iReaction];
    bf = Tcf_b[iReaction];
    ab = Tcb_a[iReaction];
    bb = Tcb_b[iReaction];
    Trxnf = pow(T, af)*pow(Tve, bf);
    Trxnb = pow(T, ab)*pow(Tve, bb);

    /*--- Calculate the modified temperature ---*/
    Thf = 0.5 * (Trxnf+T_min + sqrt((Trxnf-T_min)*(Trxnf-T_min)+epsilon*epsilon));
    Thb = 0.5 * (Trxnb+T_min + sqrt((Trxnb-T_min)*(Trxnb-T_min)+epsilon*epsilon));

    /*--- Get the Keq & Arrhenius coefficients ---*/
    ComputeKeqConstants(iReaction);

    /*--- Calculate Keq ---*/
    const su2double Keq = exp(  A[0]*(Thb/1E4) + A[1] + A[2]*log(1E4/Thb)
        + A[3]*(1E4/Thb) + A[4]*(1E4/Thb)*(1E4/Thb) );

    /*--- Calculate rate coefficients ---*/
    kf  = ArrheniusCoefficient[iReaction] * exp(ArrheniusEta[iReaction]*log(Thf)) * exp(-ArrheniusTheta[iReaction]/Thf);
    kfb = ArrheniusCoefficient[iReaction] * exp(ArrheniusEta[iReaction]*log(Thb)) * exp(-ArrheniusTheta[iReaction]/Thb);
    kb  = kfb / Keq;

    /*--- Determine production & destruction of each species ---*/
    fwdRxn = 1.0;
    bkwRxn = 1.0;
    for (ii = 0; ii < 3; ii++) {

      /*--- Reactants ---*/
      iSpecies = Reactions(iReaction,0,ii);
      if ( iSpecies != nSpecies)
        fwdRxn *= 0.001*rhos[iSpecies]/MolarMass[iSpecies];

      /*--- Products ---*/
      jSpecies = Reactions(iReaction,1,ii);
      if (jSpecies != nSpecies) {
        bkwRxn *= 0.001*rhos[jSpecies]/MolarMass[jSpecies];
      }
    }

    fwdRxn = 1000.0 * kf * fwdRxn;
    bkwRxn = 1000.0 * kb * bkwRxn;

    for (ii = 0; ii < 3; ii++) {

      /*--- Products ---*/
      iSpecies = Reactions(iReaction,1,ii);
      if (iSpecies != nSpecies)
        ws[iSpecies] += MolarMass[iSpecies] * (fwdRxn-bkwRxn);

      /*--- Reactants ---*/
      iSpecies = Reactions(iReaction,0,ii);
      if (iSpecies != nSpecies)
        ws[iSpecies] -= MolarMass[iSpecies] * (fwdRxn-bkwRxn);
    }

    if (implicit) {
      ChemistryJacobian(iReaction, V, eve, cvve, dTdU, dTvedU, val_jacobian);
    }
  } //iReaction

  return ws;
}

void CSU2TCLib::ChemistryJacobian(unsigned short iReaction, const su2double *V,
                                  const su2double* eve, const su2double *cvve,
                                  const su2double* dTdU, const su2double* dTvedU,
                                  su2double **val_jacobian) {

  unsigned short ii, iVar, jVar, iSpecies;
  unsigned short nEve = nSpecies+nDim+1;
  unsigned short nVar = nSpecies+nDim+2;

  su2double T_min   = 800.0;
  su2double epsilon = 80;

  /*--- Initializing derivative variables ---*/
  dkf.resize(nVar,0.0);      dkb.resize(nVar,0.0);
  dRfok.resize(nVar,0.0);    dRbok.resize(nVar,0.0);
  alphak.resize(nSpecies,0); betak.resize(nSpecies,0);

  for (iVar=0;iVar<nVar;iVar++){
   dkf[iVar]=0.0; dRfok[iVar]=0.0;
   dkb[iVar]=0.0; dRbok[iVar]=0.0;
  }

  for (iSpecies=0;iSpecies<nSpecies;iSpecies++){
   alphak[iSpecies]=0; betak[iSpecies]=0;
  }

  /*--- Extract additional Arrhenius information ---*/
  su2double eta   = ArrheniusEta[iReaction];
  su2double theta = ArrheniusTheta[iReaction];

  /*--- Derivative of modified temperature wrt Trxnf ---*/
  dThf = 0.5 * (1.0 + (Trxnf-T_min)/sqrt((Trxnf-T_min)*(Trxnf-T_min)
                                                  + epsilon*epsilon));
  dThb = 0.5 * (1.0 + (Trxnb-T_min)/sqrt((Trxnb-T_min)*(Trxnb-T_min)
                                                  + epsilon*epsilon));

  /*--- Fwd rate coefficient derivatives ---*/
  coeff = kf * (eta/Thf+theta/(Thf*Thf)) * dThf;
  for (iVar = 0; iVar < nVar; iVar++) {
    dkf[iVar] = coeff * ( af*Trxnf/T*dTdU[iVar] +
                          bf*Trxnf/Tve*dTvedU[iVar] );
  }

  /*--- Bkwd rate coefficient derivatives ---*/
  coeff = kb * (eta/Thb+theta/(Thb*Thb)) * dThb;
  for (iVar = 0; iVar < nVar; iVar++) {
    dkb[iVar] = coeff*( ab*Trxnb/T*dTdU[iVar] + bb*Trxnb/Tve*dTvedU[iVar])
                      - kb*((A[0]*Thb/1E4 - A[2] - A[3]*1E4/Thb
                      - 2*A[4]*(1E4/Thb)*(1E4/Thb))/Thb) * dThb *
                      ( ab*Trxnb/T*dTdU[iVar] + bb*Trxnb/Tve*dTvedU[iVar]);
  }

  /*--- Rxn rate derivatives ---*/
  for (ii = 0; ii < 3; ii++) {

    /*--- Products ---*/
    iSpecies = Reactions(iReaction,1,ii);
    if (iSpecies != nSpecies)
      betak[iSpecies]++;

    /*--- Reactants ---*/
    iSpecies = Reactions(iReaction,0,ii);
    if (iSpecies != nSpecies)
      alphak[iSpecies]++;
  }

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {

    // Fwd
    dRfok[iSpecies] =  0.001*alphak[iSpecies]/MolarMass[iSpecies] *
                       pow(0.001*rhos[iSpecies]/MolarMass[iSpecies],
                       max(0, alphak[iSpecies]-1));

    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++)
      if (jSpecies != iSpecies)
        dRfok[iSpecies] *= pow(0.001*rhos[jSpecies]/MolarMass[jSpecies],
                                                       alphak[jSpecies]);
    dRfok[iSpecies] *= 1000.0;

    // Bkw
    dRbok[iSpecies] =  0.001*betak[iSpecies]/MolarMass[iSpecies] *
                       pow(0.001*rhos[iSpecies]/MolarMass[iSpecies],
                       max(0, betak[iSpecies]-1));

    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++)
      if (jSpecies != iSpecies)
        dRbok[iSpecies] *= pow(0.001*rhos[jSpecies]/MolarMass[jSpecies],
                                                        betak[jSpecies]);
    dRbok[iSpecies] *= 1000.0;
  }

  for (ii = 0; ii < 3; ii++) {

    /*--- Products ---*/
    iSpecies = Reactions(iReaction,1,ii);
    if (iSpecies != nSpecies) {
      for (iVar = 0; iVar < nVar; iVar++) {
        val_jacobian[iSpecies][iVar] += MolarMass[iSpecies] * ( dkf[iVar]*(fwdRxn/kf) + kf*dRfok[iVar] -
                                                                dkb[iVar]*(bkwRxn/kb) - kb*dRbok[iVar]); //TODO * Volume;
        val_jacobian[nEve][iVar]     += MolarMass[iSpecies] * ( dkf[iVar]*(fwdRxn/kf) + kf*dRfok[iVar] -
                                                                dkb[iVar]*(bkwRxn/kb) - kb*dRbok[iVar]) *
                                                                eve[iSpecies];//TODO * Volume;
      }

      for (jVar = 0; jVar < nVar; jVar++) {
        val_jacobian[nEve][jVar] += MolarMass[iSpecies] * (fwdRxn-bkwRxn)* cvve[iSpecies] *
                                                                             dTvedU[jVar];//TODO * Volume;
      }
    }

    /*--- Reactants ---*/
    iSpecies = Reactions(iReaction,0,ii);
    if (iSpecies != nSpecies) {
      for (iVar = 0; iVar < nVar; iVar++) {
        val_jacobian[iSpecies][iVar] -= MolarMass[iSpecies] * ( dkf[iVar]*(fwdRxn/kf) + kf*dRfok[iVar] -
                                                                dkb[iVar]*(bkwRxn/kb) - kb*dRbok[iVar]);//TODO * Volume;
        val_jacobian[nEve][iVar] -=     MolarMass[iSpecies] * ( dkf[iVar]*(fwdRxn/kf) + kf*dRfok[iVar] -
                                                                dkb[iVar]*(bkwRxn/kb) - kb*dRbok[iVar]) *
                                                                                         eve[iSpecies];//TODO * Volume;
      }

      for (jVar = 0; jVar < nVar; jVar++) {
        val_jacobian[nEve][jVar] -= MolarMass[iSpecies] * (fwdRxn-bkwRxn) * cvve[iSpecies] *
                                                                              dTvedU[jVar];//TODO * Volume;
      }
    }
  } // ii
}

void CSU2TCLib::ComputeKeqConstants(unsigned short val_Reaction) {

  unsigned short ii;

  /*--- Acquire database constants from CConfig ---*/
  GetChemistryEquilConstants(val_Reaction);

  /*--- Calculate mixture number density ---*/
  su2double N = 0.0;
  for (iSpecies =0 ; iSpecies < nSpecies; iSpecies++) {
    N += rhos[iSpecies]/MolarMass[iSpecies]*AVOGAD_CONSTANT;
  }

  /*--- Convert number density from 1/m^3 to 1/cm^3 for table look-up ---*/
  N = N*(1E-6);

  /*--- Determine table index based on mixture N ---*/
  unsigned short tbl_offset = 14;
  unsigned short pwr        = floor(log10(N));

  /*--- Bound the interpolation to table limit values ---*/
  unsigned short iIndex = int(pwr) - tbl_offset;
  if (iIndex <= 0) {
    for (ii = 0; ii < 5; ii++)
      A[ii] = RxnConstantTable(0,ii);
    return;
  } if (iIndex >= 5) {
    for (ii = 0; ii < 5; ii++)
      A[ii] = RxnConstantTable(5,ii);
    return;
  }

  /*--- Calculate interpolation denominator terms avoiding pow() ---*/
  su2double tmp1 = 1.0;
  su2double tmp2 = 1.0;
  for (ii = 0; ii < pwr; ii++) {
    tmp1 *= 10.0;
    tmp2 *= 10.0;
  }
  tmp2 *= 10.0;

  /*--- Interpolate ---*/
  for (ii = 0; ii < 5; ii++) {
    A[ii] =  (RxnConstantTable(iIndex+1,ii) - RxnConstantTable(iIndex,ii))
        / (tmp2 - tmp1) * (N - tmp1)
        + RxnConstantTable(iIndex,ii);
  }
}

su2double CSU2TCLib::ComputeEveSourceTerm(){

  /*---                                                                    ---*/
  /*--- Trans.-rot. & vibrational energy exchange via inelastic collisions ---*/
  /*---                                                                    ---*/
  // Note: Electronic energy not implemented
  // Note: Landau-Teller formulation
  // Note: Millikan & White relaxation time (requires P in Atm.)
  // Note: Park limiting cross section

  su2activematrix mu;
  vector<su2double> MolarFrac;

  MolarFrac.resize(nSpecies,0.0);
  mu.resize(nSpecies,nSpecies)=su2double(0.0);

  su2double omegaVT = 0.0;
  su2double omegaCV = 0.0;

  /*--- Calculate mole fractions ---*/
  su2double N    = 0.0;
  su2double conc = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    conc += rhos[iSpecies] / MolarMass[iSpecies];
    N    += rhos[iSpecies] / MolarMass[iSpecies] * AVOGAD_CONSTANT;
  }

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    MolarFrac[iSpecies] = (rhos[iSpecies] / MolarMass[iSpecies]) / conc;

  /*--- Compute Eve and Eve* ---*/
  eve_eq = ComputeSpeciesEve(T, true);
  eve    = ComputeSpeciesEve(Tve, true);

  /*--- Loop over species to calculate source term --*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {

    /*--- Millikan & White relaxation time ---*/
    su2double num   = 0.0;
    su2double denom = 0.0;
    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      mu(iSpecies,jSpecies) = MolarMass[iSpecies]*MolarMass[jSpecies] / (MolarMass[iSpecies] + MolarMass[jSpecies]);
      const su2double A_sr   = 1.16 * 1E-3 * sqrt(mu(iSpecies,jSpecies)) * pow(CharVibTemp[iSpecies], 4.0/3.0);
      const su2double B_sr   = 0.015 * pow(mu(iSpecies,jSpecies), 0.25);
      const su2double tau_sr = 101325.0/Pressure * exp(A_sr*(pow(T,-1.0/3.0) - B_sr) - 18.42);

      num   += MolarFrac[jSpecies];
      denom += MolarFrac[jSpecies] / tau_sr;
    }

    const su2double tauMW = num / denom;

    /*--- Park limiting cross section ---*/
    const su2double Cs    = sqrt((8.0*Ru*T)/(PI_NUMBER*MolarMass[iSpecies]));
    const su2double sig_s = 3E-21*(2.5E9)/(T*T);

    const su2double tauP = 1/(sig_s*Cs*N);

    /*--- Species relaxation time ---*/
    taus[iSpecies] = tauMW + tauP;

    /*--- Add species contribution to residual ---*/
    omegaVT += rhos[iSpecies] * (eve_eq[iSpecies] -
                                 eve[iSpecies]) / taus[iSpecies];
  }

  /*--- Vibrational energy change due to chemical reactions ---*/
  if(!frozen){
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      omegaCV += ws[iSpecies]*eve[iSpecies];
  }

  omega = omegaVT + omegaCV;

  return omega;

}

void CSU2TCLib::GetEveSourceTermJacobian(const su2double *V, const su2double *eve, const su2double *cvve, const su2double *dTdU, const su2double* dTvedU, su2double **val_jacobian){

  unsigned short iVar;
  unsigned short nEv  = nSpecies+nDim+1;
  unsigned short nVar = nSpecies+nDim+2;

  /*--- Compute Cvvs ---*/
  const auto& cvve_eq = ComputeSpeciesCvVibEle(T);

  /*--- Loop through species ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++){

    for (iVar = 0; iVar < nVar; iVar++) {
        val_jacobian[nEv][iVar] += rhos[iSpecies]/taus[iSpecies]*(cvve_eq[iSpecies]*dTdU[iVar]-cvve[iSpecies]*dTvedU[iVar]);//TODO*Volume;
    }
  }

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      val_jacobian[nEv][iSpecies] += (eve_eq[iSpecies]-eve[iSpecies])/taus[iSpecies];//TODO *Volume;
}

vector<su2double>& CSU2TCLib::ComputeSpeciesEnthalpy(su2double val_T, su2double val_Tve, su2double *val_eves){

  vector<su2double> cvtrs;
  //TODO: ADD Electrons?
  cvtrs = GetSpeciesCvTraRot();

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++){
    eves[iSpecies] = val_eves[iSpecies];
    hs[iSpecies] = Ru/MolarMass[iSpecies]*val_T + cvtrs[iSpecies]*val_T + Enthalpy_Formation[iSpecies] + eves[iSpecies];
  }

  return hs;

}

vector<su2double>& CSU2TCLib::GetDiffusionCoeff(){

  if(Kind_TransCoeffModel == TRANSCOEFFMODEL::WILKE)
   DiffusionCoeffWBE();
  if(Kind_TransCoeffModel == TRANSCOEFFMODEL::GUPTAYOS)
   DiffusionCoeffGY();
  if(Kind_TransCoeffModel == TRANSCOEFFMODEL::SUTHERLAND)
   DiffusionCoeffWBE();

  return DiffusionCoeff;

}

su2double CSU2TCLib::GetViscosity(){

  if(Kind_TransCoeffModel == TRANSCOEFFMODEL::WILKE)
    ViscosityWBE();
  if(Kind_TransCoeffModel == TRANSCOEFFMODEL::GUPTAYOS)
    ViscosityGY();
  if(Kind_TransCoeffModel == TRANSCOEFFMODEL::SUTHERLAND)
    ViscositySuth();

  return Mu;

}

vector<su2double>& CSU2TCLib::GetThermalConductivities(){

  if(Kind_TransCoeffModel == TRANSCOEFFMODEL::WILKE)
    ThermalConductivitiesWBE();
  if(Kind_TransCoeffModel == TRANSCOEFFMODEL::GUPTAYOS)
    ThermalConductivitiesGY();
  if(Kind_TransCoeffModel == TRANSCOEFFMODEL::SUTHERLAND)
    ThermalConductivitiesSuth();

  return ThermalConductivities;

}

void CSU2TCLib::DiffusionCoeffWBE(){

  /*--- Calculate species mole fraction ---*/
  su2double conc = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    MolarFracWBE[iSpecies] = rhos[iSpecies]/MolarMass[iSpecies];
    conc += MolarFracWBE[iSpecies];
  }
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    MolarFracWBE[iSpecies] = MolarFracWBE[iSpecies]/conc;

  /*--- Calculate mixture molar mass (kg/mol) ---*/
  // Note: Species molar masses stored as kg/kmol, need 1E-3 conversion
  su2double M = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    M += MolarMass[iSpecies]*MolarFracWBE[iSpecies];
  M = M*1E-3;

  /*---+++                  +++---*/
  /*--- Diffusion coefficients ---*/
  /*---+++                  +++---*/
  /*--- Solve for binary diffusion coefficients ---*/
  // Note: Dij = Dji, so only loop through req'd indices
  // Note: Correlation requires kg/mol, hence 1E-3 conversion from kg/kmol
  su2activematrix Dij;
  Dij.resize(nSpecies, nSpecies) = su2double(0.0);

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    const su2double Mi = MolarMass[iSpecies]*1E-3;
    for (jSpecies = iSpecies; jSpecies < nSpecies; jSpecies++) {
      const su2double Mj = MolarMass[jSpecies]*1E-3;

      /*--- Calculate the Omega^(1,1)_ij collision cross section ---*/

      /*--- If collisions between electrons/ion, used Coloumb potentials ---*/
      bool coulomb = false;
      if (abs(Omega11(iSpecies, jSpecies, 0)) == 1.0 && ionization) coulomb = true;

      // Used Tve for electron collisions
      const su2double T_col = (iSpecies == 0 && ionization) ? Tve : T;

      /*--- Compute the collisional cross section (omega_ij) ---*/
      const su2double Omega_ij = ComputeCollisionCrossSection(iSpecies, jSpecies, T_col, true, coulomb) / PI_NUMBER;

      /*--- Calculate and populate diffusion coefficients ---*/
      Dij(iSpecies,jSpecies) = 7.1613E-25*M*sqrt(T*(1/Mi+1/Mj))/(Density*Omega_ij);
      Dij(jSpecies,iSpecies) = 7.1613E-25*M*sqrt(T*(1/Mi+1/Mj))/(Density*Omega_ij);
    }
  }

  /*--- Calculate species-mixture diffusion coefficient --*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    DiffusionCoeff[iSpecies] = 0.0;
    su2double denom = 0.0;
    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      if (jSpecies != iSpecies) {
        denom += MolarFracWBE[jSpecies]/Dij(iSpecies,jSpecies);
      }
    }

    if (nSpecies==1) DiffusionCoeff[0] = 0;
    else DiffusionCoeff[iSpecies] = (1-MolarFracWBE[iSpecies])/denom;
  }
}

void CSU2TCLib::ViscosityWBE(){

  /*--- Calculate species mole fraction ---*/
  su2double conc = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    MolarFracWBE[iSpecies] = rhos[iSpecies]/MolarMass[iSpecies];
    conc += MolarFracWBE[iSpecies];
  }
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    MolarFracWBE[iSpecies] = MolarFracWBE[iSpecies]/conc;

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    mus[iSpecies] = 0.1*exp((Blottner[iSpecies][0]*log(T)  +
                             Blottner[iSpecies][1])*log(T) +
                             Blottner[iSpecies][2]);

  /*--- Determine species 'phi' value for Blottner model ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    phis[iSpecies] = 0.0;
    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      const su2double tmp1 = 1.0 + sqrt(mus[iSpecies]/mus[jSpecies])*pow(MolarMass[jSpecies]/MolarMass[iSpecies], 0.25);
      const su2double tmp2 = sqrt(8.0*(1.0+MolarMass[iSpecies]/MolarMass[jSpecies]));
      phis[iSpecies] += MolarFracWBE[jSpecies]*tmp1*tmp1/tmp2;
    }
  }

  /*--- Calculate mixture laminar viscosity ---*/
  Mu = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++){
    Mu += MolarFracWBE[iSpecies]*mus[iSpecies]/phis[iSpecies];
  }
}

void CSU2TCLib::ThermalConductivitiesWBE(){

  vector<su2double> ks, kves;

  ks.resize(nSpecies,0.0);
  kves.resize(nSpecies,0.0);

  Cvves = ComputeSpeciesCvVibEle(Tve);

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    ks[iSpecies] = mus[iSpecies]*(15.0/4.0 + RotationModes[iSpecies]/2.0)*Ru/MolarMass[iSpecies];
    kves[iSpecies] = mus[iSpecies]*Cvves[iSpecies];
  }

  /*--- Calculate mixture tr & ve conductivities ---*/
  ThermalCond_tr = 0.0;
  ThermalCond_ve = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    ThermalCond_tr += MolarFracWBE[iSpecies]*ks[iSpecies]/phis[iSpecies];
    ThermalCond_ve += MolarFracWBE[iSpecies]*kves[iSpecies]/phis[iSpecies];
  }

  ThermalConductivities[0] = ThermalCond_tr;
  ThermalConductivities[1] = ThermalCond_ve;
}

su2double CSU2TCLib::ComputeCollisionCrossSection(unsigned iSpecies, unsigned jSpecies, su2double T, bool d1, bool coulomb) {

  const su2double pi = PI_NUMBER;
  const su2double Na = AVOGAD_CONSTANT;

  if (coulomb) {

    const su2double e_cgs = FUND_ELEC_CHARGE_CGS; // CGS unit of fundamental electric charge
    const su2double kb_cgs = BOLTZMANN_CONSTANT * 1E7; // CGS unit of Boltzmann Constant
    const su2double ne_cgs = Na * rhos[0] / MolarMass[0] * 1E-6; // CGS unit of electron number density

    const su2double debyeLength = sqrt(kb_cgs * T / 4 / pi / ne_cgs / pow(e_cgs,2));
    const su2double T_star = debyeLength / (pow(e_cgs,2) / (kb_cgs * T));

    /*--- Compute the collisionion cross section ---*/
    // Note: Omega11 is used for diffusion, viscosity, translational, internal, and reaction components of
    //       thermal conductivity
    //       Omega22 is used for viscosity and translational components of thermal conductivity

    if (Omega11(iSpecies, jSpecies, 0) == 1.0 && d1) {
      return 1E-20 * 5E15 * pi * pow((debyeLength / T), 2) * log(D1_a*T_star*(1 - C1_a * exp(-c1_a * T_star))+1);
    } if (Omega11(iSpecies, jSpecies, 0) == -1.0 && d1) {
      return 1E-20 * 5E15 * pi * pow((debyeLength / T), 2) * log(D1_r*T_star*(1 - C1_r * exp(-c1_r * T_star))+1);
    } else if (Omega22(iSpecies, jSpecies, 0) == 1.0 && !d1) {
      return 1E-20 * 5E15 * pi * pow((debyeLength / T), 2) * log(D2_a*T_star*(1 - C2_a * exp(-c2_a * T_star))+1);
    } else {
      return 1E-20 * 5E15 * pi * pow((debyeLength / T), 2) * log(D2_r*T_star*(1 - C2_r * exp(-c2_r * T_star))+1);
    }

  } else {
    if (d1) {
      return 1E-20 * Omega11(iSpecies,jSpecies,3) * pow(T, Omega11(iSpecies,jSpecies,0)*log(T)*log(T) + Omega11(iSpecies,jSpecies,1)*log(T) + Omega11(iSpecies,jSpecies,2));
    }       return 1E-20 * Omega22(iSpecies,jSpecies,3) * pow(T, Omega22(iSpecies,jSpecies,0)*log(T)*log(T) + Omega22(iSpecies,jSpecies,1)*log(T) + Omega22(iSpecies,jSpecies,2));

  }
}

su2double CSU2TCLib::ComputeCollisionDelta(unsigned iSpecies, unsigned jSpecies, su2double Mi, su2double Mj, su2double T, bool d1) {

  bool coulomb = false;
  if (abs(Omega11(iSpecies, jSpecies, 0)) == 1.0 && ionization) {
    coulomb = true;
  }

  const su2double Omega_ij = ComputeCollisionCrossSection(iSpecies, jSpecies, T, d1, coulomb);
  const su2double pi = PI_NUMBER;
  su2double delta = 0.0;

  if (d1) {
    delta = 8.0/3.0 * sqrt((2.0*Mi*Mj) / (pi*Ru*T*(Mi+Mj))) * Omega_ij; // d1_ij
  } else {
    delta = 16.0/5.0 * sqrt((2.0*Mi*Mj) / (pi*Ru*T*(Mi+Mj))) * Omega_ij; // d2_ij
  }
  return fmin(delta, 1E16);
}

void CSU2TCLib::DiffusionCoeffGY(){

  /*--- Calculate mixture gas constant ---*/
  su2double gam_t = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    gam_t += rhos[iSpecies] / (Density*MolarMass[iSpecies]);
  }

  /*--- Mixture thermal conductivity via Gupta-Yos approximation ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {

    /*--- Initialize the species diffusion coefficient ---*/
    DiffusionCoeff[iSpecies] = 0.0;

    /*--- Calculate molar concentration ---*/
    const su2double Mi    = (MolarMass[iSpecies] + EPS);
    const su2double gam_i = rhos[iSpecies] / (Density*Mi);
    su2double denom = 0.0;

    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      if (jSpecies != iSpecies) {

        const su2double Mj    = (MolarMass[jSpecies] + EPS);
        const su2double gam_j = rhos[jSpecies] / (Density*Mj);

        const su2double kb = BOLTZMANN_CONSTANT;

        const su2double T_col = (iSpecies == 0 && ionization) ? Tve : T;

        su2double d1_ij = ComputeCollisionDelta(iSpecies, jSpecies, Mi, Mj, T_col, true);

        const su2double D_ij = kb*T_col/(Pressure*d1_ij);
        denom += gam_j/D_ij;
      }
    }
    /*--- Calculate species diffusion coefficient ---*/
    DiffusionCoeff[iSpecies] = (gam_t*gam_t*Mi*(1-Mi*gam_i) / denom);
  }
}

void CSU2TCLib::ViscosityGY(){

  const su2double Na = AVOGAD_CONSTANT;
  Mu = 0.0;

  /*--- Mixture viscosity via Gupta-Yos approximation ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {

    su2double denom = 0.0;

    /*--- Calculate molar concentration ---*/
    const su2double Mi    = (MolarMass[iSpecies] + EPS);
    const su2double gam_i = rhos[iSpecies] / (Density*Mi);

    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      const su2double Mj    = (MolarMass[jSpecies] + EPS);
      const su2double gam_j = rhos[jSpecies] / (Density*Mj);

      const su2double T_col = (iSpecies == 0 && ionization) ? Tve : T;

      su2double d2_ij = ComputeCollisionDelta(iSpecies, jSpecies, Mi, Mj, T_col, false);

      denom += gam_j*d2_ij;
    }
    /*--- Calculate species laminar viscosity ---*/
    Mu += (Mi/Na * gam_i) / denom;
  }
}

void CSU2TCLib::ThermalConductivitiesGY(){

  const su2double Na   = AVOGAD_CONSTANT;
  const su2double kb   = BOLTZMANN_CONSTANT;

  /*--- Mixture vibrational-electronic specific heat ---*/
  const auto Cvves = ComputeSpeciesCvVibEle(Tve);
  su2double rhoCvve = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    rhoCvve += rhos[iSpecies]*Cvves[iSpecies];
  const su2double Cvve = rhoCvve/Density;

  /*--- Calculate mixture gas constant ---*/
  su2double R = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    R += Ru / MolarMass[iSpecies] * rhos[iSpecies]/Density;
  }

  /*--- Mixture thermal conductivity via Gupta-Yos approximation ---*/
  su2double ThermalCond_tr = 0.0;
  su2double ThermalCond_ve = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {

    /*--- Calculate molar concentration ---*/
    const su2double Mi    = (MolarMass[iSpecies] + EPS);
    const su2double mi    = Mi/Na;
    const su2double gam_i = rhos[iSpecies] / (Density*Mi);
    su2double denom_t = 0.0;
    su2double denom_r = 0.0;
    su2double denom_re = 0.0;

    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      const su2double Mj    = (MolarMass[jSpecies] + EPS);
      const su2double mj    = Mj/Na;
      const su2double gam_j = rhos[iSpecies] / (Density*Mj);
      const su2double a_ij = 1.0 + (1.0 - mi/mj)*(0.45 - 2.54*mi/mj) / ((1.0 + mi/mj)*(1.0 + mi/mj));

      const su2double T_col = ((iSpecies == 0 && ionization) || (jSpecies == 0 && ionization)) ? Tve : T;

      su2double d1_ij = ComputeCollisionDelta(iSpecies, jSpecies, Mi, Mj, T_col, true);
      su2double d2_ij = ComputeCollisionDelta(iSpecies, jSpecies, Mi, Mj, T_col, false);

      if (jSpecies == 0 && ionization) { denom_t += 3.54*gam_j*d2_ij; }
      else { denom_t += a_ij*gam_j*d2_ij; }

      denom_r += gam_j*d1_ij;
      denom_re += gam_j*d2_ij;
    }

    /*--- Prevent divide by 0 ---*/
    if (denom_t <= 0.0) denom_t = EPS;
    if (denom_r <= 0.0) denom_r = EPS;
    if (denom_re <= 0.0) denom_re = EPS;

    /*--- Translational contribution to thermal conductivity ---*/
    if (!ionization || iSpecies != 0) ThermalCond_tr += ((15.0/4.0)*kb*gam_i/denom_t);

    /*--- Rotational contribution to thermal conductivity ---*/
    if (RotationModes[iSpecies] != 0.0) ThermalCond_tr += (kb*gam_i/denom_r);

    /*--- Vibrational-electronic contribution to thermal conductivity ---*/
    if ((!ionization || iSpecies != 0) && RotationModes[iSpecies] != 0.0) ThermalCond_ve += (kb*Cvve/R*gam_i / denom_r);

    if (ionization && iSpecies == 0) ThermalCond_ve += ((15.0/4.0)*kb*gam_i/(1.45*denom_re));
  }

  ThermalConductivities[0] = ThermalCond_tr;
  ThermalConductivities[1] = ThermalCond_ve;
}

void CSU2TCLib::ViscositySuth(){

  su2double T_nd = T / T_ref_suth;

  /*--- Calculate mixture laminar viscosity ---*/
  Mu = mu_ref[0] * T_nd * sqrt(T_nd) * ((T_ref_suth + Sm_ref[0]) / (T + Sm_ref[0]));
}

void CSU2TCLib::ThermalConductivitiesSuth(){

  /*--- Compute mixture quantities ---*/
  su2double mass = 0.0, rho = 0.0;
  for (unsigned short ii=0; ii<nSpecies; ii++) rho  += rhos[ii];
  for (unsigned short ii=0; ii<nSpecies; ii++) mass += rhos[ii]/rho*MolarMass[ii];

  su2double Cvtr = ComputerhoCvtr()/rho;
  su2double Cvve = ComputerhoCvve()/rho;

  /*--- Compute simple Kve scaling factor ---*/
  su2double scl  = Cvve/Cvtr;

  /*--- Compute k's using Sutherland's law ---*/
  su2double T_nd = T / T_ref_suth;
  su2double k = k_ref[0] * T_nd * sqrt(T_nd) * ((T_ref_suth + Sk_ref[0]) / (T + Sk_ref[0]));
  su2double kve = scl*k;

  ThermalConductivities[0] = k;
  ThermalConductivities[1] = kve;
}

vector<su2double>& CSU2TCLib::ComputeTemperatures(vector<su2double>& val_rhos, su2double rhoE, su2double rhoEve, su2double rhoEvel, su2double Tve_old) {

  rhos = val_rhos;

  /*----------Translational temperature----------*/
  su2double rhoE_f   = 0.0;
  su2double rhoE_ref = 0.0;
  su2double rhoCvtr  = 0.0;
  for (iSpecies = nEl; iSpecies < nSpecies; iSpecies++) {
    rhoCvtr  += rhos[iSpecies] * Cvtrs[iSpecies];
    rhoE_ref += rhos[iSpecies] * Cvtrs[iSpecies] * Ref_Temperature[iSpecies];
    rhoE_f   += rhos[iSpecies] * (Enthalpy_Formation[iSpecies] - Ru/MolarMass[iSpecies]*Ref_Temperature[iSpecies]);
  }

  T = (rhoE - rhoEve - rhoE_f + rhoE_ref - rhoEvel) / rhoCvtr;

  /*--- Set temperature clipping values ---*/
  const su2double Tmin   = 50.0; const su2double Tmax   = 8E4;
  const su2double Tvemin = 50.0; const su2double Tvemax = 8E4;
  su2double Tve_o  = 50.0; su2double Tve2  = 8E4;

  /* Determine if the temperature lies within the acceptable range */
  if (Tve_old < 1) Tve_old = T;                           //For first fluid iteration
  if (T < Tmin) T = Tmin;  else if (T > Tmax) T = Tmax;
  if (Tve_old<Tvemin) Tve_old = Tvemin; else if (Tve_old>Tvemax) Tve_old = Tvemax;

  /*--- Set vibrational temperature algorithm parameters ---*/
  const su2double NRtol         = 1.0E-6;    // Tolerance for the Newton-Raphson method
  const su2double Btol          = 1.0E-6;    // Tolerance for the Bisection method
  const unsigned short maxBIter = 100;        // Maximum Bisection method iterations
  const unsigned short maxNIter = 100;        // Maximum Newton-Raphson iterations
  const su2double scale         = 0.9;       // Scaling factor for Newton-Raphson step

  /*--- Execute a Newton-Raphson root-finding method for Tve ---*/
  //Initialize solution
  Tve = Tve_old;

  bool Bconvg = false;
  bool NRconvg = false;
  su2double rhoEve_t = 0.0, rhoCvve = 0.0;

  /*--- Newton-Raphson Method --*/
  for (unsigned short iIter = 0; iIter < maxNIter; iIter++) {
    rhoEve_t = rhoCvve = 0.0;
    const auto& val_eves  = ComputeSpeciesEve(Tve);
    const auto& val_cvves = ComputeSpeciesCvVibEle(Tve);

    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++){
      rhoEve_t += rhos[iSpecies] * val_eves[iSpecies];
      rhoCvve += rhos[iSpecies] * val_cvves[iSpecies];
    }

    /*--- Find the roots ---*/
    su2double f  = rhoEve - rhoEve_t;
    su2double df = -rhoCvve;
    Tve2 = Tve - (f/df)*scale;

    /*--- Check for convergence ---*/
    if ((fabs(Tve2-Tve) < NRtol) && (Tve > Tvemin) && (Tve < Tvemax)) {
      NRconvg = true;
      Tve = Tve2;
      break;
    }       Tve = Tve2;

  }

  // If the Newton-Raphson method has converged, assign the value of Tve.
  // Otherwise, execute a bisection root-finding method
  Tve_o = Tvemin; Tve2 = Tvemax;
  if (!NRconvg) {
    for (unsigned short iIter = 0; iIter < maxBIter; iIter++) {
      Tve      = (Tve_o+Tve2)/2.0;
      const auto& val_eves = ComputeSpeciesEve(Tve);
      rhoEve_t = 0.0;
      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) rhoEve_t += rhos[iSpecies] * val_eves[iSpecies];
      if (fabs(rhoEve_t - rhoEve) < Btol) {
        Bconvg = true;
        break;
      }         if (rhoEve_t > rhoEve) Tve2 = Tve;
        else                  Tve_o = Tve;

    }
  }

  // If absolutely no convergence, then assign to the TR temperature
  if (!NRconvg && !Bconvg ) {
    Tve = T;
  }

  temperatures[0] = T;
  temperatures[1] = Tve;

  return temperatures;
}

void CSU2TCLib::GetChemistryEquilConstants(unsigned short iReaction){

  if (gas_model == "O2"){
    // THESE ARE UNUSED.  SHOULD WE KEEP????  Good for future?
    //O2 + M -> 2O + M
    RxnConstantTable(0,0) = 1.8103;  RxnConstantTable(0,1) = 1.9607;  RxnConstantTable(0,2) = 3.5716;  RxnConstantTable(0,3) = -7.3623;   RxnConstantTable(0,4) = 0.083861;
    RxnConstantTable(1,0) = 0.91354; RxnConstantTable(1,1) = 2.3160;  RxnConstantTable(1,2) = 2.2885;  RxnConstantTable(1,3) = -6.7969;   RxnConstantTable(1,4) = 0.046338;
    RxnConstantTable(2,0) = 0.64183; RxnConstantTable(2,1) = 2.4253;  RxnConstantTable(2,2) = 1.9026;  RxnConstantTable(2,3) = -6.6277;   RxnConstantTable(2,4) = 0.035151;
    RxnConstantTable(3,0) = 0.55388; RxnConstantTable(3,1) = 2.4600;  RxnConstantTable(3,2) = 1.7763;  RxnConstantTable(3,3) = -6.5720;   RxnConstantTable(3,4) = 0.031445;
    RxnConstantTable(4,0) = 0.52455; RxnConstantTable(4,1) = 2.4715;  RxnConstantTable(4,2) = 1.7342;  RxnConstantTable(4,3) = -6.55534;  RxnConstantTable(4,4) = 0.030209;
    RxnConstantTable(5,0) = 0.50989; RxnConstantTable(5,1) = 2.4773;  RxnConstantTable(5,2) = 1.7132;  RxnConstantTable(5,3) = -6.5441;   RxnConstantTable(5,4) = 0.029591;

  } else if (gas_model == "N2"){

    //N2 + M -> 2N + M
    RxnConstantTable(0,0) = 3.4907;  RxnConstantTable(0,1) = 0.83133; RxnConstantTable(0,2) = 4.0978;  RxnConstantTable(0,3) = -12.728; RxnConstantTable(0,4) = 0.07487;   //n = 1E14
    RxnConstantTable(1,0) = 2.0723;  RxnConstantTable(1,1) = 1.38970; RxnConstantTable(1,2) = 2.0617;  RxnConstantTable(1,3) = -11.828; RxnConstantTable(1,4) = 0.015105;  //n = 1E15
    RxnConstantTable(2,0) = 1.6060;  RxnConstantTable(2,1) = 1.57320; RxnConstantTable(2,2) = 1.3923;  RxnConstantTable(2,3) = -11.533; RxnConstantTable(2,4) = -0.004543; //n = 1E16
    RxnConstantTable(3,0) = 1.5351;  RxnConstantTable(3,1) = 1.60610; RxnConstantTable(3,2) = 1.2993;  RxnConstantTable(3,3) = -11.494; RxnConstantTable(3,4) = -0.00698;  //n = 1E17
    RxnConstantTable(4,0) = 1.4766;  RxnConstantTable(4,1) = 1.62910; RxnConstantTable(4,2) = 1.2153;  RxnConstantTable(4,3) = -11.457; RxnConstantTable(4,4) = -0.00944;  //n = 1E18
    RxnConstantTable(5,0) = 1.4766;  RxnConstantTable(5,1) = 1.62910; RxnConstantTable(5,2) = 1.2153;  RxnConstantTable(5,3) = -11.457; RxnConstantTable(5,4) = -0.00944;  //n = 1E19

  } else if (gas_model == "AIR-5"){

    if (iReaction <= 4) {

      //N2 + M -> 2N + M
      RxnConstantTable(0,0) = 3.4907;  RxnConstantTable(0,1) = 0.83133; RxnConstantTable(0,2) = 4.0978;  RxnConstantTable(0,3) = -12.728; RxnConstantTable(0,4) = 0.07487;   //n = 1E14
      RxnConstantTable(1,0) = 2.0723;  RxnConstantTable(1,1) = 1.38970; RxnConstantTable(1,2) = 2.0617;  RxnConstantTable(1,3) = -11.828; RxnConstantTable(1,4) = 0.015105;  //n = 1E15
      RxnConstantTable(2,0) = 1.6060;  RxnConstantTable(2,1) = 1.57320; RxnConstantTable(2,2) = 1.3923;  RxnConstantTable(2,3) = -11.533; RxnConstantTable(2,4) = -0.004543; //n = 1E16
      RxnConstantTable(3,0) = 1.5351;  RxnConstantTable(3,1) = 1.60610; RxnConstantTable(3,2) = 1.2993;  RxnConstantTable(3,3) = -11.494; RxnConstantTable(3,4) = -0.00698;  //n = 1E17
      RxnConstantTable(4,0) = 1.4766;  RxnConstantTable(4,1) = 1.62910; RxnConstantTable(4,2) = 1.2153;  RxnConstantTable(4,3) = -11.457; RxnConstantTable(4,4) = -0.00944;  //n = 1E18
      RxnConstantTable(5,0) = 1.4766;  RxnConstantTable(5,1) = 1.62910; RxnConstantTable(5,2) = 1.2153;  RxnConstantTable(5,3) = -11.457; RxnConstantTable(5,4) = -0.00944;  //n = 1E19

    } else if (iReaction > 4 && iReaction <= 9) {

      //O2 + M -> 2O + M
      RxnConstantTable(0,0) = 1.8103;  RxnConstantTable(0,1) = 1.9607;  RxnConstantTable(0,2) = 3.5716;  RxnConstantTable(0,3) = -7.3623;   RxnConstantTable(0,4) = 0.083861;
      RxnConstantTable(1,0) = 0.91354; RxnConstantTable(1,1) = 2.3160;  RxnConstantTable(1,2) = 2.2885;  RxnConstantTable(1,3) = -6.7969;   RxnConstantTable(1,4) = 0.046338;
      RxnConstantTable(2,0) = 0.64183; RxnConstantTable(2,1) = 2.4253;  RxnConstantTable(2,2) = 1.9026;  RxnConstantTable(2,3) = -6.6277;   RxnConstantTable(2,4) = 0.035151;
      RxnConstantTable(3,0) = 0.55388; RxnConstantTable(3,1) = 2.4600;  RxnConstantTable(3,2) = 1.7763;  RxnConstantTable(3,3) = -6.5720;   RxnConstantTable(3,4) = 0.031445;
      RxnConstantTable(4,0) = 0.52455; RxnConstantTable(4,1) = 2.4715;  RxnConstantTable(4,2) = 1.7342;  RxnConstantTable(4,3) = -6.55534;  RxnConstantTable(4,4) = 0.030209;
      RxnConstantTable(5,0) = 0.50989; RxnConstantTable(5,1) = 2.4773;  RxnConstantTable(5,2) = 1.7132;  RxnConstantTable(5,3) = -6.5441;   RxnConstantTable(5,4) = 0.029591;

    } else if (iReaction > 9 && iReaction <= 14) {

      //NO + M -> N + O + M
      RxnConstantTable(0,0) = 2.1649;  RxnConstantTable(0,1) = 0.078577;  RxnConstantTable(0,2) = 2.8508;  RxnConstantTable(0,3) = -8.5422; RxnConstantTable(0,4) = 0.053043;
      RxnConstantTable(1,0) = 1.0072;  RxnConstantTable(1,1) = 0.53545;   RxnConstantTable(1,2) = 1.1911;  RxnConstantTable(1,3) = -7.8098; RxnConstantTable(1,4) = 0.004394;
      RxnConstantTable(2,0) = 0.63817; RxnConstantTable(2,1) = 0.68189;   RxnConstantTable(2,2) = 0.66336; RxnConstantTable(2,3) = -7.5773; RxnConstantTable(2,4) = -0.011025;
      RxnConstantTable(3,0) = 0.55889; RxnConstantTable(3,1) = 0.71558;   RxnConstantTable(3,2) = 0.55396; RxnConstantTable(3,3) = -7.5304; RxnConstantTable(3,4) = -0.014089;
      RxnConstantTable(4,0) = 0.5150;  RxnConstantTable(4,1) = 0.73286;   RxnConstantTable(4,2) = 0.49096; RxnConstantTable(4,3) = -7.5025; RxnConstantTable(4,4) = -0.015938;
      RxnConstantTable(5,0) = 0.50765; RxnConstantTable(5,1) = 0.73575;   RxnConstantTable(5,2) = 0.48042; RxnConstantTable(5,3) = -7.4979; RxnConstantTable(5,4) = -0.016247;

    } else if (iReaction == 15) {

      //N2 + O -> NO + N
      RxnConstantTable(0,0) = 1.3261;  RxnConstantTable(0,1) = 0.75268; RxnConstantTable(0,2) = 1.2474;  RxnConstantTable(0,3) = -4.1857; RxnConstantTable(0,4) = 0.02184;
      RxnConstantTable(1,0) = 1.0653;  RxnConstantTable(1,1) = 0.85417; RxnConstantTable(1,2) = 0.87093; RxnConstantTable(1,3) = -4.0188; RxnConstantTable(1,4) = 0.010721;
      RxnConstantTable(2,0) = 0.96794; RxnConstantTable(2,1) = 0.89131; RxnConstantTable(2,2) = 0.7291;  RxnConstantTable(2,3) = -3.9555; RxnConstantTable(2,4) = 0.006488;
      RxnConstantTable(3,0) = 0.97646; RxnConstantTable(3,1) = 0.89043; RxnConstantTable(3,2) = 0.74572; RxnConstantTable(3,3) = -3.9642; RxnConstantTable(3,4) = 0.007123;
      RxnConstantTable(4,0) = 0.96188; RxnConstantTable(4,1) = 0.89617; RxnConstantTable(4,2) = 0.72479; RxnConstantTable(4,3) = -3.955;  RxnConstantTable(4,4) = 0.006509;
      RxnConstantTable(5,0) = 0.96921; RxnConstantTable(5,1) = 0.89329; RxnConstantTable(5,2) = 0.73531; RxnConstantTable(5,3) = -3.9596; RxnConstantTable(5,4) = 0.006818;

    } else if (iReaction == 16) {

      //NO + O -> O2 + N
      RxnConstantTable(0,0) = 0.35438;   RxnConstantTable(0,1) = -1.8821; RxnConstantTable(0,2) = -0.72111;  RxnConstantTable(0,3) = -1.1797;   RxnConstantTable(0,4) = -0.030831;
      RxnConstantTable(1,0) = 0.093613;  RxnConstantTable(1,1) = -1.7806; RxnConstantTable(1,2) = -1.0975;   RxnConstantTable(1,3) = -1.0128;   RxnConstantTable(1,4) = -0.041949;
      RxnConstantTable(2,0) = -0.003732; RxnConstantTable(2,1) = -1.7434; RxnConstantTable(2,2) = -1.2394;   RxnConstantTable(2,3) = -0.94952;  RxnConstantTable(2,4) = -0.046182;
      RxnConstantTable(3,0) = 0.004815;  RxnConstantTable(3,1) = -1.7443; RxnConstantTable(3,2) = -1.2227;   RxnConstantTable(3,3) = -0.95824;  RxnConstantTable(3,4) = -0.045545;
      RxnConstantTable(4,0) = -0.009758; RxnConstantTable(4,1) = -1.7386; RxnConstantTable(4,2) = -1.2436;   RxnConstantTable(4,3) = -0.949;    RxnConstantTable(4,4) = -0.046159;
      RxnConstantTable(5,0) = -0.002428; RxnConstantTable(5,1) = -1.7415; RxnConstantTable(5,2) = -1.2331;   RxnConstantTable(5,3) = -0.95365;  RxnConstantTable(5,4) = -0.04585;
    }

  } else if (gas_model == "AIR-7"){

    if (iReaction <= 5) {

      //N2 + M -> 2N + M
      RxnConstantTable(0,0) = 3.4907;  RxnConstantTable(0,1) = 0.83133; RxnConstantTable(0,2) = 4.0978;  RxnConstantTable(0,3) = -12.728; RxnConstantTable(0,4) = 0.07487;   //n = 1E14
      RxnConstantTable(1,0) = 2.0723;  RxnConstantTable(1,1) = 1.38970; RxnConstantTable(1,2) = 2.0617;  RxnConstantTable(1,3) = -11.828; RxnConstantTable(1,4) = 0.015105;  //n = 1E15
      RxnConstantTable(2,0) = 1.6060;  RxnConstantTable(2,1) = 1.57320; RxnConstantTable(2,2) = 1.3923;  RxnConstantTable(2,3) = -11.533; RxnConstantTable(2,4) = -0.004543; //n = 1E16
      RxnConstantTable(3,0) = 1.5351;  RxnConstantTable(3,1) = 1.60610; RxnConstantTable(3,2) = 1.2993;  RxnConstantTable(3,3) = -11.494; RxnConstantTable(3,4) = -0.00698;  //n = 1E17
      RxnConstantTable(4,0) = 1.4766;  RxnConstantTable(4,1) = 1.62910; RxnConstantTable(4,2) = 1.2153;  RxnConstantTable(4,3) = -11.457; RxnConstantTable(4,4) = -0.00944;  //n = 1E18
      RxnConstantTable(5,0) = 1.4766;  RxnConstantTable(5,1) = 1.62910; RxnConstantTable(5,2) = 1.2153;  RxnConstantTable(5,3) = -11.457; RxnConstantTable(5,4) = -0.00944;  //n = 1E19

    } else if (iReaction > 5 && iReaction <= 11) {

      //O2 + M -> 2O + M
      RxnConstantTable(0,0) = 1.8103;  RxnConstantTable(0,1) = 1.9607;  RxnConstantTable(0,2) = 3.5716;  RxnConstantTable(0,3) = -7.3623;   RxnConstantTable(0,4) = 0.083861;
      RxnConstantTable(1,0) = 0.91354; RxnConstantTable(1,1) = 2.3160;  RxnConstantTable(1,2) = 2.2885;  RxnConstantTable(1,3) = -6.7969;   RxnConstantTable(1,4) = 0.046338;
      RxnConstantTable(2,0) = 0.64183; RxnConstantTable(2,1) = 2.4253;  RxnConstantTable(2,2) = 1.9026;  RxnConstantTable(2,3) = -6.6277;   RxnConstantTable(2,4) = 0.035151;
      RxnConstantTable(3,0) = 0.55388; RxnConstantTable(3,1) = 2.4600;  RxnConstantTable(3,2) = 1.7763;  RxnConstantTable(3,3) = -6.5720;   RxnConstantTable(3,4) = 0.031445;
      RxnConstantTable(4,0) = 0.52455; RxnConstantTable(4,1) = 2.4715;  RxnConstantTable(4,2) = 1.7342;  RxnConstantTable(4,3) = -6.55534;  RxnConstantTable(4,4) = 0.030209;
      RxnConstantTable(5,0) = 0.50989; RxnConstantTable(5,1) = 2.4773;  RxnConstantTable(5,2) = 1.7132;  RxnConstantTable(5,3) = -6.5441;   RxnConstantTable(5,4) = 0.029591;

    } else if (iReaction > 11 && iReaction <= 17) {

      //NO + M -> N + O + M
      RxnConstantTable(0,0) = 2.1649;  RxnConstantTable(0,1) = 0.078577;  RxnConstantTable(0,2) = 2.8508;  RxnConstantTable(0,3) = -8.5422; RxnConstantTable(0,4) = 0.053043;
      RxnConstantTable(1,0) = 1.0072;  RxnConstantTable(1,1) = 0.53545;   RxnConstantTable(1,2) = 1.1911;  RxnConstantTable(1,3) = -7.8098; RxnConstantTable(1,4) = 0.004394;
      RxnConstantTable(2,0) = 0.63817; RxnConstantTable(2,1) = 0.68189;   RxnConstantTable(2,2) = 0.66336; RxnConstantTable(2,3) = -7.5773; RxnConstantTable(2,4) = -0.011025;
      RxnConstantTable(3,0) = 0.55889; RxnConstantTable(3,1) = 0.71558;   RxnConstantTable(3,2) = 0.55396; RxnConstantTable(3,3) = -7.5304; RxnConstantTable(3,4) = -0.014089;
      RxnConstantTable(4,0) = 0.5150;  RxnConstantTable(4,1) = 0.73286;   RxnConstantTable(4,2) = 0.49096; RxnConstantTable(4,3) = -7.5025; RxnConstantTable(4,4) = -0.015938;
      RxnConstantTable(5,0) = 0.50765; RxnConstantTable(5,1) = 0.73575;   RxnConstantTable(5,2) = 0.48042; RxnConstantTable(5,3) = -7.4979; RxnConstantTable(5,4) = -0.016247;

    } else if (iReaction == 18) {

      //N2 + O -> NO + N
      RxnConstantTable(0,0) = 1.3261;  RxnConstantTable(0,1) = 0.75268; RxnConstantTable(0,2) = 1.2474;  RxnConstantTable(0,3) = -4.1857; RxnConstantTable(0,4) = 0.02184;
      RxnConstantTable(1,0) = 1.0653;  RxnConstantTable(1,1) = 0.85417; RxnConstantTable(1,2) = 0.87093; RxnConstantTable(1,3) = -4.0188; RxnConstantTable(1,4) = 0.010721;
      RxnConstantTable(2,0) = 0.96794; RxnConstantTable(2,1) = 0.89131; RxnConstantTable(2,2) = 0.7291;  RxnConstantTable(2,3) = -3.9555; RxnConstantTable(2,4) = 0.006488;
      RxnConstantTable(3,0) = 0.97646; RxnConstantTable(3,1) = 0.89043; RxnConstantTable(3,2) = 0.74572; RxnConstantTable(3,3) = -3.9642; RxnConstantTable(3,4) = 0.007123;
      RxnConstantTable(4,0) = 0.96188; RxnConstantTable(4,1) = 0.89617; RxnConstantTable(4,2) = 0.72479; RxnConstantTable(4,3) = -3.955;  RxnConstantTable(4,4) = 0.006509;
      RxnConstantTable(5,0) = 0.96921; RxnConstantTable(5,1) = 0.89329; RxnConstantTable(5,2) = 0.73531; RxnConstantTable(5,3) = -3.9596; RxnConstantTable(5,4) = 0.006818;

    } else if (iReaction == 19) {

      //NO + O -> O2 + N
      RxnConstantTable(0,0) = 0.35438;   RxnConstantTable(0,1) = -1.8821; RxnConstantTable(0,2) = -0.72111;  RxnConstantTable(0,3) = -1.1797;   RxnConstantTable(0,4) = -0.030831;
      RxnConstantTable(1,0) = 0.093613;  RxnConstantTable(1,1) = -1.7806; RxnConstantTable(1,2) = -1.0975;   RxnConstantTable(1,3) = -1.0128;   RxnConstantTable(1,4) = -0.041949;
      RxnConstantTable(2,0) = -0.003732; RxnConstantTable(2,1) = -1.7434; RxnConstantTable(2,2) = -1.2394;   RxnConstantTable(2,3) = -0.94952;  RxnConstantTable(2,4) = -0.046182;
      RxnConstantTable(3,0) = 0.004815;  RxnConstantTable(3,1) = -1.7443; RxnConstantTable(3,2) = -1.2227;   RxnConstantTable(3,3) = -0.95824;  RxnConstantTable(3,4) = -0.045545;
      RxnConstantTable(4,0) = -0.009758; RxnConstantTable(4,1) = -1.7386; RxnConstantTable(4,2) = -1.2436;   RxnConstantTable(4,3) = -0.949;    RxnConstantTable(4,4) = -0.046159;
      RxnConstantTable(5,0) = -0.002428; RxnConstantTable(5,1) = -1.7415; RxnConstantTable(5,2) = -1.2331;   RxnConstantTable(5,3) = -0.95365;  RxnConstantTable(5,4) = -0.04585;

    } else if (iReaction == 20) {

      //N + O -> NO+ + e-
      RxnConstantTable(0,0) = -2.1852;   RxnConstantTable(0,1) = -6.6709; RxnConstantTable(0,2) = -4.2968; RxnConstantTable(0,3) = -2.2175; RxnConstantTable(0,4) = -0.050748;
      RxnConstantTable(1,0) = -1.0276;   RxnConstantTable(1,1) = -7.1278; RxnConstantTable(1,2) = -2.637;  RxnConstantTable(1,3) = -2.95;   RxnConstantTable(1,4) = -0.0021;
      RxnConstantTable(2,0) = -0.65871;  RxnConstantTable(2,1) = -7.2742; RxnConstantTable(2,2) = -2.1096; RxnConstantTable(2,3) = -3.1823; RxnConstantTable(2,4) = 0.01331;
      RxnConstantTable(3,0) = -0.57924;  RxnConstantTable(3,1) = -7.3079; RxnConstantTable(3,2) = -1.9999; RxnConstantTable(3,3) = -3.2294; RxnConstantTable(3,4) = 0.016382;
      RxnConstantTable(4,0) = -0.53538;  RxnConstantTable(4,1) = -7.3252; RxnConstantTable(4,2) = -1.937;  RxnConstantTable(4,3) = -3.2572; RxnConstantTable(4,4) = 0.01823;
      RxnConstantTable(5,0) = -0.52801;  RxnConstantTable(5,1) = -7.3281; RxnConstantTable(5,2) = -1.9264; RxnConstantTable(5,3) = -3.2618; RxnConstantTable(5,4) = 0.01854;

    } else if (iReaction == 21) {

      //N2 + e -> N + N + e
      RxnConstantTable(0,0) = 3.4907;  RxnConstantTable(0,1) = 0.83133; RxnConstantTable(0,2) = 4.0978; RxnConstantTable(0,3) = -12.728; RxnConstantTable(0,4) = 0.07487;
      RxnConstantTable(1,0) = 2.0723;  RxnConstantTable(1,1) = 1.3897;  RxnConstantTable(1,2) = 2.0617; RxnConstantTable(1,3) = -11.828; RxnConstantTable(1,4) = 0.015105;
      RxnConstantTable(2,0) = 1.6060;  RxnConstantTable(2,1) = 1.5732;  RxnConstantTable(2,2) = 1.3923; RxnConstantTable(2,3) = -11.533; RxnConstantTable(2,4) = -0.004543;
      RxnConstantTable(3,0) = 1.5351;  RxnConstantTable(3,1) = 1.6061;  RxnConstantTable(3,2) = 1.2993; RxnConstantTable(3,3) = -11.494; RxnConstantTable(3,4) = -0.00698;
      RxnConstantTable(4,0) = 1.4766;  RxnConstantTable(4,1) = 1.6291;  RxnConstantTable(4,2) = 1.2153; RxnConstantTable(4,3) = -11.457; RxnConstantTable(4,4) = -0.009444;
      RxnConstantTable(5,0) = 1.4766;  RxnConstantTable(5,1) = 1.6291;  RxnConstantTable(5,2) = 1.2153; RxnConstantTable(5,3) = -11.457; RxnConstantTable(5,4) = -0.009444;
    }

  } else if (gas_model == "AIR-11"){

    if (iReaction <= 9) {

      //N2 + M -> 2N + M
      RxnConstantTable(0,0) = 3.4907;  RxnConstantTable(0,1) = 0.83133; RxnConstantTable(0,2) = 4.0978;  RxnConstantTable(0,3) = -12.728; RxnConstantTable(0,4) = 0.07487;   //n = 1E14
      RxnConstantTable(1,0) = 2.0723;  RxnConstantTable(1,1) = 1.38970; RxnConstantTable(1,2) = 2.0617;  RxnConstantTable(1,3) = -11.828; RxnConstantTable(1,4) = 0.015105;  //n = 1E15
      RxnConstantTable(2,0) = 1.6060;  RxnConstantTable(2,1) = 1.57320; RxnConstantTable(2,2) = 1.3923;  RxnConstantTable(2,3) = -11.533; RxnConstantTable(2,4) = -0.004543; //n = 1E16
      RxnConstantTable(3,0) = 1.5351;  RxnConstantTable(3,1) = 1.60610; RxnConstantTable(3,2) = 1.2993;  RxnConstantTable(3,3) = -11.494; RxnConstantTable(3,4) = -0.00698;  //n = 1E17
      RxnConstantTable(4,0) = 1.4766;  RxnConstantTable(4,1) = 1.62910; RxnConstantTable(4,2) = 1.2153;  RxnConstantTable(4,3) = -11.457; RxnConstantTable(4,4) = -0.00944;  //n = 1E18
      RxnConstantTable(5,0) = 1.4766;  RxnConstantTable(5,1) = 1.62910; RxnConstantTable(5,2) = 1.2153;  RxnConstantTable(5,3) = -11.457; RxnConstantTable(5,4) = -0.00944;  //n = 1E19

    } else if (iReaction > 9 && iReaction <= 19) {

      //O2 + M -> 2O + M
      RxnConstantTable(0,0) = 1.8103;  RxnConstantTable(0,1) = 1.9607;  RxnConstantTable(0,2) = 3.5716;  RxnConstantTable(0,3) = -7.3623;   RxnConstantTable(0,4) = 0.083861;
      RxnConstantTable(1,0) = 0.91354; RxnConstantTable(1,1) = 2.3160;  RxnConstantTable(1,2) = 2.2885;  RxnConstantTable(1,3) = -6.7969;   RxnConstantTable(1,4) = 0.046338;
      RxnConstantTable(2,0) = 0.64183; RxnConstantTable(2,1) = 2.4253;  RxnConstantTable(2,2) = 1.9026;  RxnConstantTable(2,3) = -6.6277;   RxnConstantTable(2,4) = 0.035151;
      RxnConstantTable(3,0) = 0.55388; RxnConstantTable(3,1) = 2.4600;  RxnConstantTable(3,2) = 1.7763;  RxnConstantTable(3,3) = -6.5720;   RxnConstantTable(3,4) = 0.031445;
      RxnConstantTable(4,0) = 0.52455; RxnConstantTable(4,1) = 2.4715;  RxnConstantTable(4,2) = 1.7342;  RxnConstantTable(4,3) = -6.55534;  RxnConstantTable(4,4) = 0.030209;
      RxnConstantTable(5,0) = 0.50989; RxnConstantTable(5,1) = 2.4773;  RxnConstantTable(5,2) = 1.7132;  RxnConstantTable(5,3) = -6.5441;   RxnConstantTable(5,4) = 0.029591;

    } else if (iReaction > 19 && iReaction <= 29) {

      //NO + M -> N + O + M
      RxnConstantTable(0,0) = 2.1649;  RxnConstantTable(0,1) = 0.078577;  RxnConstantTable(0,2) = 2.8508;  RxnConstantTable(0,3) = -8.5422; RxnConstantTable(0,4) = 0.053043;
      RxnConstantTable(1,0) = 1.0072;  RxnConstantTable(1,1) = 0.53545;   RxnConstantTable(1,2) = 1.1911;  RxnConstantTable(1,3) = -7.8098; RxnConstantTable(1,4) = 0.004394;
      RxnConstantTable(2,0) = 0.63817; RxnConstantTable(2,1) = 0.68189;   RxnConstantTable(2,2) = 0.66336; RxnConstantTable(2,3) = -7.5773; RxnConstantTable(2,4) = -0.011025;
      RxnConstantTable(3,0) = 0.55889; RxnConstantTable(3,1) = 0.71558;   RxnConstantTable(3,2) = 0.55396; RxnConstantTable(3,3) = -7.5304; RxnConstantTable(3,4) = -0.014089;
      RxnConstantTable(4,0) = 0.5150;  RxnConstantTable(4,1) = 0.73286;   RxnConstantTable(4,2) = 0.49096; RxnConstantTable(4,3) = -7.5025; RxnConstantTable(4,4) = -0.015938;
      RxnConstantTable(5,0) = 0.50765; RxnConstantTable(5,1) = 0.73575;   RxnConstantTable(5,2) = 0.48042; RxnConstantTable(5,3) = -7.4979; RxnConstantTable(5,4) = -0.016247;

    } else if (iReaction == 30) {

      //N2 + O -> NO + N
      RxnConstantTable(0,0) = 1.3261;  RxnConstantTable(0,1) = 0.75268; RxnConstantTable(0,2) = 1.2474;  RxnConstantTable(0,3) = -4.1857; RxnConstantTable(0,4) = 0.02184;
      RxnConstantTable(1,0) = 1.0653;  RxnConstantTable(1,1) = 0.85417; RxnConstantTable(1,2) = 0.87093; RxnConstantTable(1,3) = -4.0188; RxnConstantTable(1,4) = 0.010721;
      RxnConstantTable(2,0) = 0.96794; RxnConstantTable(2,1) = 0.89131; RxnConstantTable(2,2) = 0.7291;  RxnConstantTable(2,3) = -3.9555; RxnConstantTable(2,4) = 0.006488;
      RxnConstantTable(3,0) = 0.97646; RxnConstantTable(3,1) = 0.89043; RxnConstantTable(3,2) = 0.74572; RxnConstantTable(3,3) = -3.9642; RxnConstantTable(3,4) = 0.007123;
      RxnConstantTable(4,0) = 0.96188; RxnConstantTable(4,1) = 0.89617; RxnConstantTable(4,2) = 0.72479; RxnConstantTable(4,3) = -3.955;  RxnConstantTable(4,4) = 0.006509;
      RxnConstantTable(5,0) = 0.96921; RxnConstantTable(5,1) = 0.89329; RxnConstantTable(5,2) = 0.73531; RxnConstantTable(5,3) = -3.9596; RxnConstantTable(5,4) = 0.006818;

    } else if (iReaction == 31) {

      //NO + O -> O2 + N
      RxnConstantTable(0,0) = 0.35438;   RxnConstantTable(0,1) = -1.8821; RxnConstantTable(0,2) = -0.72111;  RxnConstantTable(0,3) = -1.1797;   RxnConstantTable(0,4) = -0.030831;
      RxnConstantTable(1,0) = 0.093613;  RxnConstantTable(1,1) = -1.7806; RxnConstantTable(1,2) = -1.0975;   RxnConstantTable(1,3) = -1.0128;   RxnConstantTable(1,4) = -0.041949;
      RxnConstantTable(2,0) = -0.003732; RxnConstantTable(2,1) = -1.7434; RxnConstantTable(2,2) = -1.2394;   RxnConstantTable(2,3) = -0.94952;  RxnConstantTable(2,4) = -0.046182;
      RxnConstantTable(3,0) = 0.004815;  RxnConstantTable(3,1) = -1.7443; RxnConstantTable(3,2) = -1.2227;   RxnConstantTable(3,3) = -0.95824;  RxnConstantTable(3,4) = -0.045545;
      RxnConstantTable(4,0) = -0.009758; RxnConstantTable(4,1) = -1.7386; RxnConstantTable(4,2) = -1.2436;   RxnConstantTable(4,3) = -0.949;    RxnConstantTable(4,4) = -0.046159;
      RxnConstantTable(5,0) = -0.002428; RxnConstantTable(5,1) = -1.7415; RxnConstantTable(5,2) = -1.2331;   RxnConstantTable(5,3) = -0.95365;  RxnConstantTable(5,4) = -0.04585;

    } else if (iReaction == 32) {

      //N + O -> NO+ + e-
      RxnConstantTable(0,0) = -2.1852;   RxnConstantTable(0,1) = -6.6709; RxnConstantTable(0,2) = -4.2968; RxnConstantTable(0,3) = -2.2175; RxnConstantTable(0,4) = -0.050748;
      RxnConstantTable(1,0) = -1.0276;   RxnConstantTable(1,1) = -7.1278; RxnConstantTable(1,2) = -2.637;  RxnConstantTable(1,3) = -2.95;   RxnConstantTable(1,4) = -0.0021;
      RxnConstantTable(2,0) = -0.65871;  RxnConstantTable(2,1) = -7.2742; RxnConstantTable(2,2) = -2.1096; RxnConstantTable(2,3) = -3.1823; RxnConstantTable(2,4) = 0.01331;
      RxnConstantTable(3,0) = -0.57924;  RxnConstantTable(3,1) = -7.3079; RxnConstantTable(3,2) = -1.9999; RxnConstantTable(3,3) = -3.2294; RxnConstantTable(3,4) = 0.016382;
      RxnConstantTable(4,0) = -0.53538;  RxnConstantTable(4,1) = -7.3252; RxnConstantTable(4,2) = -1.937;  RxnConstantTable(4,3) = -3.2572; RxnConstantTable(4,4) = 0.01823;
      RxnConstantTable(5,0) = -0.52801;  RxnConstantTable(5,1) = -7.3281; RxnConstantTable(5,2) = -1.9264; RxnConstantTable(5,3) = -3.2618; RxnConstantTable(5,4) = 0.01854;

    } else if (iReaction == 33) {

      //N2 + e -> N + N + e
      RxnConstantTable(0,0) = 3.4907;  RxnConstantTable(0,1) = 0.83133; RxnConstantTable(0,2) = 4.0978; RxnConstantTable(0,3) = -12.728; RxnConstantTable(0,4) = 0.07487;
      RxnConstantTable(1,0) = 2.0723;  RxnConstantTable(1,1) = 1.3897;  RxnConstantTable(1,2) = 2.0617; RxnConstantTable(1,3) = -11.828; RxnConstantTable(1,4) = 0.015105;
      RxnConstantTable(2,0) = 1.6060;  RxnConstantTable(2,1) = 1.5732;  RxnConstantTable(2,2) = 1.3923; RxnConstantTable(2,3) = -11.533; RxnConstantTable(2,4) = -0.004543;
      RxnConstantTable(3,0) = 1.5351;  RxnConstantTable(3,1) = 1.6061;  RxnConstantTable(3,2) = 1.2993; RxnConstantTable(3,3) = -11.494; RxnConstantTable(3,4) = -0.00698;
      RxnConstantTable(4,0) = 1.4766;  RxnConstantTable(4,1) = 1.6291;  RxnConstantTable(4,2) = 1.2153; RxnConstantTable(4,3) = -11.457; RxnConstantTable(4,4) = -0.009444;
      RxnConstantTable(5,0) = 1.4766;  RxnConstantTable(5,1) = 1.6291;  RxnConstantTable(5,2) = 1.2153; RxnConstantTable(5,3) = -11.457; RxnConstantTable(5,4) = -0.009444;

  } else if (iReaction == 34) {

      // O + O -> O2+ + e-
      RxnConstantTable(0,0) = -0.11682;  RxnConstantTable(0,1) =  -7.6883;  RxnConstantTable(0,2) =  -2.2498;  RxnConstantTable(0,3) =  -7.7905;  RxnConstantTable(0,4) = -0.011079;  //n = 1E14
      RxnConstantTable(1,0) =  0.77986;  RxnConstantTable(1,1) =  -8.0436;  RxnConstantTable(1,2) = -0.96678;  RxnConstantTable(1,3) =  -8.3559;  RxnConstantTable(1,4) =  0.026440;  //n = 1E15
      RxnConstantTable(2,0) =   1.0516;  RxnConstantTable(2,1) =  -8.1530;  RxnConstantTable(2,2) = -0.58082;  RxnConstantTable(2,3) =  -8.5251;  RxnConstantTable(2,4) =  0.037629;  //n = 1E16
      RxnConstantTable(3,0) =   1.1395;  RxnConstantTable(3,1) =  -8.1876;  RxnConstantTable(3,2) = -0.45461;  RxnConstantTable(3,3) =  -8.5808;  RxnConstantTable(3,4) =  0.041333;  //n = 1E17
      RxnConstantTable(4,0) =   1.1689;  RxnConstantTable(4,1) =  -8.1991;  RxnConstantTable(4,2) = -0.41245;  RxnConstantTable(4,3) =  -8.5995;  RxnConstantTable(4,4) =  0.042571;  //n = 1E18
      RxnConstantTable(5,0) =   1.1835;  RxnConstantTable(5,1) =  -8.2049;  RxnConstantTable(5,2) = -0.39146;  RxnConstantTable(5,3) =  -8.6087;  RxnConstantTable(5,4) =  0.043187;  //n = 1E19

    } else if (iReaction == 35) {

      // N + N -> N2+ + e-
      RxnConstantTable(0,0) =  -4.3785;  RxnConstantTable(0,1) =  -4.2726;  RxnConstantTable(0,2) =  -7.8709;  RxnConstantTable(0,3) =  -4.4628;  RxnConstantTable(0,4) = -0.124020;  //n = 1E14
      RxnConstantTable(1,0) =  -2.9601;  RxnConstantTable(1,1) =  -4.8310;  RxnConstantTable(1,2) =  -5.8348;  RxnConstantTable(1,3) =  -5.3621;  RxnConstantTable(1,4) = -0.064252;  //n = 1E15
      RxnConstantTable(2,0) =  -2.4938;  RxnConstantTable(2,1) =  -5.0145;  RxnConstantTable(2,2) =  -5.1654;  RxnConstantTable(2,3) =  -5.6577;  RxnConstantTable(2,4) = -0.044602;  //n = 1E16
      RxnConstantTable(3,0) =  -2.4229;  RxnConstantTable(3,1) =  -5.0474;  RxnConstantTable(3,2) =  -5.0724;  RxnConstantTable(3,3) =  -5.6961;  RxnConstantTable(3,4) = -0.042167;  //n = 1E17
      RxnConstantTable(4,0) =  -2.3644;  RxnConstantTable(4,1) =  -5.0704;  RxnConstantTable(4,2) =  -4.9885;  RxnConstantTable(4,3) =  -5.7332;  RxnConstantTable(4,4) = -0.039703;  //n = 1E18
      RxnConstantTable(5,0) =  -2.3644;  RxnConstantTable(5,1) =  -5.0704;  RxnConstantTable(5,2) =  -4.9885;  RxnConstantTable(5,3) =  -5.7332;  RxnConstantTable(5,4) = -0.039703;  //n = 1E19

    } else if (iReaction == 36) {

      // NO+ + O -> N+ + O2
      RxnConstantTable(0,0) =  -1.5349;  RxnConstantTable(0,1) =   1.6836;  RxnConstantTable(0,2) =  -2.9690;  RxnConstantTable(0,3) =  -6.4640;  RxnConstantTable(0,4) = -0.083316;  //n = 1E14
      RxnConstantTable(1,0) =  -1.0864;  RxnConstantTable(1,1) =   1.5059;  RxnConstantTable(1,2) =  -2.3273;  RxnConstantTable(1,3) =  -6.7468;  RxnConstantTable(1,4) = -0.064551;  //n = 1E15
      RxnConstantTable(2,0) = -0.95072;  RxnConstantTable(2,1) =   1.4513;  RxnConstantTable(2,2) =  -2.1346;  RxnConstantTable(2,3) =  -6.8313;  RxnConstantTable(2,4) = -0.058964;  //n = 1E16
      RxnConstantTable(3,0) = -0.90672;  RxnConstantTable(3,1) =   1.4340;  RxnConstantTable(3,2) =  -2.0714;  RxnConstantTable(3,3) =  -6.8592;  RxnConstantTable(3,4) = -0.057110;  //n = 1E17
      RxnConstantTable(4,0) = -0.89206;  RxnConstantTable(4,1) =   1.4282;  RxnConstantTable(4,2) =  -2.0504;  RxnConstantTable(4,3) =  -6.8685;  RxnConstantTable(4,4) = -0.056493;  //n = 1E18
      RxnConstantTable(5,0) = -0.88472;  RxnConstantTable(5,1) =   1.4254;  RxnConstantTable(5,2) =  -2.0398;  RxnConstantTable(5,3) =  -6.8731;  RxnConstantTable(5,4) = -0.056184;  //n = 1E19

    //} else if (iReaction == 37) {

      // N+ + N2 -> N2+ + N  (not in Tables A13-A15, forward only)  Using Mutation++ air-11 data
      //RxnConstantTable(0,0) = -0.00271;  RxnConstantTable(0,1) = 0.02134;  RxnConstantTable(0,2) = -1.3825;  RxnConstantTable(0,3) = -0.64735;  RxnConstantTable(0,4) = -0.02472;
      //RxnConstantTable(1,0) = -0.00271;  RxnConstantTable(1,1) = 0.02134;  RxnConstantTable(1,2) = -1.3825;  RxnConstantTable(1,3) = -0.64735;  RxnConstantTable(1,4) = -0.02472;
      //RxnConstantTable(2,0) = -0.00271;  RxnConstantTable(2,1) = 0.02134;  RxnConstantTable(2,2) = -1.3825;  RxnConstantTable(2,3) = -0.64735;  RxnConstantTable(2,4) = -0.02472;
      //RxnConstantTable(3,0) = -0.00271;  RxnConstantTable(3,1) = 0.02134;  RxnConstantTable(3,2) = -1.3825;  RxnConstantTable(3,3) = -0.64735;  RxnConstantTable(3,4) = -0.02472;
      //RxnConstantTable(4,0) = -0.00271;  RxnConstantTable(4,1) = 0.02134;  RxnConstantTable(4,2) = -1.3825;  RxnConstantTable(4,3) = -0.64735;  RxnConstantTable(4,4) = -0.02472;
      //RxnConstantTable(5,0) = -0.00271;  RxnConstantTable(5,1) = 0.02134;  RxnConstantTable(5,2) = -1.3825;  RxnConstantTable(5,3) = -0.64735;  RxnConstantTable(5,4) = -0.02472;

      // N+ + N2 -> N2+ + N   (refit from NASA RP-1232, valid 2000-30000 K)
      //RxnConstantTable(0,0) = -2.85898;  RxnConstantTable(0,1) = -2.69958;  RxnConstantTable(0,2) = -13.5303; RxnConstantTable(0,3) = 6.05649;   RxnConstantTable(0,4) = -0.483583; //n = 1E14
      //RxnConstantTable(1,0) = -1.92614;  RxnConstantTable(1,1) = -2.67693;  RxnConstantTable(1,2) = -10.5971; RxnConstantTable(1,3) = 4.69519;   RxnConstantTable(1,4) = -0.402935; //n = 1E15
      //RxnConstantTable(2,0) = -1.01745;  RxnConstantTable(2,1) = -2.36419;  RxnConstantTable(2,2) = -7.47922; RxnConstantTable(2,3) = 3.14358;   RxnConstantTable(2,4) = -0.308171; //n = 1E16
      //RxnConstantTable(3,0) = -0.323931; RxnConstantTable(3,1) = -1.86662;  RxnConstantTable(3,2) = -4.81154; RxnConstantTable(3,3) = 1.73391;   RxnConstantTable(3,4) = -0.219949; //n = 1E17
      //RxnConstantTable(4,0) = 0.146971;  RxnConstantTable(4,1) = -1.18386;  RxnConstantTable(4,2) = -2.52465; RxnConstantTable(4,3) = 0.437343;  RxnConstantTable(4,4) = -0.136602; //n = 1E18
      //RxnConstantTable(5,0) = 0.230879;  RxnConstantTable(5,1) = -0.644595; RxnConstantTable(5,2) = -1.47619; RxnConstantTable(5,3) = -0.231832; RxnConstantTable(5,4) = -0.091808; //n = 1E19

    } else if (iReaction == 37) {

      // O2+ + N -> N+ + O2
      RxnConstantTable(0,0) =  -3.6030;  RxnConstantTable(0,1) =   2.7010;  RxnConstantTable(0,2) =  -5.0155;  RxnConstantTable(0,3) = -0.89125;  RxnConstantTable(0,4) = -0.122970;  //n = 1E14
      RxnConstantTable(1,0) =  -2.8938;  RxnConstantTable(1,1) =   2.4218;  RxnConstantTable(1,2) =  -3.9975;  RxnConstantTable(1,3) =  -1.3409;  RxnConstantTable(1,4) = -0.093088;  //n = 1E15
      RxnConstantTable(2,0) =  -2.6607;  RxnConstantTable(2,1) =   2.3300;  RxnConstantTable(2,2) =  -3.6628;  RxnConstantTable(2,3) =  -1.4887;  RxnConstantTable(2,4) = -0.083264;  //n = 1E16
      RxnConstantTable(3,0) =  -2.6252;  RxnConstantTable(3,1) =   2.3136;  RxnConstantTable(3,2) =  -3.6163;  RxnConstantTable(3,3) =  -1.5079;  RxnConstantTable(3,4) = -0.082048;  //n = 1E17
      RxnConstantTable(4,0) =  -2.5960;  RxnConstantTable(4,1) =   2.3021;  RxnConstantTable(4,2) =  -3.5744;  RxnConstantTable(4,3) =  -1.5264;  RxnConstantTable(4,4) = -0.080816;  //n = 1E18
      RxnConstantTable(5,0) =  -2.5960;  RxnConstantTable(5,1) =   2.3021;  RxnConstantTable(5,2) =  -3.5744;  RxnConstantTable(5,3) =  -1.5264;  RxnConstantTable(5,4) = -0.080816;  //n = 1E19

    } else if (iReaction == 38) {

      // O+ + NO -> N+ + O2
      RxnConstantTable(0,0) =  -1.6355;  RxnConstantTable(0,1) =  0.83058;  RxnConstantTable(0,2) =  -2.9952;  RxnConstantTable(0,3) =  -1.3794;  RxnConstantTable(0,4) = -0.079927;  //n = 1E14
      RxnConstantTable(1,0) =  -1.6355;  RxnConstantTable(1,1) =  0.83058;  RxnConstantTable(1,2) =  -2.9952;  RxnConstantTable(1,3) =  -1.3794;  RxnConstantTable(1,4) = -0.079927;  //n = 1E15
      RxnConstantTable(2,0) =  -1.6355;  RxnConstantTable(2,1) =  0.83058;  RxnConstantTable(2,2) =  -2.9952;  RxnConstantTable(2,3) =  -1.3794;  RxnConstantTable(2,4) = -0.079927;  //n = 1E16
      RxnConstantTable(3,0) =  -1.6355;  RxnConstantTable(3,1) =  0.83058;  RxnConstantTable(3,2) =  -2.9952;  RxnConstantTable(3,3) =  -1.3794;  RxnConstantTable(3,4) = -0.079927;  //n = 1E17
      RxnConstantTable(4,0) =  -1.6355;  RxnConstantTable(4,1) =  0.83058;  RxnConstantTable(4,2) =  -2.9952;  RxnConstantTable(4,3) =  -1.3794;  RxnConstantTable(4,4) = -0.079927;  //n = 1E18
      RxnConstantTable(5,0) =  -1.6355;  RxnConstantTable(5,1) =  0.83058;  RxnConstantTable(5,2) =  -2.9952;  RxnConstantTable(5,3) =  -1.3794;  RxnConstantTable(5,4) = -0.079927;  //n = 1E19

    } else if (iReaction == 39) {

      // O2+ + N2 -> N2+ + O2
      RxnConstantTable(0,0) =  -2.5811;  RxnConstantTable(0,1) =   2.2863;  RxnConstantTable(0,2) =  -5.0946;  RxnConstantTable(0,3) =  -2.0378;  RxnConstantTable(0,4) = -0.121920;  //n = 1E14
      RxnConstantTable(1,0) =  -2.5811;  RxnConstantTable(1,1) =   2.2863;  RxnConstantTable(1,2) =  -5.0946;  RxnConstantTable(1,3) =  -2.0378;  RxnConstantTable(1,4) = -0.121920;  //n = 1E15
      RxnConstantTable(2,0) =  -2.5811;  RxnConstantTable(2,1) =   2.2863;  RxnConstantTable(2,2) =  -5.0946;  RxnConstantTable(2,3) =  -2.0378;  RxnConstantTable(2,4) = -0.121920;  //n = 1E16
      RxnConstantTable(3,0) =  -2.5811;  RxnConstantTable(3,1) =   2.2863;  RxnConstantTable(3,2) =  -5.0946;  RxnConstantTable(3,3) =  -2.0378;  RxnConstantTable(3,4) = -0.121920;  //n = 1E17
      RxnConstantTable(4,0) =  -2.5811;  RxnConstantTable(4,1) =   2.2863;  RxnConstantTable(4,2) =  -5.0946;  RxnConstantTable(4,3) =  -2.0378;  RxnConstantTable(4,4) = -0.121920;  //n = 1E18
      RxnConstantTable(5,0) =  -2.5811;  RxnConstantTable(5,1) =   2.2863;  RxnConstantTable(5,2) =  -5.0946;  RxnConstantTable(5,3) =  -2.0378;  RxnConstantTable(5,4) = -0.121920;  //n = 1E19

    //} else if (iReaction == 41) {

      // O2+ + O -> O+ + O2  (not in Tables A13-A15, forward only)  Using Mutation++ air-11 data
      //RxnConstantTable(0,0) = -0.6531;  RxnConstantTable(0,1) = -0.36419;  RxnConstantTable(0,2) = -1.56307;  RxnConstantTable(0,3) = -1.197;  RxnConstantTable(0,4) = -0.02707;
      //RxnConstantTable(1,0) = -0.6531;  RxnConstantTable(1,1) = -0.36419;  RxnConstantTable(1,2) = -1.56307;  RxnConstantTable(1,3) = -1.197;  RxnConstantTable(1,4) = -0.02707;
      //RxnConstantTable(2,0) = -0.6531;  RxnConstantTable(2,1) = -0.36419;  RxnConstantTable(2,2) = -1.56307;  RxnConstantTable(2,3) = -1.197;  RxnConstantTable(2,4) = -0.02707;
      //RxnConstantTable(3,0) = -0.6531;  RxnConstantTable(3,1) = -0.36419;  RxnConstantTable(3,2) = -1.56307;  RxnConstantTable(3,3) = -1.197;  RxnConstantTable(3,4) = -0.02707;
      //RxnConstantTable(4,0) = -0.6531;  RxnConstantTable(4,1) = -0.36419;  RxnConstantTable(4,2) = -1.56307;  RxnConstantTable(4,3) = -1.197;  RxnConstantTable(4,4) = -0.02707;
      //RxnConstantTable(5,0) = -0.6531;  RxnConstantTable(5,1) = -0.36419;  RxnConstantTable(5,2) = -1.56307;  RxnConstantTable(5,3) = -1.197;  RxnConstantTable(5,4) = -0.02707;

      //O + O2+ -> O2 + O+   (refit from NASA RP-1232, valid 2000-30000 K)
      //RxnConstantTable(0,0) = 0.397867;  RxnConstantTable(0,1) = 6.25044;   RxnConstantTable(0,2) = 9.50945;  RxnConstantTable(0,3) = -10.8487;  RxnConstantTable(0,4) = 0.712607;  //n = 1E14
      //RxnConstantTable(1,0) = -0.456737; RxnConstantTable(1,1) = 5.83809;   RxnConstantTable(1,2) = 6.45651;  RxnConstantTable(1,3) = -9.21609;  RxnConstantTable(1,4) = 0.601659;  //n = 1E15
      //RxnConstantTable(2,0) = -1.08176;  RxnConstantTable(2,1) = 4.75297;   RxnConstantTable(2,2) = 3.44499;  RxnConstantTable(2,3) = -7.22751;  RxnConstantTable(2,4) = 0.444605;  //n = 1E16
      //RxnConstantTable(3,0) = -1.38293;  RxnConstantTable(3,1) = 3.3259;    RxnConstantTable(3,2) = 1.04416;  RxnConstantTable(3,3) = -5.33071;  RxnConstantTable(3,4) = 0.280138;  //n = 1E17
      //RxnConstantTable(4,0) = -1.38452;  RxnConstantTable(4,1) = 1.63499;   RxnConstantTable(4,2) = -0.82487; RxnConstantTable(4,3) = -3.5436;   RxnConstantTable(4,4) = 0.113324;  //n = 1E18
      //RxnConstantTable(5,0) = -1.1651;   RxnConstantTable(5,1) = 0.532367;  RxnConstantTable(5,2) = -1.5107;  RxnConstantTable(5,3) = -2.64254;  RxnConstantTable(5,4) = 0.0219787; //n = 1E19

    } else if (iReaction == 40) {

      // NO+ + N -> O+ + N2
      RxnConstantTable(0,0) =  -1.2255;  RxnConstantTable(0,1) =  0.10039;  RxnConstantTable(0,2) =  -1.2212;  RxnConstantTable(0,3) = -0.89883;  RxnConstantTable(0,4) = -0.025232;  //n = 1E14
      RxnConstantTable(1,0) = -0.51629;  RxnConstantTable(1,1) = -0.17877;  RxnConstantTable(1,2) = -0.20321;  RxnConstantTable(1,3) =  -1.3485;  RxnConstantTable(1,4) =  0.004649;  //n = 1E15
      RxnConstantTable(2,0) = -0.28311;  RxnConstantTable(2,1) = -0.27056;  RxnConstantTable(2,2) =  0.13152;  RxnConstantTable(2,3) =  -1.4963;  RxnConstantTable(2,4) =  0.014474;  //n = 1E16
      RxnConstantTable(3,0) = -0.24765;  RxnConstantTable(3,1) = -0.28699;  RxnConstantTable(3,2) =  0.17802;  RxnConstantTable(3,3) =  -1.5155;  RxnConstantTable(3,4) =  0.015692;  //n = 1E17
      RxnConstantTable(4,0) = -0.21842;  RxnConstantTable(4,1) = -0.29849;  RxnConstantTable(4,2) =  0.21998;  RxnConstantTable(4,3) =  -1.5340;  RxnConstantTable(4,4) =  0.016923;  //n = 1E18
      RxnConstantTable(5,0) = -0.21842;  RxnConstantTable(5,1) = -0.29849;  RxnConstantTable(5,2) =  0.21998;  RxnConstantTable(5,3) =  -1.5340;  RxnConstantTable(5,4) =  0.016923;  //n = 1E19

    } else if (iReaction == 41) {

      // NO+ + O2 -> O2+ + NO
      RxnConstantTable(0,0) =   1.7139;  RxnConstantTable(0,1) =  0.86469;  RxnConstantTable(0,2) =   2.7679;  RxnConstantTable(0,3) =  -4.3932;  RxnConstantTable(0,4) =  0.070493;  //n = 1E14
      RxnConstantTable(1,0) =   1.7139;  RxnConstantTable(1,1) =  0.86469;  RxnConstantTable(1,2) =   2.7679;  RxnConstantTable(1,3) =  -4.3932;  RxnConstantTable(1,4) =  0.070493;  //n = 1E15
      RxnConstantTable(2,0) =   1.7139;  RxnConstantTable(2,1) =  0.86469;  RxnConstantTable(2,2) =   2.7679;  RxnConstantTable(2,3) =  -4.3932;  RxnConstantTable(2,4) =  0.070493;  //n = 1E16
      RxnConstantTable(3,0) =   1.7139;  RxnConstantTable(3,1) =  0.86469;  RxnConstantTable(3,2) =   2.7679;  RxnConstantTable(3,3) =  -4.3932;  RxnConstantTable(3,4) =  0.070493;  //n = 1E17
      RxnConstantTable(4,0) =   1.7139;  RxnConstantTable(4,1) =  0.86469;  RxnConstantTable(4,2) =   2.7679;  RxnConstantTable(4,3) =  -4.3932;  RxnConstantTable(4,4) =  0.070493;  //n = 1E18
      RxnConstantTable(5,0) =   1.7139;  RxnConstantTable(5,1) =  0.86469;  RxnConstantTable(5,2) =   2.7679;  RxnConstantTable(5,3) =  -4.3932;  RxnConstantTable(5,4) =  0.070493;  //n = 1E19

    } else if (iReaction == 42) {

      // NO+ + O -> O2+ + N
      RxnConstantTable(0,0) =   2.0681;  RxnConstantTable(0,1) =  -1.0173;  RxnConstantTable(0,2) =   2.0466;  RxnConstantTable(0,3) =  -5.5728;  RxnConstantTable(0,4) =  0.039655;  //n = 1E14
      RxnConstantTable(1,0) =   1.8073;  RxnConstantTable(1,1) = -0.91584;  RxnConstantTable(1,2) =   1.6701;  RxnConstantTable(1,3) =  -5.4058;  RxnConstantTable(1,4) =  0.028533;  //n = 1E15
      RxnConstantTable(2,0) =   1.7100;  RxnConstantTable(2,1) = -0.87869;  RxnConstantTable(2,2) =   1.5282;  RxnConstantTable(2,3) =  -5.3426;  RxnConstantTable(2,4) =  0.024301;  //n = 1E16
      RxnConstantTable(3,0) =   1.7185;  RxnConstantTable(3,1) = -0.87958;  RxnConstantTable(3,2) =   1.5449;  RxnConstantTable(3,3) =  -5.3513;  RxnConstantTable(3,4) =  0.024936;  //n = 1E17
      RxnConstantTable(4,0) =   1.7039;  RxnConstantTable(4,1) = -0.87383;  RxnConstantTable(4,2) =   1.5239;  RxnConstantTable(4,3) =  -5.3420;  RxnConstantTable(4,4) =  0.024321;  //n = 1E18
      RxnConstantTable(5,0) =   1.7112;  RxnConstantTable(5,1) = -0.87672;  RxnConstantTable(5,2) =   1.5345;  RxnConstantTable(5,3) =  -5.3467;  RxnConstantTable(5,4) =  0.024631;  //n = 1E19

    } else if (iReaction == 43) {

      // O+ + N2 -> N2+ + O
      RxnConstantTable(0,0) = -0.96795;  RxnConstantTable(0,1) =   2.2979;  RxnConstantTable(0,2) =  -2.3531;  RxnConstantTable(0,3) =  -1.3463;  RxnConstantTable(0,4) = -0.048042;  //n = 1E14
      RxnConstantTable(1,0) =  -1.4164;  RxnConstantTable(1,1) =   2.4756;  RxnConstantTable(1,2) =  -2.9947;  RxnConstantTable(1,3) =  -1.0636;  RxnConstantTable(1,4) = -0.066805;  //n = 1E15
      RxnConstantTable(2,0) =  -1.5522;  RxnConstantTable(2,1) =   2.5303;  RxnConstantTable(2,2) =  -3.1876;  RxnConstantTable(2,3) = -0.97903;  RxnConstantTable(2,4) = -0.072396;  //n = 1E16
      RxnConstantTable(3,0) =  -1.5962;  RxnConstantTable(3,1) =   2.5476;  RxnConstantTable(3,2) =  -3.2507;  RxnConstantTable(3,3) = -0.95116;  RxnConstantTable(3,4) = -0.074249;  //n = 1E17
      RxnConstantTable(4,0) =  -1.6108;  RxnConstantTable(4,1) =   2.5533;  RxnConstantTable(4,2) =  -3.2718;  RxnConstantTable(4,3) = -0.94186;  RxnConstantTable(4,4) = -0.074867;  //n = 1E18
      RxnConstantTable(5,0) =  -1.6181;  RxnConstantTable(5,1) =   2.5562;  RxnConstantTable(5,2) =  -3.2823;  RxnConstantTable(5,3) = -0.93721;  RxnConstantTable(5,4) = -0.075176;  //n = 1E19

    } else if (iReaction == 44) {

      // NO+ + N -> N2+ + O
      RxnConstantTable(0,0) =  -2.1934;  RxnConstantTable(0,1) =   2.3983;  RxnConstantTable(0,2) =  -3.5743;  RxnConstantTable(0,3) =  -2.2452;  RxnConstantTable(0,4) = -0.073271;  //n = 1E14
      RxnConstantTable(1,0) =  -1.9325;  RxnConstantTable(1,1) =   2.2968;  RxnConstantTable(1,2) =  -3.1978;  RxnConstantTable(1,3) =  -2.4121;  RxnConstantTable(1,4) = -0.062149;  //n = 1E15
      RxnConstantTable(2,0) =  -1.8352;  RxnConstantTable(2,1) =   2.2597;  RxnConstantTable(2,2) =  -3.0560;  RxnConstantTable(2,3) =  -2.4754;  RxnConstantTable(2,4) = -0.057919;  //n = 1E16
      RxnConstantTable(3,0) =  -1.8438;  RxnConstantTable(3,1) =   2.2606;  RxnConstantTable(3,2) =  -3.0726;  RxnConstantTable(3,3) =  -2.4667;  RxnConstantTable(3,4) = -0.058554;  //n = 1E17
      RxnConstantTable(4,0) =  -1.8292;  RxnConstantTable(4,1) =   2.2548;  RxnConstantTable(4,2) =  -3.0517;  RxnConstantTable(4,3) =  -2.4759;  RxnConstantTable(4,4) = -0.057940;  //n = 1E18
      RxnConstantTable(5,0) =  -1.8365;  RxnConstantTable(5,1) =   2.2577;  RxnConstantTable(5,2) =  -3.0622;  RxnConstantTable(5,3) =  -2.4713;  RxnConstantTable(5,4) = -0.058248;  //n = 1E19

    } else if (iReaction == 45) {

      // O + e- -> O+ + e- + e-
      RxnConstantTable(0,0) =  0.08045;  RxnConstantTable(0,1) =  -5.7393;  RxnConstantTable(0,2) =  -1.4195;  RxnConstantTable(0,3) =  -15.844;  RxnConstantTable(0,4) = -0.001087;  //n = 1E14
      RxnConstantTable(1,0) =  0.52883;  RxnConstantTable(1,1) =  -5.9170;  RxnConstantTable(1,2) = -0.77795;  RxnConstantTable(1,3) =  -16.127;  RxnConstantTable(1,4) =  0.017675;  //n = 1E15
      RxnConstantTable(2,0) =  0.66478;  RxnConstantTable(2,1) =  -5.9716;  RxnConstantTable(2,2) = -0.58486;  RxnConstantTable(2,3) =  -16.212;  RxnConstantTable(2,4) =  0.023273;  //n = 1E16
      RxnConstantTable(3,0) =  0.70879;  RxnConstantTable(3,1) =  -5.9890;  RxnConstantTable(3,2) = -0.52169;  RxnConstantTable(3,3) =  -16.240;  RxnConstantTable(3,4) =  0.025127;  //n = 1E17
      RxnConstantTable(4,0) =  0.72341;  RxnConstantTable(4,1) =  -5.9947;  RxnConstantTable(4,2) =  -0.5007;  RxnConstantTable(4,3) =  -16.249;  RxnConstantTable(4,4) =  0.025743;  //n = 1E18
      RxnConstantTable(5,0) =  0.73078;  RxnConstantTable(5,1) =  -5.9976;  RxnConstantTable(5,2) = -0.49012;  RxnConstantTable(5,3) =  -16.254;  RxnConstantTable(5,4) =  0.026054;  //n = 1E19

    } else if (iReaction == 46) {

      // N + e- -> N+ + e- + e-
      RxnConstantTable(0,0) =  -1.9094;  RxnConstantTable(0,1) =  -3.0267;  RxnConstantTable(0,2) =  -3.6935;  RxnConstantTable(0,3) =  -16.044;  RxnConstantTable(0,4) = -0.050183;  //n = 1E14
      RxnConstantTable(1,0) =  -1.2002;  RxnConstantTable(1,1) =  -3.3059;  RxnConstantTable(1,2) =  -2.6755;  RxnConstantTable(1,3) =  -16.494;  RxnConstantTable(1,4) = -0.020301;  //n = 1E15
      RxnConstantTable(2,0) = -0.96709;  RxnConstantTable(2,1) =  -3.3976;  RxnConstantTable(2,2) =  -2.3408;  RxnConstantTable(2,3) =  -16.642;  RxnConstantTable(2,4) = -0.010477;  //n = 1E16
      RxnConstantTable(3,0) = -0.93184;  RxnConstantTable(3,1) =  -3.4140;  RxnConstantTable(3,2) =  -2.2946;  RxnConstantTable(3,3) =  -16.661;  RxnConstantTable(3,4) = -0.009269;  //n = 1E17
      RxnConstantTable(4,0) =  -0.9026;  RxnConstantTable(4,1) =  -3.4255;  RxnConstantTable(4,2) =  -2.2526;  RxnConstantTable(4,3) =  -16.679;  RxnConstantTable(4,4) = -0.008037;  //n = 1E18
      RxnConstantTable(5,0) =  -0.9026;  RxnConstantTable(5,1) =  -3.4255;  RxnConstantTable(5,2) =  -2.2526;  RxnConstantTable(5,3) =  -16.679;  RxnConstantTable(5,4) = -0.008037;  //n = 1E19

    }
  }

}

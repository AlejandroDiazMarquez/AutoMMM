#!/bin/bash

function main () {
MainDir=$PWD

#*******************************
# ideally put  CADSS.x  binary and shared folder in the CADSS_DIR  

CADSS_DIR=/home/adiaz/programs/CADSS
#CADSS_DIR=${HOME}/soft/CADSS_V3
#*******************************

CADSSBinary=${CADSS_DIR}/scripts/CADSS_V3.1_2021-long.X
[ ! -e ${CADSSBinary} ] &&  {
   echo" ${CADSSBinary} not found"
   exit
   }

#let's keep only one shared dir for all CADSSS calculations, we use softlinks instead of multiple copies 
shared=${CADSS_DIR}/share
[ ! -e ${shared} ] &&  {
   echo" ${shared} not found"
   exit
   } || ln -sf ${shared} .

System="Systema"
ForcefieldDir="${System}"
#ExtraTag="fixed" # When the job name becomes messy, use a simple and distinguishable one 

T=300

Adsorbate1="CO2"
MoleculeDefinition="TraPPE"

x1=1.00 # mole fraction of constituent 1
IDRW_x1=1.00000000         #IdealGasRosenbluthWeights

SimType="gcmc" # widom or gcmc or NVT
NeedMakeGrid="no"
[ "$NeedMakeGrid" == "yes" ] && SimType="widom" 

UseTabularGrid="no" 
WorkingDir="$MainDir/T${T}K-$SimType"
mkdir -p $WorkingDir


[ "$NeedMakeGrid" == "yes" ] || [ "$SimType" == 'widom' ]  && { 
    CalcDir=$WorkingDir/$System-2E6
    InputDir=$CalcDir/input
    mkdir -p $CalcDir
    mkdir -p $InputDir
    ComputeRDF="yes"
    CenterOfMassMovies="yes"
    P=0.00001
    MakeCADSSInputWidom > $InputDir/simulation.input
    JOBName="$System"
    SLURMFile="run-CADSS.slurm"
    MakeSLURM > $CalcDir/$SLURMFile
    chmod +x $CalcDir/$SLURMFile

    cd $CalcDir
    sbatch $CalcDir/$SLURMFile
    exit
} # end widom 

[ "$SimType" == 'gcmc' ]  && { 
   Point=0
   for Px in 10; do
   #for Px in $(seq  0.001 0.0025 0.01) $(seq  0.01 0.1 1) $(seq 1 0.5 10) ; do
       # for T = 298 P0 = 0.0318, for T=303 K, P0 = 0.0423? or  0.042455  or 0.0425
       #P=$( echo "$Px * 1.0 " | bc -l)
       #P=$(echo "${Px} * 0.0425" | bc -l)
   
       #[ "$NeedMakeGrid" == "no" ] && [ $Point == 1 ] && continue
       ComputeRDF="no"
       CenterOfMassMovies="no"
       Point=$((Point+1))
       [ $Point == 1 ] || [ $Point == 3 ] ||  [ $Point == 12 ] || [ $Point == 22 ] || [ $Point == 24 ] || [ $Point == 26 ] || [ $Point == 28 ] || [ $Point == 38 ] &&  {
       ComputeRDF="yes"
       CenterOfMassMovies="yes"
       } 
   
       #P=$(printf "%12.10f"  $P)
       #P_in_MPa=$(echo "$P * 0.1"  | bc -l)
       #P_in_MPa=$(printf "%14.12f"  ${P_in_MPa})
   
       #Pi=Ptotal*xi
       #P1=$(echo "$P * $x1"  | bc -l)
       #P2=$(echo "$P * $x2"  | bc -l)
       #T=$(printf "%05.1f"  $T)
       P=$(printf "%08.05f"  $Px | tr ',' '.')
       P_in_MPa=$(awk "BEGIN {printf \"%.5f\", $P * 0.1}")
       echo "$Point  GCMC run for P = $P bar,  $Adsorbate1  RDF=$ComputeRDF " | tee -a $WorkingDir/log
   
       CalcDir=$WorkingDir/P$P
       InputDir=$CalcDir/input
       mkdir -p $CalcDir
       mkdir -p $InputDir
   
       MakeCADSSInputGCMC > $InputDir/simulation.input
   
       JOBName="$System-$P"
       SLURMFile="run-CADSS.slurm"
       MakeSLURM > $CalcDir/$SLURMFile
       chmod +x $CalcDir/$SLURMFile
   
       cd $CalcDir
       sbatch $CalcDir/$SLURMFile
       [ "$NeedMakeGrid" == "yes" ] && exit 
   
       cd $WorkingDir
   done
} # end gcmc 

[ "$SimType" == 'NVT' ]  && { 
   Point=0
   for nMolecules  in 1 ; do
   
       ComputeRDF="yes"
       CenterOfMassMovies="yes"

       [ $nMolecules == 2 ] &&  {
        CenterOfMassMovies="yes"
       }

       Point=$((Point+1))
       P=0.000001 
       #P=$(printf "%12.10f"  $P)
       P_in_MPa=$(echo "$P * 0.1"  | bc -l)
       P_in_MPa=$(printf "%14.12f"  ${P_in_MPa})
   
       echo "$Point  NVT  run for nMolecules=$nMolecules   RDF=$ComputeRDF CenterOfMassMovies=$CenterOfMassMovies "
   
       CalcDir=$WorkingDir/$System-$nMolecules-2E7
       InputDir=$CalcDir/input
       mkdir -p $CalcDir
       mkdir -p $InputDir
   
       MakeCADSSInputNVT > $InputDir/simulation.input
   
       JOBName="$System-$nMolecules"
       SLURMFile="run-CADSS.slurm"
       MakeSLURM > $CalcDir/$SLURMFile
       chmod +x $CalcDir/$SLURMFile
       cd $CalcDir 
       sbatch $CalcDir/$SLURMFile
       [ "$NeedMakeGrid" == "yes" ] && exit 
   
       cd $WorkingDir
   
   done
} # end NVT  

} # end main funtion 

function MakeSLURM () {
cat << EOF
#!/bin/bash
# CADSS-V3 in chesdens2/slurm
#SBATCH --job-name=${JOBName}
#SBATCH --nodes=1
#SBATCH --partition=damp,fun
#SBATCH --cpus-per-task=1
#SBATCH --exclude=node35
#SBATCH --output=job-%n-%A.out
#_SBATCH --time=120:00:00

export OMP_NUM_THREADS=1
cd \$SLURM_SUBMIT_DIR  
echo \$SLURM_JOBID > jobid  
echo \$SLURM_JOB_NODELIST >> jobid  

${CADSSBinary}  >  ${JOBName}.log 2>&1

exit 0

EOF
}

function MakeCADSSInputGCMC {
    cat << EOF 
#HT-CADSS version3.1_2021
#copyright            : (C) 2021 by Dr. Qingyuan Yang
#Unit of Varibales---Pressure[input]---Temperature[K]---Distance[Angstrom]---Time[ps] 

#-----------------General Simulation Information------------------------------------------
ScreenPurpose                       no
SimulationType                      MonteCarlo
Ensemble                            MuVT
ExternalTemperature_[K]             $T
ExternalPressure_Unit               $P bar
CutOff_[A]                          12.0
CutoffCoulomb_[A]                   12.0
RandomSeedOfGenerator               no
NumberOfComponents                  1
#------------------Equation of State Information (calculate fugacity)------------------------
UseEOStoCalculateFugacity           yes
NameOfEquationOfState               Peng_Robinson
ConsiderGasNonIdeality              yes
#-----------------Simulation Step Controlling Information----------------------------------
MonteCarloStepStyle                 InnerCycleNotUsed                          
NumberOfEquilibrationCycles         15E6
NumberOfProductionCycles            1E5
MinimumInnerCycles                  20
PrintEvery                          1E5
Output                              file     
#-----------------Simulation Cell Information---------------------------------------------
Framework                           ${System}
IonicTypeFramework                  no
InputFileType                       car
FlexibleFramework                   no
#-----------------Framework Charge Equilibration Information-------------------------------
FrameworkAtomChargeSource           	frameworkstructurefile
#---------------Ewald Summation Controlling Information-----------------------------------
ChargeMethod                        Ewald
ChargeFrameworkMethod               Ewald
EwaldAutomatic                      yes   
EwaldPrecision                      1e-6
#-----------------Force field Information-------------------------------------------------
Forcefield                          ${ForcefieldDir}
AdsorbateAdsorbateInteraction       yes
#-----------------Potential Grids Tabulating Information----------------------------------
UseTabularGrid                      ${UseTabularGrid}
NeedMakeGrid                        ${NeedMakeGrid}
SpacingVDWGrid_[A]                  0.15
SpacingCoulombGrid_[A]              0.15
#-------------------------------- Make Center of Mass trajectory ---------------------------
CenterOfMassMovies                  ${CenterOfMassMovies}
WriteComMoviesEvery                 2E5
#-------------- Radial Distribution Function Information -----------------------------------
ComputeRDF                          ${ComputeRDF}
ComputeRDFEvery                     1E4
PrintRDFEvery                       2E4
#---------------Components Change Probability in Monte Carlo Moves-------------------------
  component 1 MoleculeName                     ${Adsorbate1}
              MoleculeDefinition               ${MoleculeDefinition}
              BulkMolarFraction                1.0                     
              StartingBead                     1
              IdealGasRosenbluthWeight         1.0
              BlockPockets                     no
                BlockPocketsFilename           ${System}
              TranslationProbability           0.15
                 TranslationDirection          XYZ
              RotationProbability              0.15
              RegrowFullProbability            0.15
                 RegrowFullMethod              RandomGrowth
              RegrowPartialProbability         0.0
              SwapProbability                  0.55
                 SwapMethod                    RandomSwap
              IdentityChangeProbability        0.0
                 NumberOfIdentityChanges       0
                 IdentityChangesList           2
                 IndentityChangeMethod         RandomIdentityChange
              WidomProbability                 0.0
              ExtraFrameworkMolecule           no
              CreateNumberOfMolecules_Unit     0  unitcell
              UsingOldConfigurations           no

EOF
}

function MakeCADSSInputWidom {
    cat << EOF 
#HT-CADSS version3.1_2021
#copyright            : (C) 2021 by Dr. Qingyuan Yang
#Unit of Varibales---Pressure[input]---Temperature[K]---Distance[Angstrom]---Time[ps] 

#-----------------General Simulation Information------------------------------------------
ScreenPurpose                       no
SimulationType                      MonteCarlo
Ensemble                            MuVT
ExternalTemperature_[K]             $T
ExternalPressure_Unit               $P bar
CutOff_[A]                          12.0
CutoffCoulomb_[A]                   12.0
RandomSeedOfGenerator               no
NumberOfComponents                  1
#------------------Equation of State Information (calculate fugacity)------------------------
UseEOStoCalculateFugacity           yes
NameOfEquationOfState               Peng_Robinson
ConsiderGasNonIdeality              yes
#-----------------Simulation Step Controlling Information----------------------------------
MonteCarloStepStyle                 InnerCycleNotUsed                          
NumberOfEquilibrationCycles         1E5
NumberOfProductionCycles            2E5
MinimumInnerCycles                  20
PrintEvery                          1000
Output                              file     
#-----------------Simulation Cell Information---------------------------------------------
Framework                           ${System}
IonicTypeFramework                  no
InputFileType                       car
FlexibleFramework                   no
#-----------------Framework Charge Equilibration Information-------------------------------
FrameworkAtomChargeSource           forcefieldfile
#---------------Ewald Summation Controlling Information-----------------------------------
ChargeMethod                        Ewald
ChargeFrameworkMethod               Ewald
EwaldAutomatic                      yes   
EwaldPrecision                      1e-6
#-----------------Force field Information-------------------------------------------------
Forcefield                          ${ForcefieldDir}
AdsorbateAdsorbateInteraction       yes
#-----------------Potential Grids Tabulating Information----------------------------------
UseTabularGrid                      ${UseTabularGrid}
NeedMakeGrid                        ${NeedMakeGrid}
SpacingVDWGrid_[A]                  0.15
SpacingCoulombGrid_[A]              0.15
#-------------------------------- Make Center of Mass trajectory ---------------------------
CenterOfMassMovies                  ${CenterOfMassMovies}
WriteComMoviesEvery                 2E5
#-------------- Radial Distribution Function Information -----------------------------------
ComputeRDF                          ${ComputeRDF}
ComputeRDFEvery                     1E4
PrintRDFEvery                       2E4
#-------------- Poistion Historgram Information --------------------------------------------
ComputePositionHistograms           no
PositionHistogramMappingType        ABC
WritePositionHistogramEvery         1E4
#-------------- free energy profiles            --------------------------------------------
ComputeFreeEnergyProfileWPIMethod   yes
FreeEnergyProfileMappingType        ABC
WriteFreeEnergyProfileEvery         1E4
#-------------- PSD ------------------------------------------------------------------------
ComputePoreSizeDistribution         yes
PSDCubicGridSize_[A]                0.10
PSDMaxPoreSize_[A]                  10
PSDBinSize_[A]                      0.20
PSDNitrogenDiameter_[A]             3.681
PSDUseCSDVdwRadiiForFrameworkAtoms  yes
#----------------- Framework Surface Area Calculation Information --------------------------
ComputeSurfaceArea                  no
SurfaceAreaProbeAtomName            N2
SurfaceAreaSamplingPointsPerShere   200
SurfaceAreaProbeDiameter_[A]        3.681
SurfaceAreaProbeDistance            Sigma
#----------------- Framework Void Volume Calculation Information --------------------------
ComputeVoidVolume                   no
VoidVolumeProbeAtomName             VOL_He
HeliumVoidFraction                  0.81
#-------------- END to END Distance --------------------------------------------------------
ComputeEndToEndDistanceHistograms   no
ComputeEnergyHistograms             no
ComputeHistogramsReweightingGCMC    yes
#---------------Components Change Probability in Monte Carlo Moves-------------------------
  component 1 MoleculeName                     ${Adsorbate1}
              MoleculeDefinition               ${MoleculeDefinition}
              BulkMolarFraction                1.0                     
              IdealGasRosenbluthWeight         1.0
              BlockPockets                     no
                BlockPocketsFilename           ${System} 
              WidomProbability                 1.0

EOF
}

function MakeCADSSInputNVT {
    cat << EOF 
#HT-CADSS version3.1_2021
#copyright            : (C) 2021 by Dr. Qingyuan Yang
#Unit of Varibales---Pressure[input]---Temperature[K]---Distance[Angstrom]---Time[ps] 

#-----------------General Simulation Information------------------------------------------
ScreenPurpose                       no
SimulationType                      MonteCarlo
Ensemble                            NVT
ExternalTemperature_[K]             $T
ExternalPressure_Unit               $P bar
CutOff_[A]                          12.0
CutoffCoulomb_[A]                   12.0
RandomSeedOfGenerator               no
NumberOfComponents                  1
#------------------Equation of State Information (calculate fugacity)------------------------
UseEOStoCalculateFugacity           yes
NameOfEquationOfState               Peng_Robinson
ConsiderGasNonIdeality              yes
#-----------------Simulation Step Controlling Information----------------------------------
MonteCarloStepStyle                 InnerCycleNotUsed                          
NumberOfEquilibrationCycles         1E7
NumberOfProductionCycles            2E7
MinimumInnerCycles                  20
PrintEvery                          20000
Output                              file     
#-----------------Simulation Cell Information---------------------------------------------
Framework                           ${System}
IonicTypeFramework                  no
InputFileType                       car
FlexibleFramework                   no
#-----------------Framework Charge Equilibration Information-------------------------------
FrameworkAtomChargeSource           forcefieldfile
#---------------Ewald Summation Controlling Information-----------------------------------
ChargeMethod                        Ewald
ChargeFrameworkMethod               Ewald
EwaldAutomatic                      yes   
EwaldPrecision                      1e-6
#-----------------Force field Information-------------------------------------------------
Forcefield                          ${ForcefieldDir}
AdsorbateAdsorbateInteraction       yes
#-----------------Potential Grids Tabulating Information----------------------------------
UseTabularGrid                      ${UseTabularGrid}
NeedMakeGrid                        ${NeedMakeGrid}
SpacingVDWGrid_[A]                  0.15
SpacingCoulombGrid_[A]              0.15
#-------------------------------- Make Center of Mass trajectory ---------------------------
CenterOfMassMovies                  ${CenterOfMassMovies}
WriteComMoviesEvery                 2E5
#-------------- Radial Distribution Function Information -----------------------------------
ComputeRDF                          ${ComputeRDF}
ComputeRDFEvery                     1E4
PrintRDFEvery                       2E4
#-------------- END to END Distance --------------------------------------------------------
ComputeEndToEndDistanceHistograms   no
ComputeEnergyHistograms             no
ComputeHistogramsReweightingGCMC    yes
#---------------Components Change Probability in Monte Carlo Moves-------------------------
  component 1 MoleculeName                     ${Adsorbate1}
              MoleculeDefinition               ${MoleculeDefinition}
              BulkMolarFraction                1.0                     
              StartingBead                     1
              IdealGasRosenbluthWeight         1.0
              BlockPockets                     no
                BlockPocketsFilename           ${System}
              TranslationProbability           0.5
                 TranslationDirection          XYZ
              RotationProbability              0.5
              ExtraFrameworkMolecule           no
              CreateNumberOfMolecules_Unit     ${nMolecules}  box
              UsingOldConfigurations           no

EOF
}


#script executation starts...  
main 
exit

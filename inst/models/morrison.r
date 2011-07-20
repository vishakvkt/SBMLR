model=list( 

id="MorrisonAllegra",

notes=c(
"This is a folate model that includes MTX polyglutamation.", 
"Morrison and Allegra, JBC:264,10552-10566 (1989)",
"Folate cycle kinetics in breast cancer cells.",
"o                                             ",
"Two flow BCs were converted into two downstream",
"concentration BCs, thus removing the GAR and dUMP state variables.",
"This dropped the number of ODEs from 21 to 19.",
"o                                               ",
"The subscript b following some species names implies bound to DHFR.",
"FH2 free and bound to DHFR are in rapid equilibrium.",
"DHFRsyn(thesis) parameters were picked to get SS concentration ",
"right with kdeg=.03 and also to get 1uM MTX to ramp DHFR up by a factor of 3.3.",
"MTX1 efflux is strong enough that GGH21 overwhelms MTX2 export.",
"Keq is a global parameter; this wouldn't be needed if rules had local parameters!",
"DHFR protein degradation reactions are independent of folate binding.",
"In DHFRdeg, since mass cannot come out of the BC it was taken form free DHFR",
"In FH2bdeg we create back the FH2 missing from DHFR protein degradation.",
"The units in this model are concentrations in micromolar and time in hours.",
"Volume sizes of 1 were picked to equate concentrations to mass or moles"
),

compartments=list(
list(id="cell", size = 1),
list(id="ext", size = 1)
),

species=list(
list(id="FH2f"		,ic=  0.0012,	 compartment="cell",	 bc=FALSE),
list(id="FH2b"		,ic=  0.0024,	 compartment="cell",	 bc=TRUE),
list(id="DHFRf"		,ic=    0.64,	 compartment="cell",	 bc=FALSE),
list(id="DHFRtot"	,ic=    0.64,	 compartment="cell",	 bc=TRUE),
list(id="FH4"		,ic=    0.46,	 compartment="cell",	 bc=FALSE),
list(id="CH2FH4"	,ic=    0.26,	 compartment="cell",	 bc=FALSE),
list(id="CH3FH4"	,ic=    1.63,	 compartment="cell",	 bc=FALSE),
list(id="CHOFH4"	,ic=       1,	 compartment="cell",	 bc=FALSE),
list(id="FFH2"		,ic=0.000332,	 compartment="cell",	 bc=FALSE),
list(id="HCHO"		,ic=  0.0074,	 compartment="cell",	 bc=FALSE),
list(id="FGAR"		,ic=   16.49,	 compartment="cell",	 bc=FALSE),
list(id="AICAR"		,ic=   3.695,	 compartment="cell",	 bc=FALSE),
list(id="MTX1"		,ic=       0,	 compartment="cell",	 bc=FALSE),
list(id="MTX2"		,ic=       0,	 compartment="cell",	 bc=FALSE),
list(id="MTX3"		,ic=       0,	 compartment="cell",	 bc=FALSE),
list(id="MTX4"		,ic=       0,	 compartment="cell",	 bc=FALSE),
list(id="MTX5"		,ic=       0,	 compartment="cell",	 bc=FALSE),
list(id="MTX1b"		,ic=       0,	 compartment="cell",	 bc=FALSE),
list(id="MTX2b"		,ic=       0,	 compartment="cell",	 bc=FALSE),
list(id="MTX3b"		,ic=       0,	 compartment="cell",	 bc=FALSE),
list(id="MTX4b"		,ic=       0,	 compartment="cell",	 bc=FALSE),
list(id="MTX5b"		,ic=       0,	 compartment="cell",	 bc=FALSE),
list(id="EMTX"		,ic=       0,	 compartment="ext",	 bc=TRUE),
list(id="dUMP"		,ic=   20.76,	 compartment="cell",	 bc=TRUE),
list(id="GAR"		,ic=   689.6,	 compartment="cell",	 bc=TRUE),
list(id="serine"	,ic=   123.3,	 compartment="cell",	 bc=TRUE),
list(id="formate"	,ic=     500,	 compartment="cell",	 bc=TRUE),
list(id="ATP"		,ic=    2980,	 compartment="cell",	 bc=TRUE),
list(id="glutamine"	,ic=    7170,	 compartment="cell",	 bc=TRUE),
list(id="glycine"	,ic=    1600,	 compartment="cell",	 bc=TRUE),
list(id="NADP"		,ic=    6.73,	 compartment="cell",	 bc=TRUE),
list(id="NADPH"		,ic=     294,	 compartment="cell",	 bc=TRUE),
list(id="homocysteine"	,ic=      10,	 compartment="cell",	 bc=TRUE)
),

globalParameters=list(Keq=0.32), 

rules=list(
list(idOutput=c("FH2b"),
inputs=c("FH2f","DHFRf"),
strLaw="FH2f*DHFRf/Keq"
),

list(idOutput=c("DHFRtot"),
inputs=c("FH2b","DHFRf","MTX1b","MTX2b","MTX3b","MTX4b","MTX5b"),
strLaw="FH2b+DHFRf+MTX1b+MTX2b+MTX3b+MTX4b+MTX5b"
)),
   

reactions=list(
list( id="SHMT", reversible=FALSE,
reactants=c("FH4","serine"),
products=c("CH2FH4"),
parameters=c(Vm = 18330,Km1 = 1.7,Km2 = 210),
strLaw="Vm*(serine/Km2/(1+serine/Km2))*(FH4/Km1/(1+FH4/Km1))"
),

list( id="SHMTr", reversible=FALSE,
reactants=c("CH2FH4"),
modifiers=c("glycine"),
products=c("FH4"),
parameters=c(Vm = 1.22e+007,Km1 = 3200,Km2 = 10000),
strLaw="Vm*(glycine/Km2/(1+glycine/Km2))*(CH2FH4/Km1/(1+CH2FH4/Km1))"
),

list( id="HCHOtoCH2FH4", reversible=FALSE,
reactants=c("FH4","HCHO"),
products=c("CH2FH4"),
parameters=c(hp = 23.2),
strLaw="hp*FH4*HCHO"
),

list( id="CH2FH4toHCHO", reversible=FALSE,
reactants=c("CH2FH4"),
products=c("FH4","HCHO"),
parameters=c(hl = 0.3),
strLaw="hl*CH2FH4"
),

list( id="MTHFR", reversible=FALSE,
reactants=c("CH2FH4","NADPH"),
modifiers=c("FH2f","MTX1","MTX2","MTX3","MTX4","MTX5"),
products=c("CH3FH4"),
parameters=c(Vm = 224.8,Km1 = 50,Km2 = 50,Ki1 = 0.4,Ki21 = 59,Ki22 = 21.3,Ki23 = 7.68,Ki24 = 2.77,Ki25 = 1),
strLaw="Vm*(NADPH/Km2/(1+NADPH/Km2))*(CH2FH4/(CH2FH4+Km1*(1+MTX1/Ki21+MTX2/Ki22+MTX3/Ki23+MTX4/Ki24+MTX5/Ki25+FH2f/Ki1)))"
),

list( id="MTR", reversible=FALSE,
reactants=c("CH3FH4","homocysteine"),
products=c("FH4"),
parameters=c(Vm = 22600,Km1 = 125,Km2 = 2900),
strLaw="Vm*(homocysteine/Km2/(1+homocysteine/Km2))*(CH3FH4/Km1/(1+CH3FH4/Km1))"
),

list( id="HCOOHtoCHOFH4", reversible=FALSE,
reactants=c("FH4","formate","ATP"),
products=c("CHOFH4"),
parameters=c(Vm = 3600,Km1 = 230,Km2 = 56,Km3 = 1600),
strLaw="Vm*(FH4/Km1/(1+FH4/Km1))*(ATP/Km2/(1+ATP/Km2))*(formate/Km3/(1+formate/Km3))"
),

list( id="GARFT", reversible=FALSE,
reactants=c("CHOFH4","GAR"),
modifiers=c("FH2f","FFH2","MTX1","MTX2","MTX3","MTX4","MTX5"),
products=c("FGAR","FH4"),
parameters=c(Vm = 4126,Km1 = 4.9,Km2 = 52,Ki1 = 5,Ki1f = 1,Ki21 = 84,Ki22 = 60,Ki23 = 43,Ki24 = 31,Ki25 = 22),
strLaw="Vm*(GAR/Km2/(1+GAR/Km2))*(CHOFH4/(CHOFH4+Km1*(1+MTX1/Ki21+MTX2/Ki22+MTX3/Ki23+MTX4/Ki24+MTX5/Ki25+FH2f/Ki1+FFH2/Ki1f)))"
),

list( id="ATIC7", reversible=FALSE,
reactants=c("CHOFH4","AICAR"),
modifiers=c("FH2f","FFH2","MTX1","MTX2","MTX3","MTX4","MTX5"),
products=c("FH4"),
parameters=c(Vm = 31675,Km1 = 5.5,Km2 = 24,Ki1 = 2.89,Ki1f = 5.3,Ki21 = 40,Ki22 = 31.5,Ki23 = 2.33,Ki24 = 3.61,Ki25 = 5.89),
strLaw="Vm*(AICAR/Km2/(1+AICAR/Km2))*(CHOFH4/(CHOFH4+Km1*(1+MTX1/Ki21+MTX2/Ki22+MTX3/Ki23+MTX4/Ki24+MTX5/Ki25+FH2f/Ki1+FFH2/Ki1f)))"
),

list( id="MTHFD", reversible=FALSE,
reactants=c("CH2FH4","NADP"),
products=c("CHOFH4"),
parameters=c(Vm = 68500,Km1 = 3,Km2 = 21.8),
strLaw="Vm*(CH2FH4/Km1/(1+CH2FH4/Km1))*(NADP/Km2/(1+NADP/Km2))"
),

list( id="TYMS", reversible=FALSE,
reactants=c("CH2FH4","dUMP"),
modifiers=c("FH2f","FFH2","MTX1","MTX2","MTX3","MTX4","MTX5"),
products=c("FH2f"),
parameters=c(Vm = 58,Km1 = 2.5,Km2 = 1.8,Ki1 = 3,Ki1f = 1.6,Ki21 = 13,Ki22 = 0.08,Ki23 = 0.07,Ki24 = 0.065,Ki25 = 0.047),
strLaw="Vm*(CH2FH4/Km1)*(dUMP/Km2)/((CH2FH4/Km1)*(dUMP/Km2)*(1+MTX1/Ki21+MTX2/Ki22+MTX3/Ki23+MTX4/Ki24+MTX5/Ki25+FH2f/Ki1)+(dUMP/Km2)*((FFH2/Ki1f)*(MTX1/Ki21)+(1+FFH2/Ki1f)*(1+MTX2/Ki22+MTX3/Ki23+MTX4/Ki24+MTX5/Ki25+FH2f/Ki1))+(1+MTX2/Ki22+MTX3/Ki23+MTX4/Ki24+MTX5/Ki25+FH2f/Ki1))"
),

list( id="DHFReductase", reversible=FALSE,
reactants=c("FH2f"),
modifiers=c("FH2b"),
products=c("FH4"),
parameters=c(kter = 2109.4),
strLaw="kter*FH2b"
),

list( id="FFH2syn", reversible=FALSE,
reactants=c("FH2f"),
products=c("FFH2"),
parameters=c(Vm = 65),
strLaw="Vm*FH2f"
),

list( id="ATIC12", reversible=FALSE,
reactants=c("FFH2","AICAR"),
modifiers=c("CHOFH4","FH2f","MTX1","MTX2","MTX3","MTX4","MTX5"),
products=c("FH2f"),
parameters=c(Vm = 9503,Km1 = 5.3,Km2 = 24,Ki1 = 2.89,Ki1f = 5.5,Ki21 = 40,Ki22 = 31.5,Ki23 = 2.33,Ki24 = 3.61,Ki25 = 5.89),
strLaw="Vm*(AICAR/Km2/(1+AICAR/Km2))*(FFH2/(FFH2+Km1*(1+MTX1/Ki21+MTX2/Ki22+MTX3/Ki23+MTX4/Ki24+MTX5/Ki25+FH2f/Ki1+CHOFH4/Ki1f)))"
),

list( id="AICARsyn", reversible=FALSE,
reactants=c("FGAR"),
modifiers=c("glutamine"),
products=c("AICAR"),
parameters=c(Vm = 4656,Km1 = 100,Km2 = 100),
strLaw="Vm*(glutamine/Km1/(1+glutamine/Km1))*(FGAR/Km2/(1+FGAR/Km2))"
),

list( id="FPGS12", reversible=FALSE,
reactants=c("MTX1"),
products=c("MTX2"),
parameters=c(Vm = 0.129),
strLaw="Vm*MTX1"
),

list( id="FPGS23", reversible=FALSE,
reactants=c("MTX2"),
products=c("MTX3"),
parameters=c(Vm = 0.369),
strLaw="Vm*MTX2"
),

list( id="FPGS34", reversible=FALSE,
reactants=c("MTX3"),
products=c("MTX4"),
parameters=c(Vm = 0.118),
strLaw="Vm*MTX3"
),

list( id="FPGS45", reversible=FALSE,
reactants=c("MTX4"),
products=c("MTX5"),
parameters=c(Vm = 0.185),
strLaw="Vm*MTX4"
),

list( id="GGH21", reversible=FALSE,
reactants=c("MTX2"),
products=c("MTX1"),
parameters=c(Vm = 0.195),
strLaw="Vm*MTX2"
),

list( id="GGH32", reversible=FALSE,
reactants=c("MTX3"),
products=c("MTX2"),
parameters=c(Vm = 0.025),
strLaw="Vm*MTX3"
),

list( id="GGH43", reversible=FALSE,
reactants=c("MTX4"),
products=c("MTX3"),
parameters=c(Vm = 0.031),
strLaw="Vm*MTX4"
),

list( id="GGH54", reversible=FALSE,
reactants=c("MTX5"),
products=c("MTX4"),
parameters=c(Vm = 0.191),
strLaw="Vm*MTX5"
),

list( id="RFC", reversible=FALSE,
reactants=c("EMTX"),
products=c("MTX1"),
parameters=c(Vm = 82.2,Km = 8.2),
strLaw="Vm*EMTX/Km/(1+EMTX/Km)"
),

list( id="MTX1export", reversible=FALSE,
reactants=c("MTX1"),
parameters=c(Vm = 4.65),
strLaw="Vm*MTX1"
),

list( id="MTX2export", reversible=FALSE,
reactants=c("MTX2"),
parameters=c(Vm = 0),
strLaw="Vm*MTX2"
),

list( id="MTX3export", reversible=FALSE,
reactants=c("MTX3"),
parameters=c(Vm = 0.063),
strLaw="Vm*MTX3"
),

list( id="MTX4export", reversible=FALSE,
reactants=c("MTX4"),
parameters=c(Vm = 0.063),
strLaw="Vm*MTX4"
),

list( id="MTX5export", reversible=FALSE,
reactants=c("MTX5"),
parameters=c(Vm = 0.063),
strLaw="Vm*MTX5"
),

list( id="MTX1on", reversible=FALSE,
reactants=c("MTX1","DHFRf"),
products=c("MTX1b"),
parameters=c(Vm = 23100),
strLaw="Vm*DHFRf*MTX1"
),

list( id="MTX2on", reversible=FALSE,
reactants=c("MTX2","DHFRf"),
products=c("MTX2b"),
parameters=c(Vm = 44300),
strLaw="Vm*DHFRf*MTX2"
),

list( id="MTX3on", reversible=FALSE,
reactants=c("MTX3","DHFRf"),
products=c("MTX3b"),
parameters=c(Vm = 85100),
strLaw="Vm*DHFRf*MTX3"
),

list( id="MTX4on", reversible=FALSE,
reactants=c("MTX4","DHFRf"),
products=c("MTX4b"),
parameters=c(Vm = 163000),
strLaw="Vm*DHFRf*MTX4"
),

list( id="MTX5on", reversible=FALSE,
reactants=c("MTX5","DHFRf"),
products=c("MTX5b"),
parameters=c(Vm = 314000),
strLaw="Vm*DHFRf*MTX5"
),

list( id="MTX1off", reversible=FALSE,
reactants=c("MTX1b"),
products=c("MTX1","DHFRf"),
parameters=c(Vm = 0.42),
strLaw="Vm*MTX1b"
),

list( id="MTX2off", reversible=FALSE,
reactants=c("MTX2b"),
products=c("MTX2","DHFRf"),
parameters=c(Vm = 0.42),
strLaw="Vm*MTX2b"
),

list( id="MTX3off", reversible=FALSE,
reactants=c("MTX3b"),
products=c("MTX3","DHFRf"),
parameters=c(Vm = 0.42),
strLaw="Vm*MTX3b"
),

list( id="MTX4off", reversible=FALSE,
reactants=c("MTX4b"),
products=c("MTX4","DHFRf"),
parameters=c(Vm = 0.42),
strLaw="Vm*MTX4b"
),

list( id="MTX5off", reversible=FALSE,
reactants=c("MTX5b"),
products=c("MTX5","DHFRf"),
parameters=c(Vm = 0.42),
strLaw="Vm*MTX5b"
),

list( id="DHFRfsyn", reversible=FALSE,
modifiers=c("EMTX"),
products=c("DHFRf"),
parameters=c(k0 = 0.0192,k1 = 0.04416),
strLaw="k0+k1*EMTX"
),

list( id="DHFRdeg", reversible=FALSE,
reactants=c("DHFRf"),
modifiers=c("FH2b"),
parameters=c(Vm = 0.03),
strLaw="Vm*(DHFRf+FH2b)"
),

list( id="FH2bdeg", reversible=FALSE,
modifiers=c("FH2b"),
products=c("FH2f"),
parameters=c(Vm = 0.03),
strLaw="Vm*FH2b"
),

list( id="MTX1deg", reversible=FALSE,
reactants=c("MTX1b"),
products=c("MTX1"),
parameters=c(Vm = 0.03),
strLaw="Vm*MTX1b"
),

list( id="MTX2deg", reversible=FALSE,
reactants=c("MTX2b"),
products=c("MTX2"),
parameters=c(Vm = 0.03),
strLaw="Vm*MTX2b"
),

list( id="MTX3deg", reversible=FALSE,
reactants=c("MTX3b"),
products=c("MTX3"),
parameters=c(Vm = 0.03),
strLaw="Vm*MTX3b"
),

list( id="MTX4deg", reversible=FALSE,
reactants=c("MTX4b"),
products=c("MTX4"),
parameters=c(Vm = 0.03),
strLaw="Vm*MTX4b"
),

list( id="MTX5deg", reversible=FALSE,
reactants=c("MTX5b"),
products=c("MTX5"),
parameters=c(Vm = 0.03),
strLaw="Vm*MTX5b"
)) 

)  # end model

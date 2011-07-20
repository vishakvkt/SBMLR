model=list( 

id="curto",

notes=c(
"This is a purine metabolism model that is geared toward studies of gout.", 
"The model is fully described in Curto et al., MBSC 151 (1998) pp 1-49",
"The model uses Generalized Mass Action (GMA;i.e. power law) descriptions of reaction rate laws.",
"Such descriptions are local approximations that assume independent substrate binding. ",
"                    ",
"The de novo purine flux vden= 2.39 is in umole/min/KG, i.e. 2.4*60=144 uM/h if we let each Kg be a",
"liter of water. Morrison and Allegra (JBC, 1989) have vden at 650 uM/h (model) and 415 (exp).",
"The IC's below have been set to the system's steady state.",
"The units in this model are micromolar(uM) and minutes.",
"A cell volume of 1 is used so that amounts and concentrations are the same thing."
),

compartments=list(
list(id="cell", size = 1)
),

species=list(
list(id="PRPP"		,ic= 5.01742,	 compartment="cell",	 bc=FALSE),
list(id="IMP"		,ic= 98.2634,	 compartment="cell",	 bc=FALSE),
list(id="SAMP"		,ic=0.198189,	 compartment="cell",	 bc=FALSE),
list(id="ATP"		,ic= 2475.35,	 compartment="cell",	 bc=FALSE),
list(id="SAM"		,ic= 3.99187,	 compartment="cell",	 bc=FALSE),
list(id="Ade"		,ic= 0.98473,	 compartment="cell",	 bc=FALSE),
list(id="XMP"		,ic=  24.793,	 compartment="cell",	 bc=FALSE),
list(id="GTP"		,ic= 410.223,	 compartment="cell",	 bc=FALSE),
list(id="dATP"		,ic= 6.01413,	 compartment="cell",	 bc=FALSE),
list(id="dGTP"		,ic= 3.02581,	 compartment="cell",	 bc=FALSE),
list(id="RNA"		,ic= 28680.5,	 compartment="cell",	 bc=FALSE),
list(id="DNA"		,ic= 5179.34,	 compartment="cell",	 bc=FALSE),
list(id="HX"		,ic= 9.51785,	 compartment="cell",	 bc=FALSE),
list(id="Xa"		,ic= 5.05941,	 compartment="cell",	 bc=FALSE),
list(id="Gua"		,ic= 5.50638,	 compartment="cell",	 bc=FALSE),
list(id="UA"		,ic= 100.293,	 compartment="cell",	 bc=FALSE),
list(id="R5P"		,ic=      18,	 compartment="cell",	 bc=TRUE),
list(id="Pi"		,ic=    1400,	 compartment="cell",	 bc=TRUE)
),

reactions=list(
list( id="ada", reversible=FALSE,
reactants=c("ATP"),
products=c("HX"),
parameters=c(aada = 0.001062,fada4 = 0.97),
strLaw="aada*ATP^fada4"
),

list( id="ade", reversible=FALSE,
reactants=c("Ade"),
parameters=c(aade = 0.01,fade6 = 0.55),
strLaw="aade*Ade^fade6"
),

list( id="adna", reversible=FALSE,
reactants=c("dATP"),
modifiers=c("dGTP"),
products=c("DNA"),
parameters=c(aadna = 3.2789,fdnap9 = 0.42,fdnap10 = 0.33),
strLaw="aadna*dATP^fdnap9*dGTP^fdnap10"
),

list( id="adrnr", reversible=FALSE,
reactants=c("ATP"),
modifiers=c("dGTP","dATP"),
products=c("dATP"),
parameters=c(aadrnr = 0.0602,fadrnr4 = 0.1,fadrnr9 = -0.3,fadrnr10 = 0.87),
strLaw="aadrnr*ATP^fadrnr4*dATP^fadrnr9*dGTP^fadrnr10"
),

list( id="ampd", reversible=FALSE,
reactants=c("ATP"),
modifiers=c("GTP","Pi"),
products=c("IMP"),
parameters=c(aampd = 0.02688,fampd4 = 0.8,fampd8 = -0.03,fampd18 = -0.1),
strLaw="aampd*ATP^fampd4*GTP^fampd8*Pi^fampd18"
),

list( id="aprt", reversible=FALSE,
reactants=c("PRPP","Ade"),
modifiers=c("ATP"),
products=c("ATP"),
parameters=c(aaprt = 233.8,faprt1 = 0.5,faprt4 = -0.8,faprt6 = 0.75),
strLaw="aaprt*PRPP^faprt1*ATP^faprt4*Ade^faprt6"
),

list( id="arna", reversible=FALSE,
reactants=c("ATP"),
modifiers=c("GTP"),
products=c("RNA"),
parameters=c(aarna = 614.5,frnap4 = 0.05,frnap8 = 0.13),
strLaw="aarna*ATP^frnap4*GTP^frnap8"
),

list( id="asuc", reversible=FALSE,
reactants=c("IMP"),
modifiers=c("ATP","GTP","Pi"),
products=c("SAMP"),
parameters=c(aasuc = 3.5932,fasuc2 = 0.4,fasuc4 = -0.24,fasuc8 = 0.2,fasuc18 = -0.05),
strLaw="aasuc*IMP^fasuc2*ATP^fasuc4*GTP^fasuc8*Pi^fasuc18"
),

list( id="asli", reversible=FALSE,
reactants=c("SAMP"),
modifiers=c("ATP"),
products=c("ATP"),
parameters=c(aasli = 66544,fasli3 = 0.99,fasli4 = -0.95),
strLaw="aasli*SAMP^fasli3*ATP^fasli4"
),

list( id="dada", reversible=FALSE,
reactants=c("dATP"),
products=c("HX"),
parameters=c(adada = 0.03333,fdada9 = 1),
strLaw="adada*dATP^fdada9"
),

list( id="den", reversible=FALSE,
reactants=c("PRPP"),
modifiers=c("dGTP","IMP","ATP","GTP","Pi"),
products=c("IMP"),
parameters=c(aden = 5.2728,fden1 = 2,fden2 = -0.06,fden4 = -0.25,fden8 = -0.2,fden18 = -0.08),
strLaw="aden*PRPP^fden1*IMP^fden2*ATP^fden4*GTP^fden8*Pi^fden18"
),

list( id="dgnuc", reversible=FALSE,
reactants=c("dGTP"),
products=c("Gua"),
parameters=c(adgnuc = 0.03333,fdgnuc10 = 1),
strLaw="adgnuc*dGTP^fdgnuc10"
),

list( id="dnaa", reversible=FALSE,
reactants=c("DNA"),
products=c("dATP"),
parameters=c(adnaa = 0.001938,fdnan12 = 1),
strLaw="adnaa*DNA^fdnan12"
),

list( id="dnag", reversible=FALSE,
reactants=c("DNA"),
products=c("dGTP"),
parameters=c(adnag = 0.001318,fdnan12 = 1),
strLaw="adnag*DNA^fdnan12"
),

list( id="gdna", reversible=FALSE,
reactants=c("dGTP"),
modifiers=c("dATP"),
products=c("DNA"),
parameters=c(agdna = 2.2296,fdnap9 = 0.42,fdnap10 = 0.33),
strLaw="agdna*dATP^fdnap9*dGTP^fdnap10"
),

list( id="gdrnr", reversible=FALSE,
reactants=c("GTP"),
modifiers=c("dATP","dGTP"),
products=c("dGTP"),
parameters=c(agdrnr = 0.1199,fgdrnr8 = 0.4,fgdrnr9 = -1.2,fgdrnr10 = -0.39),
strLaw="agdrnr*GTP^fgdrnr8*dATP^fgdrnr9*dGTP^fgdrnr10"
),

list( id="gmpr", reversible=FALSE,
reactants=c("GTP"),
modifiers=c("XMP","ATP","IMP"),
products=c("IMP"),
parameters=c(agmpr = 0.3005,fgmpr2 = -0.15,fgmpr4 = -0.07,fgmpr7 = -0.76,fgmpr8 = 0.7),
strLaw="agmpr*IMP^fgmpr2*ATP^fgmpr4*XMP^fgmpr7*GTP^fgmpr8"
),

list( id="gmps", reversible=FALSE,
reactants=c("XMP"),
modifiers=c("ATP"),
products=c("GTP"),
parameters=c(agmps = 0.3738,fgmps4 = 0.12,fgmps7 = 0.16),
strLaw="agmps*ATP^fgmps4*XMP^fgmps7"
),

list( id="gnuc", reversible=FALSE,
reactants=c("GTP"),
modifiers=c("Pi"),
products=c("Gua"),
parameters=c(agnuc = 0.2511,fgnuc8 = 0.9,fgnuc18 = -0.34),
strLaw="agnuc*GTP^fgnuc8*Pi^fgnuc18"
),

list( id="gprt", reversible=FALSE,
reactants=c("Gua","PRPP"),
modifiers=c("GTP"),
products=c("GTP"),
parameters=c(agprt = 361.69,fgprt1 = 1.2,fgprt8 = -1.2,fgprt15 = 0.42),
strLaw="agprt*PRPP^fgprt1*GTP^fgprt8*Gua^fgprt15"
),

list( id="grna", reversible=FALSE,
reactants=c("GTP"),
modifiers=c("ATP"),
products=c("RNA"),
parameters=c(agrna = 409.6,frnap4 = 0.05,frnap8 = 0.13),
strLaw="agrna*ATP^frnap4*GTP^frnap8"
),

list( id="gua", reversible=FALSE,
reactants=c("Gua"),
products=c("Xa"),
parameters=c(agua = 0.4919,fgua15 = 0.5),
strLaw="agua*Gua^fgua15"
),

list( id="hprt", reversible=FALSE,
reactants=c("HX","PRPP"),
modifiers=c("IMP"),
products=c("IMP"),
parameters=c(ahprt = 12.569,fhprt1 = 1.1,fhprt2 = -0.89,fhprt13 = 0.48),
strLaw="ahprt*PRPP^fhprt1*IMP^fhprt2*HX^fhprt13"
),

list( id="hx", reversible=FALSE,
reactants=c("HX"),
parameters=c(ahx = 0.003793,fhx13 = 1.12),
strLaw="ahx*HX^fhx13"
),

list( id="hxd", reversible=FALSE,
reactants=c("HX"),
products=c("Xa"),
parameters=c(ahxd = 0.2754,fhxd13 = 0.65),
strLaw="ahxd*HX^fhxd13"
),

list( id="impd", reversible=FALSE,
reactants=c("IMP"),
modifiers=c("GTP","XMP"),
products=c("XMP"),
parameters=c(aimpd = 1.2823,fimpd2 = 0.15,fimpd7 = -0.09,fimpd8 = -0.03),
strLaw="aimpd*IMP^fimpd2*XMP^fimpd7*GTP^fimpd8"
),

list( id="inuc", reversible=FALSE,
reactants=c("IMP"),
modifiers=c("Pi"),
products=c("HX"),
parameters=c(ainuc = 0.9135,finuc2 = 0.8,finuc18 = -0.36),
strLaw="ainuc*IMP^finuc2*Pi^finuc18"
),

list( id="mat", reversible=FALSE,
reactants=c("ATP"),
modifiers=c("SAM"),
products=c("SAM"),
parameters=c(amat = 7.2067,fmat4 = 0.2,fmat5 = -0.6),
strLaw="amat*ATP^fmat4*SAM^fmat5"
),

list( id="polyam", reversible=FALSE,
reactants=c("SAM"),
products=c("Ade"),
parameters=c(apolyam = 0.29,fpolyam5 = 0.9),
strLaw="apolyam*SAM^fpolyam5"
),

list( id="prpps", reversible=FALSE,
reactants=c("R5P"),
modifiers=c("ATP","GTP","Pi","PRPP"),
products=c("PRPP"),
parameters=c(aprpps = 0.9,fprpps1 = -0.03,fprpps4 = -0.45,fprpps8 = -0.04,fprpps17 = 0.65,fprpps18 = 0.7),
strLaw="aprpps*PRPP^fprpps1*ATP^fprpps4*GTP^fprpps8*R5P^fprpps17*Pi^fprpps18"
),

list( id="pyr", reversible=FALSE,
reactants=c("PRPP"),
parameters=c(apyr = 1.2951,fpyr1 = 1.27),
strLaw="apyr*PRPP^fpyr1"
),

list( id="rnaa", reversible=FALSE,
reactants=c("RNA"),
products=c("ATP"),
parameters=c(arnaa = 0.06923,frnan11 = 1),
strLaw="arnaa*RNA^frnan11"
),

list( id="rnag", reversible=FALSE,
reactants=c("RNA"),
products=c("GTP"),
parameters=c(arnag = 0.04615,frnan11 = 1),
strLaw="arnag*RNA^frnan11"
),

list( id="trans", reversible=FALSE,
reactants=c("SAM"),
products=c("ATP"),
parameters=c(atrans = 8.8539,ftrans5 = 0.33),
strLaw="atrans*SAM^ftrans5"
),

list( id="ua", reversible=FALSE,
reactants=c("UA"),
parameters=c(aua = 8.744e-005,fua16 = 2.21),
strLaw="aua*UA^fua16"
),

list( id="x", reversible=FALSE,
reactants=c("Xa"),
parameters=c(ax = 0.0012,fx14 = 2),
strLaw="ax*Xa^fx14"
),

list( id="xd", reversible=FALSE,
reactants=c("Xa"),
products=c("UA"),
parameters=c(axd = 0.949,fxd14 = 0.55),
strLaw="axd*Xa^fxd14"
))

)  

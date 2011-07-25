
After installing XML from http://www.omegahat.org/RSXML, 
an error message from library(XML)can be resolved 
by copying the *.dll files of the XML package 
(e.g. in C:\Program Files\R\rw2000\library\XML\libs) 
to the C:\windows directory (where Windows can find them). 

version changes
1.00 	initial release (11/05/2004)

1.01 	initial release trivially fixed to make SBML code pass through 
	the SBML.org validator by switching the order of prods and mods 

1.10 	substantive upgrade of fderiv to handle steady state and transient 
	microarray data perturbations (11/27/2004)

1.12 	Added the BMC Cancer 2004 example directory and removed similar 
	preliminary scripts from the demo directory of 1.10 (1/04/2005). 
     	Manual.doc was shrunken down to avoid redundancy with publication-based documentation

1.15 	Upgraded manual.doc to include a new fderiv description since the initial 
	release publication documentation was outdated.(1/17/2005) 

1.16 	Minor further upgrade of manual.doc to include figure of functions and objects (2/1/2005). 

1.20    Major upgrade. Model object of class SBML defined and summary and == methods defined for it. 
	Read and writes to files are now to and from this SBML model object. 
        The R model definition file has been maintained for model editing, but it has been streamlined.
        For example, functions definitions now only need the bottom line string expression. Conversions 
        to and from SBML and SBMLR are no longer lossy in appeal since expressions are used instead of 
        string with too many parentheses. Many new functions defined to make script writing clear. (5/4/2005)
        
1.21    Cleaning. The model field names rxns, prods, mods, reacts => reactions, products, modifiers and reactants.
	Better names for functions by starting them with action verb, e.g. readSBML, readSBMLR, saveSBML, saveSBMLR.  
	getIncidenceMatrix was incoporated into getModelInfo. (5/5/2005)

1.22    Simulate fixed to have vector rather than matrix passed for mod=1 with "Control" column now automatic for t<0. Simulate
	was also modified to take an initial state vector override of the one in the model. (5/28/2005)
	
1.25.1	getModelInfo was replaced by an expanded summary method output. The manual and BMCcancerFolates were updated (10/17/2005)

1.25.2	globalParameters now attached in summary method to sync changes in parameters in reactions (10/30/2005)

1.31.0  converted == method on SBML object to equateModels(model1,model2) function to remove check warnings (9/14/2007) 

1.37.1  Added namespaces file (11/20/2008)

1.42.0 Fixed makeLaw example to take names of p. Replaced recursive function calls in saveSBML with saveXML calls.    
       Removed species renaming block in readSBML (since already done in eventHandler), simplified readSBMLR species 
       block (and readSBML rule and reaction name blocks) using sapply rather than a for loop. 
       Used StatET's autoindent feature to format the codes better. Moved makeLaw() into end of readSBML.R to
       avoid the need to document it: it's only needed in readSBML* files (not by users) and it is
       too trivial to justify documentation maintenance for maintainer use (I also removed it form export). 
       Replaced dQuote in rule writing in saveSBMLR with escaped double quotes in format (else they come out slanted).
       Created Ops.SBML to handle overloaded "==" to work as before (see 1.31.0 above), replacing equateModels().
       A quick-start Sweave vignette was also added. (2/24/2010)

          
          
 

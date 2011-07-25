"readSBML" <-
    function(filename)
{  # takes SBML in filename.xml and maps it to a SBML class model 
# using both Sax and DOM (for mathml) based parsing.
  
  #Function to handle the SBML document - XML element analyzer to create the necessary Objects
  sbmlHandler <- function () 
  {
    sbml<-"x"
    modelid<-"x"					#storing model id
    modelname<-"x"					#storing the model name if separately given

    lnotes<-NULL					# model notes
	
    compartments <- list()				#list of compartments
    reactLaws <- list()					#list of reaction laws
    species <- list()					#list of species
    rules<- list()					#list of rules
    reactions<- list()					#list of reactions
    globalParameters<-list()			#list of Global parameters. Local params are declared locally(within a reaction) 
    ParametersList<- list()				#list of parameter objects	
    reactants=NULL				
    products=NULL
    modifiers=NULL
    currRxnID=NULL
    parameters=NULL   					# local to rate law
    parameterIDs=NULL   				# local to rate law
      
    globalParameterIDs=NULL   
    
    notes=FALSE; reactant=FALSE; product=FALSE
    law=FALSE; parameter=FALSE; math=FALSE
    
     #Sub function to handle first element of the XML document

    .startElement <- function(name, atts, ...) {
     #cat("Start: Name =",name," ",paste(names(atts),atts,sep=" = "),"\n")	
      if(name=="sbml")  
	{
		sbml<<-atts
	}
      if(name=="model")  
	{
		numitems <- length(atts)
		
		
		if(numitems < 1)			#if model does not contain a name/id, we give it an arbitrary one.
		{
			modelid[[1]]<<-"BioModel"	
		}

		else if(numitems == 1)				# if only one attribute supplied
			{
					if(is.character(atts[[1]]))	#if the attribute is a string, it must be model name. Copy to also id.
					{
						modelname<<-atts[[1]]   # store model name
						modelid<<-atts[[1]]     # store as model id (copy Essentially)
					}
			}

			else if(numitems > 1)			#both Id and name of model are supplied, read both
				{
					modelname<<-atts[["name"]]   # first element is model name
					modelid <<-atts[["id"]]    # second element has to be model id
				}
				
	}

      if(name=="compartment")  
	{
		compartments[[atts["id"]]]<<-atts
	}

      if(name=="species")  
	{
			species[[atts["id"]]]<<-atts 
	}

      if(name=="assignmentRule")
	{
		rules[[atts["variable"]]]$idOutput<<-atts[["variable"]] 
	}
      if(name=="reaction")  
	{
		lstnames <- names(atts)
		numitems <- length(lstnames)
		nameslist <- list()
		id <- "x"
		reverse <- FALSE
		name <- "x"
		count <- 1
		while( count <= numitems )
		{
			switch(lstnames[[count]],
		        	"id" = { id = atts[[count]]; nameslist[[length(nameslist)+1]] <- "id"}, 
				"reversible" = { reverse = as.logical(atts[[count]]) ;nameslist[[length(nameslist)+1]] <- "reversible" },
				"name" = { name = as.character(atts[[count]]); nameslist[[length(nameslist)+1]] <- "name"}
   			      )
			count <- count + 1
		}
		reactions[[atts["id"]]]$id<<-id
		reactions[[atts["id"]]]$reversible<<-reverse
		currRxnID<<-atts["id"]
		#reactions[[atts["id"]]]$id<<-atts[["id"]]
	        #reactions[[atts[1]]]$reversible<<-as.logical(atts[[2]])
	        #currRxnID<<-atts[1]
	}
      if(name=="listOfReactants")  
	{
		reactant<<-TRUE
	}
      if(name=="listOfProducts")  
	{
		product<<-TRUE
	}
      if(name=="kineticLaw")  
	{
		law<<-TRUE
	}
		
      if(name=="math")  
	{
		math<<-TRUE
	}
      if((name=="speciesReference")&reactant)			
	{
        	reactants<<-c(reactants,species=atts[["species"]])
	}
      if((name=="speciesReference")&product)
	{
	        products<<-c(products,species=atts[["species"]])
	}
      if(name=="modifierSpeciesReference")
	{
	        modifiers<<-c(modifiers,species=atts[["species"]])
	}
      if((name=="parameter")&law)			#parameter encountered within a kinetic law definition
	{
	        #parameterIDs<<-c(parameterIDs,atts[[1]])
	        #parameters<<-c(parameters,atts[[2]])
		parameterIDs<<-c(parameterIDs,atts[["id"]])
	        parameters<<-c(parameters,atts[["value"]])
	}
      if((name=="parameter")&!law)			#parameter encountered outside a kinetic law definition - So in globalparamslist
	{
		#old code
        	#globalParameterIDs<<-c(globalParameterIDs,atts[[1]])
        	#globalParameters<<-c(globalParameters,as.numeric(atts[[2]]))
		#cat("within parameters:", atts[["id"]], atts[["value"]], "\n")

		#parameterIDs<<-c(parameterIDs,atts[["id"]])
	        #parameters<<-c(parameters,atts[["value"]])

		globalParameterIDs<<-c(globalParameterIDs,atts[["id"]])
        	globalParameters<<-c(globalParameters,as.numeric(atts[["value"]]))

		ParametersList[[atts["id"]]] <<- atts		#our new list of Parameter Objects
	}
}# .startelement function ends  
    
    .endElement <- function(name) 			#function for the closing tags
   {
      if(name=="listOfReactants")  
	{
		reactant<<-FALSE  
	}
      if(name=="listOfProducts")  
	{
		product<<-FALSE   
	}
      if(name=="kineticLaw")  
	{
		law<<-FALSE
	}
      if(name=="math")  
	{		
		math<<-FALSE
	}
      if((name=="listOfParameters")&(!law)) 
	{
		names(globalParameters)<<-globalParameterIDs 
	} 
      if(name=="reaction")  
	{
	        names(reactants)<<-NULL
	        names(modifiers)<<-NULL
	        names(products)<<-NULL
	        reactions[[currRxnID]]$reactants<<-reactants
	        reactions[[currRxnID]]$modifiers<<-modifiers
	        reactions[[currRxnID]]$products<<-products
	        parameters<<-as.numeric(parameters)
	        names(parameters)<<-parameterIDs
	        reactions[[currRxnID]]$parameters<<-parameters
	        reactants<<-NULL;products<<-NULL
	        modifiers<<-NULL;
		parameters<<-NULL;parameterIDs<<-NULL     
	}
    }
    
    .text <- function(x, ...) 
	{
      		if (!math) lnotes<<-c(lnotes,x)
		#  cat("Txt:", x,"\n")
	}
    
    getModel<- function() 
    { 
      fixComps=function(x) 
	{
		lstnames <- names(x)
		count <- 1
		numit <- length(lstnames)
		id <- "x"
		size <- 0
		name <- "x"
		nameslist <- list()
		while( count <= numit )
		{
			switch(lstnames[[count]],
		        	"id" = { id = x[[count]]; nameslist[[length(nameslist)+1]] <- "id"}, 
				"size" = { size = as.numeric(x[[count]]) ;nameslist[[length(nameslist)+1]] <- "size" },
				"name" = { name = as.character(x[[count]]); nameslist[[length(nameslist)+1]] <- "name"}
			      )
			count = count + 1
		}

		#---DEBUG----
		#cat("Compartment id:", id,"\n", "Compartment size:" , size, "\n","Compartment name:", name, "\n")
		#-------------
		if(numit == 2)						# only 2 attributes present. We need to find them.
		{
			

			if(id == "x")		#id not set but name and size are.
			{
				id <- "default"		
			}
			else if(name == "x")    #name not set, we copy the id.
			{
				name <- id
			}
			else if(size== "0")	#size not set
			{
				size <- 1 	#arbitrary setting as 1
			}
			lst = list(id,size,name)	
			names(lst)<-c("id","size","name")
			lst
		}
		else if(numit == 3)					# 3 attributes present.
		{
			lst = list(id,size,name)
			names(lst)<-c("id","size","name")
			lst 
		}
		#old-------
		#lst=list(x[[1]],as.numeric(x[[2]]) ); 	#compartment id and compartment size
		#names(lst) <- names(x) 	
		#---------
	}
      
      fixSpecies=function(x) 
	{
		#cat (names(x), "\n")
		#cat(toString(x) , "\n")
		numitems <- length(x)
		lstnames <- names(x)
		count <-1
		id <- "x"			#species Id
		ic <- 0				#species initial concentration
		compart <- "def"		#species compartment
		bc <- FALSE			#species boundary condition
		name <- "def"
		nameslist <- list()
		while( count <= numitems)
		{
			switch(lstnames[[count]],
		        	"id" = { id <- x[[count]]; nameslist[[length(nameslist)+1]] <- "id"},
				"name" = { name <- x[[count]]; nameslist[[length(nameslist)+1]] <- "name"},
				"initialConcentration" = { ic <- as.numeric(x[[count]]) ;nameslist[[length(nameslist)+1]] <- "ic" },
				"compartment" = { compart <- as.character(x[[count]]); nameslist[[length(nameslist)+1]] <- "compartment"},



				"boundaryCondition" = { bc <- as.logical(x[[count]]); nameslist[[length(nameslist)+1]] <- "bc"}
			      )
			count = count + 1
		}
		#old lst=list(x[[1]],as.numeric(x[[2]]),x[[3]],as.logical(x[[4]])); 
		#lst = list(id,ic,compart,bc, name)
		lst = list(id,as.numeric(ic), compart, as.logical(bc))
		names(lst) <- c("id","ic","compartment","bc")
	        #names(lst)<-c("id","ic","compartment","bc", "name"); 
	        lst 
	}


	   fixParams=function(x) 
	{
		numitems <- length(x)
		lstnames <- names(x)
		count <-1
		id <- "x"			#Parameter Id
		value <- 0			#Parameter value
		name <- "def"
		constant <- FALSE
		nameslist <- list()
		while( count <= numitems)
		{
			switch(lstnames[[count]],
		        	"id" = { id <- x[[count]]; nameslist[[length(nameslist)+1]] <- "id"},
				"name" = { name <- x[[count]]; nameslist[[length(nameslist)+1]] <- "name"},
				"value" = { value <- as.numeric(x[[count]]) ;nameslist[[length(nameslist)+1]] <- "value" },
			"constant" = { constant <- as.logical(x[[count]]) ; nameslist[[length(nameslist)+1]] <- "constant"}
			     )
			count = count + 1
		}
		
		lst = list(id,as.numeric(value))
		names(lst) <- c("id","value")
	        lst 
	}		
#
      compartments=sapply(compartments,fixComps, simplify = FALSE)
      species=sapply(species,fixSpecies, simplify = FALSE)     			#this keeps the better looks in the SBMLR model definition file

      ParametersList = sapply(ParametersList, fixParams, simplify = FALSE)	#building params list
      
      list(sbml=sbml,id=modelid[[1]], notes=lnotes,compartments=compartments,
          species=species,globalParameters=globalParameters, ParametersList=ParametersList, rules=rules,reactions=reactions)
    }#end of getmodel function
    
    list(.startElement = .startElement, .endElement = .endElement, .text = .text,   # , dom = function() {con}
        getModel = getModel     
    )
  }
#  END handler definition
  
# NOTE: though handlers are neat, one must question if the added baggage is worth it, i.e. compare to read.SBML of older versions 
  
# *********************************************************************************
# The next block of three functions converts mathML XMLnode objects into R expression objects
# This approach is better than the old read.SBML approach in that the parenthesis overkill is avoided!

  mathml2R <-function(node)  
	{
		UseMethod("mathml2R", node)
	}
  
  mathml2R.XMLDocument <-function(doc) 
	{
		return(mathml2R(doc$doc$children))
	}
  
  mathml2R.default<-function(children) 
  {  expr <- expression()  # this gets used when a "list" of children nodes are sent in
    for(i in children) {    expr <- c(expr, mathml2R(i))  }
	#cat("Expression main: " ,toString(expr), "\n")
    return(expr)
  }
  
  mathml2R.XMLNode <-function(node)
   {
    		nm <- xmlName(node) 
		#DEBUG------------------- 
		#cat("xml node:", toString(nm), "\n")
		#------------------------

	   # cat("XMLNode: node name is ",nm," and the node class is",class(node),"\n")	#was commented

	    if(nm=="power"||nm == "divide"||nm =="times"||nm=="plus"||nm=="minus")
		 {
	      		op <- switch(nm, power="^", divide="/",times="*",plus="+",minus="-")
			      val <- as.name(op)	
		 }
	 else if((nm == "ci")|(nm == "cn")) 
	{
		#cat("within mathml function ", node$children[[1]]$value, "\n")
      		if(nm == "ci") val <- as.name(node$children[[1]]$value)
		if(nm == "cn") val <- as.numeric(node$children[[1]]$value)
	} 
	else if(nm == "apply") 
	{
      		val <- mathml2R(node$children)
		mode(val) <- "call"
 	}
	else  
	{
		cat("error: nm =",nm," not in set!\n")
	}
	#cat("Expression value: " , toString(val), "\n")
    return(as.expression(val))
  }
# ********** END the mathML2R block of method based on node type codes  *************************
  
  
  
# The next two functions are used by rules and were taken straight from read.SBML
# The idea is that SBML doesn't provide a list of atoms/leaves with rules, so we have to create them
# to place them in their model slots, and to use them to create the R function definition for the rule
# using makeLaw with a null for parameters, since they are passed global for rules.
  ML2R<- function(type)   # map MathML operator symbols into R symbols
    switch(type,
        "times" = "*",
        "divide" = "/",
        "plus" = "+",
        "minus" = "-",
        "power" = "^",
        "exp" = "exp",
        "ln" = "log",
        "not found") # end definition of ML2R
  
  #function used to retrieve the xml law elements and create the strlaw equation.
  getRuleLeaves<-function(math) 
  { n=length(math)
    S=c(NULL)
    op=ML2R(xmlName(math[[1]]))

    for (j in 2:n )		#reading from 2 to avoid the operator
      if ((xmlName(math[[j]])=="ci")|(xmlName(math[[j]])=="cn"))  S=c(S,as.character(xmlValue(math[[j]]))) else 
        S=c(S,Recall(math[[j]])  ) 
    S			#return value of the function
  } 

 #function used to retrieve the xml law elements and create the strlaw equation for rules.
  getRuleLeavesnew<-function(math) 
  { n=length(math)

    S=c(NULL)
    op=ML2R(xmlName(math[[1]]))

    for (j in 2:n )
      if ((xmlName(math[[j]])=="ci")|(xmlName(math[[j]])=="cn"))  S=c(S,as.character(xmlValue(math[[j]]))) else 
        S=c(S,Recall(math[[j]])  ) 
    S			#return value of the function
  } 
  
  
  if(!require(XML)) print("Error in Read.SBML(): First Install the XML package http://www.omegahat.org/RSXML")
  
  edoc <- xmlEventParse(filename,handlers=sbmlHandler(),ignoreBlanks = TRUE)	#object which is created by readsbml first.
										#from xmlparser package. arguments(file_to_parse, handler, handling 											#spaces boolean)
  model=edoc$getModel()								#model object created by taking previous and calling getModel()

  doc <- xmlTreeParse(filename,ignoreBlanks = TRUE) 				#creating XML dom tree

  model$htmlNotes=doc$doc$children$sbml[["model"]][["notes"]] 			# Traverse DOM and retrieve the notes child value
  rules=doc$doc$children$sbml[["model"]][["listOfRules"]]			# Traverse the DOM and retrieve the list of rules
  reactions=doc$doc$children$sbml[["model"]][["listOfReactions"]]		# Traverse the DOm and retrieve the list of reactions
  
  globalParameters=names(model$globalParameters)				# Set name values for the global parameters list
  
  
  nRules=length(rules)								# Number of rules in the sbml file(assuming all assignment)

  if (nRules>0)
	{	for( i in 1:(nRules-1)) {  
			mathml<-rules[[i]][["math"]][[1]]
		
		      #cat( "mathml:", toString(rules[[i]][["math"]][[1]]), "\n") gives you the entire math block
		      model$rules[[i]]$mathmlLaw=mathml
			
		      e<-mathml2R(mathml)					   #creates the expression law given all the child nodes under <math>.
			
		      model$rules[[i]]$exprLaw<-e[[1]]
		      model$rules[[i]]$strLaw<-gsub(" ","",toString(e[1]))
		      leaves<-getRuleLeaves(mathml)
		      r<-model$rules[[i]]$inputs<-setdiff(leaves,globalParameters) # must deduce inputs by substracting global params
			
		      model$rules[[i]]$law=makeLaw(r,NULL,model$rules[[i]]$exprLaw)
		}
  	} 


  nReactions=length(reactions)
  if (nReactions>0)
	{
	#    rIDs=NULL;  
	    for (i in 1:nReactions)
	    {
	      model$reactions[[i]]$mathmlLaw=reactions[[i]][["kineticLaw"]][["math"]][[1]]
	      e=mathml2R(reactions[[i]][["kineticLaw"]][["math"]][[1]])
	      model$reactions[[i]]$exprLaw=e[[1]]
	      model$reactions[[i]]$strLaw=gsub(" ","",toString(e[1]))
	      #paste(as.character(model$reactions[[i]]$expr)[c(2,1,3)],collapse=""))
      	      r=c(model$reactions[[i]]$reactants, model$reactions[[i]]$products)	#build using both reactants and products objects. TODO - add compartments
	      p=names(model$reactions[[i]]$parameters)
	      m=model$reactions[[i]]$modifiers
             
	      r=c(r,m)
	      e=model$reactions[[i]]$exprLaw
	      model$reactions[[i]]$law=makeLaw(r,p,e)
		#      rIDs[i]<-model$reactions[[i]]$id
	    }
		# This is for indexing by names/IDs of reactions
		#    names(model$reactions)<-rIDs
		#    names(model$reactions)<-sapply(model$reactions,function(x) x$id)
	}
  
	#---DEBUG CODE
	cat("Number of species: " , length(model$species), "\n")
	cat("Number of rules: ", nRules, "\n")
	cat("Number of Global Parameters: " , length(globalParameters), "\n")
	cat("Number of reactions: " , nReactions, "\n")
	cat("Parsing Successful !" , "\n")
	#-----------


  class(model)<-"SBML"
  model
}

# the following is called by both readSBML and readSBMLR so it outside where both can reach it.
# Note that keeing it here instead of in a separate file => no need to document it
"makeLaw"<-function(r,p,e){
# takes reactant list r, parameter list p and rate law R expression e 
# and makes a reaction rate law function out of them.
  lawTempl=function(r,p=NULL){ }
  i=2
  if(!is.null(p))
    for (j in 1:length(p)){
      body(lawTempl)[[i]]<-call("=",as.name(p[j]),call("[",as.name("p"),p[j]))
      i=i+1}
  for (j in 1:length(r)){ 
    body(lawTempl)[[i]]<-call("=",as.name(r[j]),call("[",as.name("r"),r[j]))
    i=i+1}
  body(lawTempl)[[i]]<-e
  lawTempl
}


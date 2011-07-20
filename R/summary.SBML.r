"summary.SBML"<-function(object, ...)
{
  model=object
  sIDs=names(model$species)
  pIDs=names(model$ParametersList)
  rIDs=names(model$reactions)
  ruleIDs=names(model$rules)
  nReactions=length(model$reactions);nSpecies=length(model$species);nRules=length(model$rules) 
  nParams = length(model$ParametersList)
  
#Parameters
  P0 = NULL; VP = NULL
 if(nParams > 0)
	{
		  for (i in 1:nParams){
		    P0[i]=model$ParametersList[[i]]$id; 
		    VP[i]=model$ParametersList[[i]]$value
		  }
	}
	names(VP)<-pIDs
# Species
  S0=NULL;BC=NULL # initialize 
  for (i in 1:nSpecies){
    BC[i]=model$species[[i]]$bc; 
    S0[i]=model$species[[i]]$ic
  }
  names(S0)<-sIDs 
  names(BC)<-sIDs 
  y0=S0[BC==FALSE]
  nStates=length(y0)
  
  attach(model$globalParameters)  # e.g. for Vmax global coordination in RNR

# Reactions
  rLaws=NULL;V0=NULL # initialize
  for (j in 1:nReactions) {
    rLaws[j]<-model$reactions[[j]]$strLaw		#this gives you null which is wrong
	
    V0[j]=model$reactions[[j]]$law(S0[c(model$reactions[[j]]$reactants,model$reactions[[j]]$modifiers,model$reactions[[j]]$products)], model$reactions[[j]]$parameters)
  }
  names(rLaws)<-rIDs
  names(V0)<-rIDs
  detach(model$globalParameters)  
  
 
# Incidence Matrix
  incid=matrix(rep(0,nStates*nReactions),nrow=nStates)
  indx=(1:nSpecies)[BC==FALSE]
  for (i in 1:nStates)
    for (j in 1:nReactions)
    {if ( is.element(model$species[[indx[i]]]$id, model$reactions[[j]]$products)) 
        incid[i,j] = summary(factor(model$reactions[[j]]$products))[[model$species[[indx[i]]]$id]]
      if ( is.element(model$species[[indx[i]]]$id, model$reactions[[j]]$reactants)) 
        incid[i,j] = incid[i,j]-summary(factor(model$reactions[[j]]$reactants))[[model$species[[indx[i]]]$id]]  }     
  rownames(incid)<-names(y0)
  
# return a list of model information
  
  options(stringsAsFactors=FALSE)
  DFs=data.frame(index=1:nSpecies,initialConcentrations=S0,boundaryConditions=BC); row.names(DFs)<-sIDs

  DFr=data.frame(index=1:nReactions,Laws=rLaws,initialFluxes=V0);   row.names(DFr)<-rIDs
  list(nSpecies=nSpecies,sIDs=sIDs,S0=S0,BC=BC,nStates=nStates,y0=y0,nReactions=nReactions,rIDs=rIDs,rLaws=rLaws, V0=V0, 
     incid=incid,nRules=nRules,ruleIDs=ruleIDs,species=DFs,reactions=DFr, VP = VP, P0=P0)

	 # list(nSpecies=nSpecies,sIDs=sIDs,S0=S0,BC=BC,nStates=nStates,y0=y0,nReactions=nReactions,rIDs=rIDs,rLaws=rLaws, V0=V0, 
      #incid=incid,nRules=nRules,ruleIDs=ruleIDs,species=DFs)
	#}
	# initial out
	# DFs=data.frame(index=1:nSpecies,initialConcentrations=S0,boundaryConditions=BC); row.names(DFs)<-sIDs
	#DFp=data.frame(index=1:nParams,id=P0,value=VP); row.names(DFp)<-pIDs
	#list(nSpecies=nSpecies,sIDs=sIDs,S0=S0,BC=BC,species=DFs,nParams=nParams,pIDs=pIDs,P0=P0,VP=VP, Parameters=DFp)
}




"Ops.SBML"<-function(e1,e2) {
# The polynom package showed me the way. It seems e1 and e2 must be the names of the arguments.
# Apparently .Generic is passed globally from the command line.  
  m1=summary(e1)
  m2=summary(e2)
  return(switch(.Generic,
          "==" = list(species=m1$species==m2$species, reactions=m1$reactions==m2$reactions),
          stop("unsupported operation")))
}




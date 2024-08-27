require(expm)

StatesSCOpti<-function(k, n) {
  ## This function returns the states of the structured coalescent of Herbots in a matrix E
  ## k is the number of lineages in the initial sample
  ## n is the number of islands
  E<-c()
  ncolE<-0
  RecursionMig<-function(p, l) {
    ## function creating all the states of the chain l islands
    ## having exactly p lineages

    M<-c() ## Initialisation of return
    if (l==1) return(matrix(p,nrow=1,ncol=1))
    if (p==0) return(matrix(0,nrow=l,ncol=1))
    if (p==1) return(diag(1,l))
    ncolM<-0
    for (i in p:0) {
      temp1<-RecursionMig(p-i, l-1)
      temp2<-matrix(c(rep(i,ncol(temp1)),t(temp1)),nrow=nrow(temp1)+1,ncol=ncol(temp1),byrow=TRUE)
      M<-matrix(c(M,temp2),nrow=nrow(temp2),ncol=ncolM+ncol(temp2))
      ncolM<-ncol(M)
      }
    return(M)
  }
  for (i in k:1) {
    temp<-RecursionMig(i, n)
    E<-matrix(c(E,temp),nrow=n,ncol=ncolE+ncol(temp))
    ncolE<-ncol(E)
  }
  return(E) }

IsCoalescence<-function(i, j, E) {
  ## Function returning a vector of size 2 :
  ## first element is 1 if the coalescence of i and j is possible, 0 otherwise
  ## second element is the number of the island in which i and j will coalesce (0 if no coalescence)
  if (sum(E[,i]!= E[,j])==1&& E[which(E[,i]!=E[,j]),i]>=2&&E[which(E[,i]!=E[,j]),i]==(E[which(E[,i]!= E[,j]),j]+1)) {
    return(c(TRUE,which(E[,i]!=E[, j])))
    }
  else {
    return(c(FALSE,0))
  }
}

IsMigration<-function(i, j, E) {
  ## Function returning a vector of size 3 :
  ## First element is 1 if the transition from i to j is a migration, 0 otherwise
  ## Two next elements are the number of the island of departure and the number of
  ## the island of arrival of the migration (0 if no migration)

  if ((sum(E[,i]!=E[,j])==2) && sum(E[,i])!=1 && sum(E[,j])!=1){
    indices<-which(E[,i]!=E[,j])
    if (E[indices[1],i]==(E[indices[1] ,j]-1)&&E[indices[2],i]==(E[indices[2], j]+1)) {
      return(c(1,indices[2],indices[1])) }
    else{
      if (E[indices[1],i]==(E[indices[1] ,j]+1)&&E[indices[2],i]==(E[indices[2], j]-1)) {
        return(c(1,indices[1],indices[2])) }
      else { return(c(0,0,0)) }
    } }
  else { return( c(0,0,0) ) } }

QmatrixSC<-function(M,c,k=2) {
  ## This function returns the Q-matrix corresponding to the structured coalescent
  ## of Herbots. M contains the migration rates between the islands, c contains
  ## the sizes of the islands. k is the number of lineages in the initial sample
  ## The function is usable even if the data do not respect the condition
  ## of constant size of the islands (relation between M and c)
  states<-StatesSCOpti(k, ncol(M))
  nbStates<-ncol(states)
  Q<-matrix(0, ncol=nbStates, nrow=nbStates)
  for (alpha in 1:nbStates) {
  for (beta in subset(1:nbStates, 1:nbStates!=alpha)) {
    testMig<-IsMigration(alpha, beta, states)
    testCoa<-IsCoalescence(alpha, beta, states)
    if (testMig[1] ) {
      Q[alpha,beta]<-states[testMig[2], alpha]*M[testMig[2], testMig[3]]/2
    }
    if (testCoa[1] ) {
      Q[alpha,beta]<-(states[testCoa[2], alpha]-1)*states[testCoa[2], alpha]/(2*c[testCoa[2]])
    }
  }
  Q[alpha,alpha]<--sum(Q[alpha,] ) }
return(Q) }

TransitionSC<-function(Q) {
  ## This function returns the transition semi-group from the infinitesimal generator Q
  P<-function(t) {return(expm(t*Q))}
  return(P) }

MaskCoaSC<-function(E, k=1) {
  ## This function returns a matrix of size sum(E[,1]) (corresponding to the
  ## number of lineages initially present). This matrix is used as a logical mask
  ## to select which states correspond to a k-th coalescence (by default k=1)
  nbLineages<-colSums(E)
  mask<-matrix(0, nrow=length(nbLineages), ncol=length(nbLineages))
  for (i in 1:length(nbLineages)) {
    for (j in 1:length(nbLineages)) {
      if (nbLineages[i]>=nbLineages[j]+k ) {mask[i, j]<-1 } } }
  return(mask) }

DistribSC<-function(P, E, states, weights=rep(1, length(states)), k=1) {
  ## This function returns the distribution function of the k-th coalescence
  ## conditionally to the states from which we can start.
  ## P is the transition semi-group of the chain
  ## E is the matrix containing the states of the chain
  ## states is a vector containing the initial states
  ## weights contains the corresponding probabilities (uniform by default)
  ## k is the number of coalescences
  F<-function(t) {
    temp<-P(t) *MaskCoaSC(E, k)
    return(sum(temp[states,]*(weights/sum(weights)))) }
  return(F) }

DensitySC<-function(P, Q, E, states, weights=rep(1, length(states)), k=1) {
  ## This function returns the density function of the k-th coalescence
  ## conditionally to the states from which we can start.
  ## P is the transition semi-group of the chain
  ## Q is the infinitesimal generator of the chain
  ## E is the matrix containing the states of the chain
  ## states is a vector containing the initial states
  ## weights contains the corresponding probabilities (uniform by default)
  ## k is the number of coalescences
  f<-function(t) {
    temp<-(Q%*%P(t))*MaskCoaSC(E, k)
    return(sum(temp[states,]*(weights/sum(weights))))}
  return(f) }

RateSC<-function(P, Q, E, states, weights=rep(1, length(states)), k=1) {
  ## This function returns the rate of the k-th coalescence
  ## conditionally to the states from which we can start.
  ## P is the transition semi-group of the chain
  ## Q is the infinitesimal generator of the chain
  ## E is the matrix containing the states of the chain
  ## states is a vector containing the initial states
  ## weights contains the corresponding probabilities (uniform by default)
  ## k is the number of coalescences
  lambda<-function(t) {
    return((1-DistribSC(P, E, states, weights, k)(t))/DensitySC(P, Q, E, states, weights, k)(t)) }
  return(lambda) }

IICRvector <- function(M,E,e,c,range,k=2){
  Q<-QmatrixSC(M,c,k)
  P<-TransitionSC(Q)
  lambda<-RateSC(P, Q, E, e)
  vec1 <- c()
  for (p in range) {
    vec1<-c(vec1,lambda(p))
  }
  return(vec1)
}

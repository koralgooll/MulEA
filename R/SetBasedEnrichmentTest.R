#####################################################################################################
# Function for FDR corrected hypergeometric enrichment test
#####################################################################################################

# snow and rlecuyer packages should be installed


## Arguments of HyperGeomFDR function:
#	steps:		the rounds of simulations (a single number)
#	pool:		background genes (character vector)
#	select:		genes to investigate (character vector)
#	DB: 		the genes set used for enrichment analysis (character list)
#	nthreads:	number of threads to use (a single number)

## Description of the hypergeometric test
# phyper(q, m, n, k, lower.tail = TRUE, log.p = FALSE)
# Arguments:
#        q: vector of quantiles representing the number of white balls
#           drawn without replacement from an urn which contains both
#           black and white balls.
#        m: the number of white balls in the urn.
#        n: the number of black balls in the urn.
#        k: the number of balls drawn from the urn.
#
# x=length(intersect(select,DB_i))    	#Number of common genes between DB and select
# m=length(intersect(pool,DB_i))        #Number of common genes between DB and pool
# n=length(pool)-length(intersect(pool,DB_i))     #Number of non-pool genes among DB (setdiff)
# k=length(select)                    	#Number of genes in select
# P_val=dhyper(length(intersect(select,DB_i)), length(intersect(pool,DB_i)), length(pool)-length(intersect(pool,DB_i)), length(select))
#
# wikipedia
# N = length(pool)
# K = length(intersect(pool,DB_i))
# n = length(select)
# k = length(intersect(select,DB_i))
#
# (choose(length(intersect(pool,DB_i)),length(intersect(select,DB_i)))*
#      choose(length(pool)-length(intersect(pool,DB_i)),length(select)-length(intersect(select,DB_i))))/
#      choose(length(pool),length(select))

# m= intersect(select,DB_i)
# n= intersect(bg,DB_i)-m
# k=length (DB_i)


# P_val=(choose(length(intersect(BG,DB_i)),length(intersect(IG,DB_i)))*choose(length(BG)-length(intersect(BG,DB_i)),length(IG)-length(intersect(IG,DB_i))))/choose(length(BG),length(IG))

#' PRIVATE class : An S4 class to represent a Hypergeometric tests in Mulea.
#'
#' @slot gmt A data.frame representing GMT's reprezentation of model.
#' @slot testData A data from expeciment to analize accross model.
#' @slot pool A background data to count test.
#' @return SetBasedEnrichmentTest object. Used as private function.
#' @examples
#' \dontrun{
#' #It is a private s4 object. Look at SetBasedTest's examples.
#' }
SetBasedEnrichmentTest <- setClass("SetBasedEnrichmentTest",
                       slots = list(
                         gmt = "data.frame",
                         testData = "character",
                         pool = "character",
                         numberOfPermutations = "numeric",
                         test = "function"
                       ))

setMethod("initialize", "SetBasedEnrichmentTest",
          function(.Object,
                   gmt = data.frame(),
                   testData = character(),
                   pool = character(),
                   numberOfPermutations = 1000,
                   test = NULL,
                   ...) {

            .Object@gmt <- gmt
            .Object@testData <- testData
            .Object@pool <- pool
            .Object@numberOfPermutations = numberOfPermutations

            .Object@test <- function(testObject) {

                pool <- NULL
                if (0 == length(testObject@pool)) {
                    pool <- unique(unlist(.Object@gmt[, 'listOfValues']))
                } else {
                    pool <- unique(testObject@pool)
                }

                DB <- .Object@gmt[, 'listOfValues']
                names(DB) <- .Object@gmt$ontologyId

                testResults <- set.based.enrichment.test(steps = .Object@numberOfPermutations, pool = pool,
                                                        select = .Object@testData, DB = DB)
                testResults
            }

            .Object
          })

#' @describeIn MuleaHypergeometricTest runs test calculations.
#' @param testObject Object of s4 class represents Mulea Test.
#' @return runTest method for MuleaHypergeometricTest object. Used as private function.
#' @examples
#' \dontrun{
#' #It is a private method. Look at runTest of SetBasedTest's examples.
#' }
setMethod("runTest",
          signature(testObject = "SetBasedEnrichmentTest"),
          function(testObject) {
            testObject@test(testObject)
          })




set.based.enrichment.test=function(steps, pool, select, DB, nthreads=4) {


  # csinalok egy masolatot a konverzio elott, mert ervin kodjanak az kell
  DB0<-DB
  select0<-select
  pool0<-pool

  ## convert the database and the select and pollt to sorted integer lists
  gene.list<-unique(unlist(DB))
  gene.string.to.integer<-data.frame(string.id=gene.list ,integer.id=seq(gene.list))

  DB<-lapply(DB, FUN = function(list1) which(gene.string.to.integer$string.id %in% list1))
  pool <- which(gene.string.to.integer$string.id %in% pool)
  select <- which(gene.string.to.integer$string.id %in% select)

  ############

  DB_names=names(DB)
  num_DB=length(DB)
  size_pool=length(pool)
  size_select=length(select)


  DB_in_select=integer(num_DB)
  DB_in_pool=integer(num_DB)
  Genes_in_DB=integer(num_DB)
  P_val=double(num_DB)
  R_obs=integer(num_DB)

  # for every DB entity in the DB list
  for (i in 1:num_DB) {
    # create a vector of genes connected to the i-th DB category
    DB_i=DB[[i]]
    # hypergometric test
    DB_in_select[i]=int.list.intersect(select,DB_i)	#q: number of common genes between a DBterm and select
    DB_in_pool[i]=int.list.intersect(pool,DB_i)	#m: number of common genes between DBterm and BackGround
    Genes_in_DB[i]=length(DB_i)

    #n:  number of non-pool genes among DB
    #k: number of genes in select
    # phyper(q, m, n, k, lower.tail = TRUE, log.p = FALSE)
    P_val[i]=1-phyper(DB_in_select[i]-1, DB_in_pool[i], size_pool-DB_in_pool[i], size_select)  ## TODO itt a size_select kell?
  }

  P_val_round=round(P_val, digits=15) ## can change the digits, this is important for the precision of '0' is R
  for (i in 1:num_DB) { # TODO ez egy rang szamitas.
    R_obs[i]=sum(P_val_round<=P_val_round[i])
  }
  P_val_df=data.frame(DB_names, DB_in_select, DB_in_pool, Genes_in_DB, P=P_val, P_adj_Bonf=p.adjust(P_val, method="bonferroni"), P_adj_BH=p.adjust(P_val, method="BH"), R_obs)

  ######
  # simualtion
  ######
  # random sampling from pool (background genes)
  # The time consuming step. The simulation here can be parallelized


  #   	require(snow)
  #	require(rlecuyer)
  #
  #	seeds=sample(seq(1e4,1e6),6) # max number of seeds for RNGstream is 6 TODO he???
  #	cl=makeCluster(nthreads, type="SOCK")
  #	clusterExport(cl,"int.list.intersect")
  #	clusterSetupRNG(cl, type='RNGstream', seed=seeds)
  #	simulation.result=clusterApply(cl, rep(ceiling(steps/nthreads), nthreads), sim_hyperGeom_FG, pool, length(select),DB, P_val_df$DB_in_pool ) # return a list
  #	stopCluster(cl)
  #
  #	print("DEBUG: tobb szalrol jovo tablazatok aggregalasa")
  #	simulation.result<-do.call("rbind",simulation.result)
  #	simulation.result<-aggregate(simulation.result$x,simulation.result[,c("intersect.size", "DB_in_pool" )],sum )
  #

  #	set.seed(1234567)
  #	simulation.result=sim_hyperGeom_FG(steps, pool, length(select),DB, P_val_df$DB_in_pool ) # return a list
  #
  #	print("DEBUG: phyper() szamitas")
  #	simulation.result$p <- 1-phyper(simulation.result$intersect.size-1, simulation.result$DB_in_pool, length(pool)-simulation.result$DB_in_pool,  length(select))
  #	simulation.result<-simulation.result[order(simulation.result$p),]
  #	simulation.result<-rbind(data.frame(intersect.size=-1 ,DB_in_pool=-1, x=0,p=0),simulation.result) # a tabla elejere rakok egy dumy sort x=0,p=0 -val; ez segit a binaris keresesnek
  #	simulation.result$cum.sum.x<-cumsum(simulation.result$x)
  #
  #	simulation.result$p_round<-round(simulation.result$p, digits=15)
  #
  #	simulation.result1<-simulation.result

  if(length(select)==0 )
  {
    P_val_df$R_exp=NaN
    P_val_df$FDR=NaN

    return(P_val_df)
  }

  flush.console(); # flush-olok, hoga a debug uzenetek ne keveredjenek ossze
  ## Most nézzük Ervin kódját
  names(DB0) <- NULL
  #steps = 10000
  simulation.result =  tryCatch(
    trial(
      DB0,
      unique(c(unlist(DB0), pool)),
      pool0,
      length(select0),
      steps,
      1234567)
    , error = print)  # RCPP hívás
  names(simulation.result)<-c("DB_in_pool","intersect.size","x")
  simulation.result$p <- 1-phyper(simulation.result$intersect.size-1, simulation.result$DB_in_pool, length(pool)-simulation.result$DB_in_pool,  length(select))
  if(! all(is.finite(simulation.result$p)) ) {stop("ERROR_002")} # TODO handle somehow the exception

  simulation.result<-simulation.result[order(simulation.result$p),]
  simulation.result<-rbind(data.frame(intersect.size=-1 ,DB_in_pool=-1, x=0,p=0),simulation.result) # a tabla elejere rakok egy dumy sort x=0,p=0 -val; ez segit a binaris keresesnek
  simulation.result$cum.sum.x<-cumsum(simulation.result$x)
  simulation.result$p_round<-round(simulation.result$p, digits=15)



  R_exp=integer(num_DB)

  #P_Sim_vec=as.vector(unlist(P_Sim_vec))
  #P_Sim_round=round(P_Sim_vec, digits=15)
  #P_Sim_round=sort(P_Sim_round)
  cnt.of.ones<-sum(simulation.result$p_round==1)
  NN<-simulation.result$cum.sum.x[nrow(simulation.result)] # ez a cumsum.x-bol az utolso/legnagyobb
  for (i in 1:num_DB) { # TODO ezt is betenni a parhuzamos szamitasba
    if(P_val_round[i]>=1)
    {
      R_exp[i]= NN  # ez a resz gyorsit sokat
    }else{
      target<-P_val_round[i]

      # binary search: az a1 az alsó határ, az a2 a felső, és addig kozelitem oket egymashoz amig osszeernek
      a1 <- 1
      a2 <- nrow(simulation.result) #-cnt.of.ones
      #P_Sim_round[a2+1] # ez mar 1
      while(a1+1 < a2)
      {

        a3<-floor((a1+a2)/2)
        current <- simulation.result$p_round[a3]
        if( current <= target)
        {
          a1<-a3
        }
        else
        {
          a2<-a3
        }
      }
      #		 		R_exp[l]<-a2
      R_exp[i]<-simulation.result$cum.sum.x[[a1]]
    }
  }
  P_val_df$R_exp=R_exp/steps
  P_val_df$FDR=P_val_df$R_exp/R_obs


  return(P_val_df)
}

######################################

sim_hyperGeom_FG=function(steps, pool, size_select, DB,DB_in_pool) {
  num_DB=length(DB)
  P_Sim_mat=matrix(numeric(num_DB*steps), ncol=steps)
  size_pool=length(pool)

  matrix_intersect_size=matrix(numeric(num_DB*steps), ncol=steps)


  #tmp<-table(matrix_intersect_size,matrix_2)


  for (j in 1:steps) {
    Rand.select=sort(sample(pool, size_select))
    for (i in 1:num_DB) {
      matrix_intersect_size[i,j]=int.list.intersect(Rand.select, DB[[i]])
    }
  }


  parameter.df<-data.frame()
  for(v in sort(unique(DB_in_pool)))
  {
    idx <- DB_in_pool==v
    tmp <- matrix_intersect_size[idx,]
    tmp <- as.integer(tmp)

    t1<-aggregate(tmp,by=list(intersect.size=tmp),FUN=length )
    t1$DB_in_pool=v

    parameter.df<-rbind(parameter.df,t1)
  }


  return(parameter.df)
}

######################################

# x and y are two sorted int[]
# The function count the number of common numbers.
int.list.intersect <- function(x,y)
{
  n <- length(x)
  m <- length(y)

  i=1
  j=1

  cnt<-0
  while(TRUE)
  {
    if( i>n | j>m ) return(cnt);
    if(x[[i]] == y[[j]])
    {
      cnt <- cnt+1
      i<- i+1
      j<- j+1
    }else if (x[[i]] < y[[j]])
    {
      i<- i+1
    }else
    {
      j<- j+1
    }

  }
  return(cnt)

}

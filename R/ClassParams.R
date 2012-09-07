setClass(Class="Params",
         representation(
                        ontology  ="character",
                        organism  ="character",
                        method    ="character",
                        combine   ="character",
                        dropCodes ="character"),
         prototype=prototype (dropCodes="NULL")
         )

setValidity("Params",
            function(object) {
		if(!exists("GOSemSimEnv")) .initial()
		organs   <- get("SupportedSpecies",envir=GOSemSimEnv)
		onts     <- c("MF", "BP", "CC")
		mets     <- c("Resnik", "Jiang", "Lin", "Rel", "Wang")
		combines <- c("max", "avg", "rcmax", "rcmax.avg", "BMA")
		if (!(object@ontology    %in% onts)) {
                    stop("*ontology* must be one of ", paste(onts, collapse=","))
		}
		if(length(object@organism != 0)) {
                    if (object@organism  %in% organs) {
                        ##loadAnnoPkg(object) ##load annotation pkg.
                    } else {
                        stop("*organism* must be one of ", paste(organs, collapse=","))
                    }
		}
		if (!(object@method      %in% mets)) {
                    stop("*method* must be one of ", paste(mets, collapse=","))
		}
		if(length(object@combine) != 0) {
                    if (!(object@combine %in% combines)) {
                        stop("*combine* must be one of ", paste(combines, collapse=","))
                    }
		}
		return (TRUE)
            }
            )

setMethod(
          f          = "computeIC",
          signature  = "Params",
          definition =function(params){
              gomap   <- get("gomap", envir=GOSemSimEnv)

              ## all GO terms appearing in an given ontology ###########
              goterms <- unlist(sapply(gomap, function(x) names(x[x == params@ontology])), use.names=FALSE)

              require(GO.db)
              if ( !exists("ALLGOID", envir=GOSemSimEnv) ) {
                  assign("ALLGOID", toTable(GOTERM), envir=GOSemSimEnv )
              }
              goids   <- get("ALLGOID", envir=GOSemSimEnv)
              ##goids <- toTable(GOTERM)

              ## all go terms which belong to the corresponding ontology..
              goids   <- unique(goids[goids[,"Ontology"] == params@ontology, "go_id"])
              gocount <- table(goterms)
              goname  <- names(gocount) #goid of specific organism and selected category.

              ## ensure goterms not appearing in the specific annotation have 0 frequency..
              go.diff        <- setdiff(goids, goname)
              m              <- double(length(go.diff))
              names(m)       <- go.diff
              gocount        <- as.vector(gocount)
              names(gocount) <- goname
              gocount        <- c(gocount, m)

              Offsprings.name <- switch(params@ontology,
                                        MF = "MFOffsprings",
                                        BP = "BPOffsprings",
                                        CC = "CCOffsprings"
                                        )
              if (!exists(Offsprings.name, envir=GOSemSimEnv)) {
                  .getOffsprings(params@ontology)
              }
              Offsprings <- get(Offsprings.name, envir=GOSemSimEnv)
              cnt        <- sapply(goids,function(x){ n=gocount[ Offsprings[[x]] ]; gocount[x]+sum(n[!is.na(n)])})
              names(cnt) <- goids
              ## the probabilities of occurrence of GO terms in a specific corpus.
              p          <- cnt/sum(gocount)
              ## IC of GO terms was quantified as the negative log likelihood.
              IC         <- -log(p)

              save(IC, file=paste(paste("Info_Contents", params@organism, params@ontology, sep="_"), ".rda", sep=""), compress="xz")
          }
          )


#########################################################
## setter functions for modifing parameters.
setReplaceMethod(
                 f          ="setOntology",
                 signature  ="Params",
                 definition =function(object, value){
                     object@ontology <- value
                     return(object)
                 }
                 )

setReplaceMethod(
                 f          ="setOrganism",
                 signature  ="Params",
                 definition =function(object, value){
                     object@organism <- value
                     ##loadAnnoPkg(object)
                     loadGOMap(object)  ## update GOMap to the corresponding species...
                     return(object)
                 }
                 )

setReplaceMethod(
                 f          ="setMethod",
                 signature  ="Params",
                 definition =function(object, value){
                     object@method <- value
                     return(object)
                 }
                 )

setReplaceMethod(
                 f          ="setCombineMethod",
                 signature  ="Params",
                 definition =function(object, value){
                     object@combine <- value
                     return(object)
                 }
                 )

############

setMethod(
          f                 ="[",
          signature         =signature(x = "Params",i  = "character"),
          definition        =function(x,i,j, ..., drop = TRUE) {
              if(i == "ontology")
                  return(x@ontology)
              if(i == "organism")
                  return(x@organism)
              if(i == "method")
                  return(x@method)
              if(i == "combine")
                  return(x@combine)
          }
          )


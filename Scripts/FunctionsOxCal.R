#####################################################################
#### Functions for grouping C14 dates based on material and taxa ####
#####################################################################

classify <- function(material, material.species, material.notes){
  p <- NA
  if(material %in% material_groups$short){
    p <- 1
    if(material %in% c("bone", 
                       "tooth", 
                      #  "charred plant macrofossils", 
                      #  "plant macrofossils",
                       "bone / tooth") && material.species %in% species_groups$medium){
      p <- 2
    }
  }
  else if(material %in% material_groups$medium){
    p <- 2
  }
  else if(material %in% material_groups$long){
    p <- 3
    if(material.species %in% species_groups$medium){
      p <- 2
    }
    if(is.na(material.species)){
      if(sum(unlist(lapply(species_groups$medium, function(x){
        grep(x, material.notes)
      }))) > 0){
        p <- 2
      }
      if(sum(unlist(lapply(species_groups$long, function(x){
        grep(x, material.notes)
      }))) > 0){
        p <- 3
      }
    }
  }
  return(p)
}

##########################################
#### Functions to prepare OxCal input ####
##########################################
Sequence <- function(name = NULL){
  prefix <- "\t"
  if(is.null(name)){
    seq <- paste(prefix, 'Sequence("%s"){')
  } else {
    seq <- paste(prefix, sprintf('Sequence("%s"){', name))
  }
  return(seq)
}

Boundary <- function(p, sequential = TRUE){
  b <- lapply(1:length(p), function(i){
    c(
      if (sequential | i == 1) paste(sprintf('Boundary("B_Start %s");', p[i])),
      if (!sequential & i > 1) paste(sprintf('Boundary("B_Transition %s/%s");', p[i-1], p[i])),
      if (sequential | i == length(p)) paste(sprintf('Boundary("B_End %s");', p[i]))
    )
  })
  return(b)
}

RDates <- function(d, outlier = FALSE, outlier.prob = NULL){
  if(outlier){
    ages <- sapply(1:nrow(d), function(r){
      if(outlier.prob[r] == 0){
        with(d[r,],
        sprintf('R_Date("%s",%d,%d);',
        LabID, C14.Age, C14.SD))
      } else {
        with(d[r,],
          sprintf('R_Date("%s",%d,%d) {Outlier(%d, "%s");};',
          LabID, C14.Age, C14.SD, OutlierProb, ifelse(is.na(TaxonCode), Material, TaxonCode)))
      }
    })
  } else {
    ages <- with(d,
    sprintf('R_Date("%s",%d,%d);', LabID, C14.Age, C14.SD)
    )
  }
  return(ages)
}

Curve <- function(calcurve = "IntCal20", calcurvepath = NULL){
  paste0('Curve("', calcurve, '","', calcurvepath, '");')
}

Phases <- function(p){
  paste(sprintf('Phase("%s"){', p))
}

Date <- function(p){
  paste(sprintf('Date("%s");', p))
}

OutlierModels <- function(d, exp.species, exp.material){
  require(data.table)
  if(sum(class(d) %in% "data.table") < 1){
    d <- setDT(d)
  }
  if(sum(class(exp.species) %in% "data.table") < 1){
    exp.species <- setDT(exp.species)
  }
  if(sum(class(exp.material) %in% "data.table") < 1){
    exp.material <- setDT(exp.material)
  }
  g <- d$Group
  if(any(g > 1)){
    sp <- d[Group > 1][, unique(TaxonCode)]
    om <- lapply(sp, function(x){
      if(!is.na(x)){
        m <- exp.species[TaxonCode == x]
        return(with(m, 
          sprintf('Outlier_Model("%s",Exp(%f,-%d,0),U(0,%d),"t");',
          TaxonCode, Tau, (Trunc/(10^Scale)), Scale)))
      } else {
        m <- exp.material[Material %in% d[Group > 1 & is.na(TaxonCode)][, unique(Material)]]
        return(with(m, 
          sprintf('Outlier_Model("%s",Exp(%f,-%d,0),U(0,%d),"t");',
          Material, Tau, (Trunc/(10^Scale)), Scale)))
      }
    })
    return(unlist(om))
  } else {
    return(NULL)
  }
}

OxCal_input <- function(curve, out.model, bound, phases, rdates, date_p, phase_type = NULL){

  inner <- lapply(1:length(phases), function(i){
    prefix <- "\t\t"
    paste(
      paste(prefix, bound[[i]][1]),
      paste(prefix, phases[[i]]),
      paste(prefix, "\t", rdates[[i]], collapse = "\n"),
      paste(prefix, "\t", date_p[[i]]),
      paste(prefix, '};'),
      if(length(bound[[i]]) > 1) paste(prefix, bound[[i]][2]),
      sep = "\n"
    )
  })

  prefix <- "\t"

  if(length(phase_type) > 0 && phase_type == 'overlapping'){
    overlaps <- unlist(lapply(inner, function(x){
        paste(
          paste(x, collapse = "\n"),
          sep = "\n"
        )
      }))
    command <- paste(
      'Plot(){',
      if(length(out.model) > 0) paste(prefix, out.model, collapse = "\n"),
      paste(prefix, curve),
      paste(prefix, 'Sequence(){'),
      paste(overlaps, collapse = "\n"),
      paste(prefix, '};'),
      '};',
      sep = "\n"
    )
  } else {
    command <- paste(
      'Plot(){',
      if(length(out.model) > 0) paste(prefix, out.model, collapse = "\n"),
      paste(prefix, curve),
      paste(prefix, 'Sequence(){'),
      paste(inner, collapse = "\n"),
      paste(prefix, '};'),
      '};',
      sep = "\n"
  )
  }
  return(command)
}

###############################################
#### Functions for extracting OxCal output ####
###############################################

get.values <- function(js, target, posterior = TRUE, numeric = TRUE){
  if(posterior){
    what <- "posterior"
  }else{
    what <- "likelihood"
    if(target == "agreement"){
      stop("agreement values are only available with posterior = TRUE")
    }
  }
  x <- js[grep(paste0(".",what,".",target), js)]
  val <- gsub(".+=(.+);", "\\1", x)
  lab <- gsub("^(ocd\\[\\d+\\]).*", "\\1", x)
  if(numeric){
    val <- as.numeric(val)
  } else {
    val <- as.character(val)
  }
  res <- data.frame(lab, val)
  setNames(res, c("OxCalID", target))
}

get.model.agreement <- function(js, target, posterior = TRUE, numeric = TRUE){
  if(posterior){
    what <- "posterior"
  } else {
    what <- "model"
  }
  if(target == "model"){
    x <- js[grep(paste0(".",what,".modelAgreement="), js)]
  } else if (target == "overall") {
     x <- js[grep(paste0(".",what,".overallAgreement="), js)]
  } else {
    stop("target can only be either \"model\" or \"overall\"")
  }
  val <- gsub(".+=(.+);", "\\1", x)
  if(numeric){
    val <- as.numeric(val)
  }
  return(val)
}

get.lab.codes <- function(js, outlier = FALSE){
  if(outlier){
    x <- js[grep(paste(c(" R_Date", " Outlier"), collapse = "|"), js)]
  } else {
    x <- js[grep(" R_Date", js)]
  }
  string <- gsub(".+=(.+);", "\\1", x)
  string <- gsub(" .*", "\\1", string)
  string <- gsub("[^A-Za-z0-9-]", "", string)
  return(string)
}

get.outlier <- function(js, prior = T){
  x <- js[grep("outlier_post", js)]
  if(length(x) > 0){
    val <- as.numeric(gsub(".+=(.+);", "\\1", x))
  if(prior){
    x <- js[grep("outlier_prior", js)]
    pr <- gsub(".+outlier_prior:(.+)", "\\1", x)
    pr <- as.numeric(gsub(",.+", "\\1", pr))
    return(data.frame(PosteriorProb = val, PriorProb = pr[1:length(val)]))
  }else{
    return(data.frame(PosteriorProb = val))
  }
  } else {return(NULL)}
}

get.probs <- function(js, posterior = TRUE, calBP = TRUE){
  if(posterior){
    what <- "posterior"
  }else{
    what <- "likelihood"
  }
  x <- js[grep(paste0(".",what,".prob="), js)]
  probs <- sapply(regmatches(x,gregexpr("\\[.+?\\]", x)), "[", 2)
  probs <- regmatches(probs,gregexpr("[[:digit:]]+\\.*[[:digit:]]*", probs))
  probs <- sapply(probs, as.numeric)
  ids <- gsub("^(ocd\\[\\d+\\]).*", "\\1", x)
  names(probs) <- ids
  # start <- get.values(js, "start", "OxCalID", posterior = posterior)
  start <- get.values(js, "start", posterior = posterior)
  # resolutions <- get.values(js, "resolution", "OxCalID", posterior = posterior)
  resolutions <- get.values(js, "resolution", posterior = posterior)
  # rmv <- c("ocd[1]", "ocd[2]", "ocd[3]")
  # start <- start[-which(start$OxCalID %in% rmv),]
  # resolutions <- resolutions[-which(resolutions$OxCalID %in% rmv),]
  dates <- lapply(seq_len(nrow(start)), function(r){
    pb <- probs[[start$OxCalID[r]]]
    date <- seq(start$start[r], by = resolutions$resolution[r], length.out = length(pb))
    data.frame(lab = start$OxCalID[r], 
                date = date, 
                prob = pb)
  })
  dates <- do.call(rbind, dates)
  if(calBP){
    dates$date <- 1950 - dates$date
  }
  return(dates)
}

get.results <- function(filepath, outpath, site_id = NULL, phase_id = NULL, phase_dating = TRUE, outlier = TRUE, output.list = FALSE, calBP = TRUE){
  require(data.table)
  require(stringr)
  require(gtools)
  ## Read results table
  results <- read.delim(
    list.files(filepath, pattern = paste0(site_id, ".txt"), full.names = T),
    header = F
  )
  results <- results[, 1:5]
  names(results) <- c("Variable", "From_68.3", "To_68.3", "From_95.4", "To_95.4")
  rows <-  grep("@", results$Variable)
  
  # Stop if modelling failed
  # if(length(rows) < 1) {cat(site_id); stop("Model failed")}
  if(length(rows) < 1) { 
    warning(paste ("Model", site_id, "failed"))
    return(NULL)
  } else {
    # Skip if there are no posterior values
  if(sum(is.na(results[rows, -1]))/(nrow(results[rows, -1])*ncol(results[rows,-1])) == 1) {
    return(NULL)
  } else {
    # Build from results table
  results$Type <- NA
  results[rows, "Type"] <- "Modelled"
  results[-(rows), "Type"] <- "Unmodelled"

  ## Read json file
  js <- readLines(list.files(filepath, pattern = paste0(site_id, ".js"), full.names = T))

  ## Extract parameter names
  vars <- results$Variable[grep("@", results$Variable)]
  if(length(grep(paste(c("Exp\\(", "U\\("), collapse = "|"), results$Variable)) > 0){
    miss <- grep(paste(c("Exp\\(", "U\\("), collapse = "|"), results$Variable, value = T)
    vars[which(vars == "@")] <- miss
  }
  vars <- gsub("@", "", vars)
  
  ## Create table for modelled data
  tables <- split(results, as.factor(results$Type))
  modelled <- tables$Modelled
  modelled$Variable <- vars
  # if(sum(is.na(modelled$To_95.4)) > 0){
  #   modelled <- modelled[which(is.na(modelled$To_95.4))[1]:nrow(modelled),]
  # }
  unmodelled <- setDT(tables$Unmodelled)
  rm(tables); gc()

  ## Get codes
  codes <- js[grep("^(ocd\\[\\d+\\]).op",  js)]
  unmodelled$OxCalID <- stringr::str_extract(codes, "([a-z]+\\[\\d+\\])")
  unmodelled$Variable_type <- gsub("\"", "", stringr::str_extract(codes, "(\"[A-z]+\")"))
  unmodelled <- unmodelled[!(Variable_type %in% c("Phase", "Curve", "Sequence")),]
  unmodelled <- as.data.frame(unmodelled)

  ix <- c("Outlier_Model", "Exp", "U")
  unmodelled$Parameter <- NA
  if(length(which(unmodelled$Variable_type %in% ix)) > 0){
    unmodelled[which(unmodelled$Variable_type %in% ix),]$Parameter <- "Model"
  }
  unmodelled$Parameter[is.na(unmodelled$Parameter)] <- "Date"

  # ix <- grep("Start", vars)[1]
  # names(vars) <- if(ix > 1){
  #   c(rep("Model", length(1:(ix-1))), rep("Date", length(ix:length(vars))))
  # } else {
  #   rep("Date", length(vars))
  # }

  ## Get values for modelled and unmodelled data
  oxcal_ids <- get.values(js, "mean", posterior = T)[, 1]
  modelled$OxCalID <- oxcal_ids
  modelled <- modelled[gtools::mixedorder(modelled$OxCalID),]
  if(sum(duplicated(modelled$OxCalID)) > 0){
    modelled <- modelled[!duplicated(modelled$OxCalID),]
  }

  # unmodelled <- unmodelled[which(unmodelled$Variable %in% modelled$Variable),]
  # unmodelled <- unmodelled[!duplicated(unmodelled$Variable),]
  # unmodelled <- unmodelled[order(match(unmodelled$Variable, modelled$Variable)),]
  # modelled$Parameter <- sapply(modelled$Variable, function(x) unique(names(which(vars == x))))
  modelled <- merge(modelled, unmodelled[, c("Variable_type", "Parameter", "OxCalID")], by = "OxCalID")
  modelled <- modelled[gtools::mixedorder(modelled$OxCalID),]
  # unmodelled$Parameter <- sapply(unmodelled$Variable, function(x) modelled[which(modelled$Variable == x), "Parameter"][1])
  # if(length(ix) > 1){
  #   unmodelled <- rbind(unmodelled[Parameter == "Model"], unmodelled[Parameter == "Date"][!duplicated(Variable)])
  # } else {
  #   unmodelled <- unmodelled[!duplicated(Variable)]
  # }

  ## Extract posterior values for the modelled data
  pars <- c("mean", "sigma", "median")
  posteriors <- lapply(pars, function(x){
    get.values(js, x)[,2]
  })
  names(posteriors) <- pars
  posteriors <- as.data.frame(do.call(cbind, posteriors))
  posteriors$OxCalID <- oxcal_ids
  posteriors <- posteriors[!duplicated(posteriors$OxCalID),]

  likelihood <- lapply(pars, function(x){
    get.values(js, x, posterior = F)[,2]
  })
  names(likelihood) <- pars
  likelihood <- as.data.frame(do.call(cbind, likelihood))
  likelihood$OxCalID <- get.values(js, "mean", posterior = F)[,1]

  ## Add to tables
  modelled <- merge(modelled, posteriors, by = "OxCalID")
  modelled <- modelled[gtools::mixedorder(modelled$OxCalID),]

  unmodelled <- merge(unmodelled, likelihood, by = "OxCalID", all.x = T)
  unmodelled <- unmodelled[gtools::mixedorder(unmodelled$OxCalID),]
  
  # modelled[, pars[1]] <- posteriors[[pars[1]]]
  # modelled[, pars[2]] <- posteriors[[pars[2]]]
  # modelled[, pars[3]] <- posteriors[[pars[3]]]

  # unmodelled[which(unmodelled$OxCalID %in% likelihood[[pars[1]]][,1]), pars[1]] <- likelihood[[pars[1]]][,2]
  # unmodelled[which(unmodelled$OxCalID %in% likelihood[[pars[2]]][,1]), pars[2]] <- likelihood[[pars[2]]][,2]
  # unmodelled[which(unmodelled$OxCalID %in% likelihood[[pars[3]]][,1]), pars[3]] <- likelihood[[pars[3]]][,2]

  ## Combine in one table
  # modelled <- modelled[, which(colnames(modelled) != "OxCalID")]
  combined <- rbind(unmodelled, modelled)
  combined <- combined[,-1]

  ## Split into a model table and a dates table
  tbl <- split(combined, as.factor(combined$Parameter))
  model_tbl <- tbl$Model
  dates_tbl <- tbl$Date
  rm(tbl); gc()
  dates_tbl$YearType <- "BCE"
  dates_tbl[, 1:ncol(dates_tbl)]  <- lapply(1:ncol(dates_tbl),function(x){
    tryCatch(
      {as.integer(dates_tbl[[x]])},
      warning = function(w){dates_tbl[[x]]}
      )
  })
  if(calBP){
    transcols <- c("From_68.3", "To_68.3", "From_95.4", "To_95.4", "mean", "median")
    # num <- which(sapply(1:ncol(dates_tbl), function(n){is.numeric(dates_tbl[,n])}))
    dates_tbl[, transcols]  <- lapply(transcols, function(x){
      1950 - dates_tbl[[x]]
    })
    dates_tbl$YearType <- "calBP"
  }

  ## Get model agreement
  dates_agreement <- get.values(js, "agreement")
  dates_agreement <- dates_agreement[!duplicated(dates_agreement$OxCalID),]
  dates_agreement$Variable <- modelled[which(modelled$OxCalID %in% dates_agreement$OxCalID), "Variable"]
  model_agreement <- rbind(dates_agreement[, -1],
                            data.frame(Variable = c("A_overall", "A_model"), 
                                      agreement = sapply(c("overall", "model"), 
                                      function(x){
                                        get.model.agreement(js, x,)
                                        })))
  names(model_agreement) <- gsub("agreement", "A", names(model_agreement))
  model_agreement <- model_agreement[, c(2,1)]
  rownames(model_agreement) <- NULL

  write.table(model_tbl, paste0(outpath, site_id, "_model.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(dates_tbl, paste0(outpath, site_id, "_dates.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(model_agreement, paste0(outpath, site_id, "_agreement.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

  if(phase_dating){
    if(is.null(phase_id)){
      stop("Phase IDs not provided")
    }
    phdate <- rbindlist(lapply(phase_id, function(pid){
      period <- dates_tbl[which(dates_tbl$Variable == pid & dates_tbl$Type == "Modelled"),]
      names(period)[1] <- "PhaseID"
      colremove <- c("Type", "Parameter", "Variable_type")
      colnames <- names(period)[-(which(names(period) %in% colremove))]
      period <- period[, colnames]
      period[, 1:ncol(period)] <- lapply(1:ncol(period),function(x) {
        tryCatch({
          as.integer(period[[x]])
          },warning = function(w) {
            period[[x]]}
        )}
      )
      return(period)
    }))
  write.table(phdate, paste0(outpath, site_id, "_phase.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  }

  if(outlier){
    outlier_tbl <- cbind(data.table(Outlier_Model = get.lab.codes(js, outlier = TRUE)),
      get.outlier(js))
    write.table(outlier_tbl, paste0(outpath, site_id, "_outlier.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  }

  if(output.list){
    return(lapply(list.files(outpath, pattern = site_id, full.names = T), read.delim))
  }
  }
  }
}

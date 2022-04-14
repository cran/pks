## Mar/30/2022: as.binmat() gains as.logical arg


## Convert binary matrix to vector of response patterns
as.pattern <- function(R, freq = FALSE, as.letters = FALSE, as.set = FALSE){
  if(freq){
    N.R <- table(apply(R, 1, paste, collapse=""))
    setNames(as.integer(N.R), names(N.R))          # convert to named int
  }else
    if(as.letters | as.set){
      nitems <- ncol(R)
      item.names <- 
       make.unique(c("a", letters[(seq_len(nitems) %% 26) + 1])[-(nitems + 1)],
                     sep="")
      lett <- apply(R, 1, function(r) paste(item.names[which(r == 1)],
                    collapse=""))
      lett[lett == ""] <- "0"

      if(as.set){
        # Separate elements in lett by "_", remove leading "_",
        # then strsplit along "_" (trailing "_" are ignored by strsplit)
        setfam <- as.set(lapply(strsplit(
          gsub("^_(.+)", "\\1", gsub("([0-9]*)", "\\1_", unname(lett))),
          "_"), as.set))
        if (set_contains_element(setfam, set("0")))
          setfam[[set("0")]] <- set()  # proper empty set
        setfam  # return family of sets, class set
      }else
        lett    # return letters, class character
    }else
      unname(apply(R, 1, paste, collapse=""))
}


## Convert vector of response patterns to named binary matrix
as.binmat <- function(N.R, uniq = TRUE, col.names = NULL, as.logical = FALSE){
  if (is.set(N.R)) {
    states <- lapply(N.R, as.character)
    items <- sort(unique(unlist(states)))
    R <- matrix(0, length(N.R), length(items),
                dimnames=list(NULL,
                              if(is.null(col.names)) items else col.names))
    for (i in seq_len(nrow(R))) R[i, states[[i]]] <- 1
  } else {
    pat <- if(is.null(names(N.R))) N.R else names(N.R)
    R   <- if(uniq) strsplit(pat, "") else strsplit(rep(pat, N.R), "")
    R   <- do.call(rbind, R)

    colnames(R) <- 
      if(is.null(col.names)){
        nitems <- ncol(R)
        make.unique(c("a", letters[(seq_len(nitems) %% 26) + 1])[-(nitems + 1)],
          sep="")
      }else
        col.names
  }
  storage.mode(R) <- "integer"
  if(as.logical) storage.mode(R) <- "logical"
  R
}


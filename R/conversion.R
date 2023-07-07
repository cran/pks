# Mar/30/2022 as.binmat() gains as.logical arg
# Jun/02/2023 useNames argument for as.pattern()
# Jun/06/2023 is.subset()


## Subset relation incidence matrix
is.subset <- function(R) {
  I <- t(apply(R, 1, function(r) apply(r * t(R) == r, 2, all)))
  if(!is.null(dimnames(I))) names(dimnames(I)) <- c("<", ">")
  I
}


## Convert binary matrix to vector of response patterns
as.pattern <- function(R, freq = FALSE, useNames = FALSE, as.set = FALSE,
                       sep = "", emptyset = "{}", as.letters = NULL) {
  if(!is.null(as.letters)) {
    warning("as.letters argument is deprecated, use useNames instead")
    useNames <- as.letters
  }
  if(freq) {
    N.R <- table(apply(R, 1, paste, collapse=""))
    setNames(as.integer(N.R), names(N.R))          # convert to named int
  } else
    if(useNames | as.set) {
      nitems <- ncol(R)
      item.names <- colnames(R)
      if(is.null(item.names))
        item.names <- make.unique(
          c("a", letters[(seq_len(nitems) %% 26) + 1])[-(nitems + 1)],
          sep = ""
        )
      lett <- unname(apply(R == TRUE, 1, function(r) item.names[r]))

      if(as.set) {
        as.set(sapply(lett, as.set))  # return family of sets, class set
      } else {
        lett <- sapply(lett, paste, collapse = sep)
        lett[lett == ""] <- emptyset
        lett                          # return letters, class character
      }
    } else
      unname(apply(R, 1, paste, collapse = ""))
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


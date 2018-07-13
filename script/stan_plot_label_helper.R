.make_plot_data <- function(object, pars, include = TRUE,
                            inc_warmup = FALSE, unconstrain = FALSE) {
  
  window <- NULL
  if (is.stanreg(object)) {
    sim <- object$stanfit@sim
  }
  else sim <- object@sim
  
  nopars <- missing(pars)
  if (is.stanfit(object) && !nopars) {
    if ("log-posterior" %in% pars)
      pars[which(pars == "log-posterior")] <- "lp__"
  }
  if (!include) {
    if (nopars) stop("pars must be specified if include=FALSE.", call. = FALSE)
    else {
      if (is.stanreg(object)) 
        pars <- setdiff(sim$fnames_oi, pars)
      else 
        pars <- setdiff(sim$pars_oi, pars)
    }
  }
  if (nopars) {
    if (is.stanreg(object)) 
      pars <- names(object$coefficients)
    else 
      pars <- setdiff(sim$pars_oi, c("lp__", "log-posterior"))
  }
  else {
    if (!is.stanreg(object)) 
      pars <- check_pars_second(sim, pars)
  }
  
  pars <- remove_empty_pars(pars, sim$dims_oi)
  if (unconstrain && "lp__" %in% pars) {
    if (length(pars) == 1L) stop("Can't unconstrain lp__", call. = FALSE)
    pars <- pars[-which(pars == "lp__")]
  }
  tidx <- pars_total_indexes(sim$pars_oi,
                             sim$dims_oi,
                             sim$fnames_oi,
                             pars)
  tidx <- lapply(tidx, function(x) attr(x, "row_major_idx"))
  tidx <- unlist(tidx, use.names = FALSE)
  num_plots <- length(tidx)
  
  if (nopars && num_plots > 10) {
    # if pars is not specified then default to showing 10 of the parameters
    tidx <- tidx[1:10]
    message("'pars' not specified. Showing first 10 parameters by default.")
  }
  
  if (!is.null(window)) {
    window <- sort(window)
    if (window[1] < 1) window[1] <- 1
    if (window[1] > sim$iter[1])
      stop("wrong specification of argument window", call. = FALSE)
    if (is.na(window[2]) || window[2] > sim$iter[1])
      window[2] <- sim$iter[1]
  } else {
    window <- c(1, sim$iter[1])
  }
  if ((sim$warmup2 == 0 || !inc_warmup) && window[1] <= sim$warmup[1]) {
    window[1] <- sim$warmup[1] + 1
  }
  if (window[1] > window[2]) {
    stop("the given window does not include sample")
  }
  if (window[1] > sim$warmup[1]) inc_warmup <- FALSE
  
  thin <- sim$thin
  warmup2 <- sim$warmup2[1]
  warmup <- sim$warmup
  start_i <- window[1]
  window_size <- (window[2] - start_i) %/% thin
  id <- seq.int(start_i, by = thin, length.out = window_size)
  start_idx <- (if (warmup2 == 0) (start_i - warmup) else start_i) %/% thin
  if (start_idx < 1)  start_idx <- 1
  idx <- seq.int(start_idx, by = 1, length.out = window_size)
  
  if (unconstrain) {
    sim$samples <- .upars(object)
    sel <- grep(paste(pars, collapse ="|"), names(sim$samples[[1]]), value = TRUE)
    samp_use <- lapply(sim$samples, function(chain) {
      out <- lapply(chain[names(chain) %in% sel], function(x) x[idx])
      names(out) <- sel
      out
    })
  } else {
    samp_use <- lapply(sim$samples, function(chain) {
      out <- lapply(chain[tidx], function(x) x[idx])
      names(out) <- sim$fnames_oi[tidx]
      out
    })
  }
  nchains <- length(samp_use)
  
  if (unconstrain) {
    if (is.stanreg(object)) 
      object <- object$stanfit
    pblock <- .get_stan_params(object)
    pars2 <- unlist(lapply(strsplit(pars, "\\["), "[[", 1))
    not_pblock <- length(setdiff(unique(pars2), pblock))
    if (not_pblock)
      stop("If 'unconstrain' is TRUE only variables declared in the ",
           "'parameters' block can be included in 'pars'.", 
           call. = FALSE)
  }
  
  dat <- .reshape_sample(samp_use)
  dat$iteration <- idx
  dat$chain <- factor(dat$chain)
  fnames <- if (unconstrain) 
    names(samp_use[[1]]) else sim$fnames_oi[tidx]
  lp <- which(dat$parameter == "lp__")
  if (!identical(lp, integer(0))) {
    dat$parameter[lp] <- "log-posterior"
    fnames[which(fnames == "lp__")] <- "log-posterior"
  }
  dat$parameter <- factor(dat$parameter, levels = fnames)
  list(samp = dat,
       nchains = nchains,
       nparams = length(fnames),
       warmup = warmup)
}
.ac_plot_data <- function(dat, lags, partial = FALSE) {
  ch <- dat[, grep("chain", colnames(dat))]
  nc <- length(unique(ch))
  ac_list <- tapply(dat$value, INDEX = ch, FUN = .ac_fun, lag.max = lags,
                    partial = partial, simplify = FALSE)
  nl <- if (partial) lags else lags + 1
  ch <- factor(rep(1:nc, each = nl), labels = paste0("chain:", 1:nc))
  ll <- rep(seq(if (partial) 1 else 0, lags), nc)
  data.frame(chains = ch, ac = do.call(c, ac_list), lag = ll)
}
.ac_plot_data_multi <- function(dat, lags, partial = FALSE) {
  ch <- dat[, grep("chain", colnames(dat))]
  nc <- length(unique(ch))
  pa <- factor(dat[, grep("parameter", colnames(dat))])
  np <- length(unique(pa))
  ac_list <- tapply(dat$value, INDEX = list(ch, pa),
                    FUN = .ac_fun, lag.max = lags,
                    partial = partial, simplify = FALSE)
  nl <- if (partial) lags else lags + 1
  ch <- factor(rep(rep(1:nc, each = nl), np), labels = paste0("chain:", 1:nc))
  ll <- rep(seq(if (partial) 1 else 0, lags), nc * np)
  pp <- factor(rep(1:np, each = nc * nl), labels = levels(pa))
  data.frame(parameters = pp, chains = ch, ac = do.call(c, ac_list), lag = ll)
}
is.stanreg <- function(x) inherits(x, "stanreg")
is.stanfit <- function(x) inherits(x, "stanfit")
check_pars_second <- function(sim, pars) {
  if (missing(pars)) return(sim$pars_oi)
  allpars <- c(sim$pars_oi, sim$fnames_oi)
  check_pars(allpars, pars)
}

# autocorrelation ---------------------------------------------------------
.ac_fun <- function(x, lag.max, partial = FALSE) {
  if (!partial)
    acf(x, lag.max = lag.max, plot = FALSE)$acf[,, 1L]
  else
    pacf(x, lag.max = lag.max, plot = FALSE)$acf[,, 1L]
}
.ac_plot_data <- function(dat, lags, partial = FALSE) {
  ch <- dat[, grep("chain", colnames(dat))]
  nc <- length(unique(ch))
  ac_list <- tapply(dat$value, INDEX = ch, FUN = .ac_fun, lag.max = lags,
                    partial = partial, simplify = FALSE)
  nl <- if (partial) lags else lags + 1
  ch <- factor(rep(1:nc, each = nl), labels = paste0("chain:", 1:nc))
  ll <- rep(seq(if (partial) 1 else 0, lags), nc)
  data.frame(chains = ch, ac = do.call(c, ac_list), lag = ll)
}
.ac_plot_data_multi <- function(dat, lags, partial = FALSE) {
  ch <- dat[, grep("chain", colnames(dat))]
  nc <- length(unique(ch))
  pa <- factor(dat[, grep("parameter", colnames(dat))])
  np <- length(unique(pa))
  ac_list <- tapply(dat$value, INDEX = list(ch, pa),
                    FUN = .ac_fun, lag.max = lags,
                    partial = partial, simplify = FALSE)
  nl <- if (partial) lags else lags + 1
  ch <- factor(rep(rep(1:nc, each = nl), np), labels = paste0("chain:", 1:nc))
  ll <- rep(seq(if (partial) 1 else 0, lags), nc * np)
  pp <- factor(rep(1:np, each = nc * nl), labels = levels(pa))
  data.frame(parameters = pp, chains = ch, ac = do.call(c, ac_list), lag = ll)
}
.reshape_sample <- function(x) {
  res <- lapply(seq_along(x), function(i) {
    data.frame(value = unlist(x[[i]], use.names = FALSE),
               parameter = rep(names(x[[i]]), each = length(x[[i]][[1L]])),
               chain = i)
  })
  res <- do.call(rbind, res)
  res$parameter <- as.character(res$parameter)
  res
}
# defaults ----------------------------------------------------------------
.rstanvis_defaults <- new.env(parent = emptyenv())
.rstanvis_defaults$theme <-
  theme_classic(base_size = 11) +
  theme(axis.line = element_line(color = "#222222"),
        axis.line.y = element_line(size = .5),
        axis.line.x = element_line(size = 1),
        axis.title = element_text(face = "bold", size = 13),
        strip.background = element_blank(),
        strip.text = element_text(color = "black", face = "bold"),
        legend.title = element_text(size = 11),
        plot.title = element_text(size = 18))

for (j in seq_along(.rstanvis_defaults$theme)) {
  if ("element_text" %in% class(.rstanvis_defaults$theme[[j]])) {
    .rstanvis_defaults$theme[[j]][["debug"]] <- FALSE
  }
}

.rstanvis_defaults$hist_theme <-
  .rstanvis_defaults$theme +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank())

.rstanvis_defaults$multiparam_theme <-
  .rstanvis_defaults$theme +
  theme(axis.title = element_blank(),
        axis.text.y = element_text(face = "bold", size = 13, debug = FALSE),
        legend.position = "none",
        panel.grid.major =  element_line(size = 0.1, color = "darkgray"))

.rstanvis_defaults$pt_color <- "#B2001D"
.rstanvis_defaults$alpha <- 0.5
.rstanvis_defaults$shape <- 19
.rstanvis_defaults$fill <-  "#B2001D"
.rstanvis_defaults$color <- "black" #"#590815"
.rstanvis_defaults$size <- 0.5
.rstanvis_defaults$pt_size <- 3

.rstanvis_defaults$chain_colors <- rgb(matrix(c(230, 97, 1,
                                                153, 142, 195,
                                                84, 39, 136,
                                                241, 163, 64,
                                                216, 218, 235,
                                                254, 224, 182),
                                              byrow = TRUE, ncol = 3),
                                       names = paste(1:6), maxColorValue = 255)
.rstanvis_defaults$grays <- rgb(matrix(c(247, 247, 247, 204, 204, 204,
                                         150, 150, 150, 82, 82, 82),
                                       byrow = TRUE, ncol = 3),
                                alpha = 100, names = paste(1:4), maxColorValue = 255)
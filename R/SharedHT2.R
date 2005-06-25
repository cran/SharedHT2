"EB.Anova" <-
function(data, labels, H0 = "equal.means", Var.Struct = "general",
         verbose = TRUE, subset, theta0 = NULL, gradient = FALSE,
         fit.only = FALSE)
{
  m <- .call. <- match.call()
  if(missing(Var.Struct)) {
    Var.Struct <- "general"
    .call.$Var.Struct <- Var.Struct
  }
  vs <- charmatch(Var.Struct, c("simple", "general"), 0)
  if(vs == 0)
    stop("Argument 'Var.Struct' must be \"simple\" or \"general\"")
  cat("Var.Struct = ",Var.Struct,"\n")
  m[[1]] <- as.name(c("SharedVar","SharedHT2")[vs])
  nms <- names(data)  
  clnum <- length(labels)
  d <- clnum
  n <- length(grep(labels[1], nms))
  if(is.numeric(H0)){
    M <- m$M <- rbind(H0)
    H0 <- m$H0 <- "user"
    r <- qr(M)$rank
    d.M <- dim(M)
    if(d.M[2] != d | r > d)
      stop("User supplied contrast matrix must be of dimension " %,%
           "q x (#groups) where q <= (#groups) and of rank <= (#groups)")
  }
  h0 <- charmatch(H0, c("no.trend","equal.means","zero.means"), 0)
  if(h0 == 0 & H0 != "user")
    stop("Argument 'H0' must be one of c(\"zero.means\",\"equal.means\",\"no.trend\")")
  if(vs==2){
    if(h0 == 3 & n <= d)
      stop("Var.Struct \"general\" and H0=\"zero.means\" requires that " %,%
           "(#reps) > (#groups)")
    if(h0 == 2 & n <= d-1)
      stop("Var.Struct \"general\" and H0=\"equal.means\" requires that " %,%
           "(#reps) > (#groups)-1")
    if(h0 == 1 & n <= 1)
      stop("H0=\"no.trend\" requires at least two replicates per group")
  }
  if(vs==1 && n <= 1)
    stop("Var.Struct==\"simple\" requires at least two replicates per group")
 
  m$Var.Struct <- NULL
  ans <- eval(m)
  if(!fit.only) 
    ans$EBfit$call <- .call.
  else ans$call <- .call.
  ans
} 

"SharedHT2" <- 
function(data, labels, H0 = "equal.means", M = NULL, verbose = TRUE, subset,
         theta0 = NULL, gradient = FALSE, fit.only = FALSE)
{
  nms <- names(data)
  m <- .call. <- match.call()
  m[[1]] <- as.name("model.frame")
  idx <- c(sapply(labels,FUN=grep, nms))
  form <- make.form("", nms[idx])
  m$formula <- form
  m$labels <- m$H0 <- m$M <- m$verbose <- m$theta0 <- m$gradient <- m$fit.only <- NULL
  m <- eval(m, sys.parent())
  id <- dimnames(m)[[1]]
  gnum <- dim(m)[1]
  clnum <- length(labels)
  N <- gnum
  p <- clnum
  n <- length(grep(labels[1], nms))
  mu.g <- matrix(NA, N, p)
  var.g <- matrix(NA, N, p^2)
  #find within group mean variance and sample size
  all.cols <- NULL
  Terms <- attr(m, "terms")
  Y <- model.matrix(Terms, m)[, -1, drop = F]
  for(i in 1:p) {
    mu.g[, i] <- rowMeans(Y[, n * (i - 1) + (1:n)], na.rm = T)
    for(j in unique(1:i)) {
      Z.i <- Y[, n * (i - 1) + (1:n)] - mu.g[, i]
      Z.j <- Y[, n * (j - 1) + (1:n)] - mu.g[, j]
      var.g[, p * (i - 1) + j] <- var.g[, p * (j - 1) + i] <-
        n/(n - 1) * rowMeans(Z.i * Z.j, na.rm = T)
    }
  }
  if(H0 == "zero.means") {
    M <- diag(p)
    d <- p
  }  
  if(H0 == "equal.means") {
    b <- diag(p) - 1/p * matrix(1, p, p)
    b.5 <- chol(b)
    b.5 <- b.5[ - p,  ]
    d <- p - 1
    M <- rbind(b.5)
  }
  if(H0 == "no.trend") {
    x <- cbind(rep(1, p), (1:p))
    PI <- solve(t(x) %*% x) %*% t(x)
    M <- PI[2,  , drop = F]
    d <- 1
  }
  if(H0 == "user") d <- qr(M)$rank
  if(verbose) {
    cat("You have chosen H0: ", H0, "\n")
    cat("N = ", N, " d = ", d, "\n")
  }
  Ybar <- matrix(mu.g %*% t(M), N, d)
  MVM <- t((M %K% diag(d)) %*% ((diag(p) %K% M) %*% t(var.g)))   
  dimMVM <- dim(MVM)
  npar <- d*(d+1)/2 + 1
  if(missing(theta0)) theta0 <- rep(0, npar)
  nreps <- rep(n, N)
  fit <- .C("FitInvWish",
            ptheta0 = as.double(theta0),
            MVM = as.double(t(MVM)),
            pN = as.integer(N),
            pd = as.integer(d),
            pnreps = as.integer(nreps),
            pverbose = as.integer(verbose),
            objval = as.double(0),
            estimate = as.double(rep(0,npar)),
            fail = as.integer(0),
            fncnt = as.integer(0),
            grcnt = as.integer(0),
            mask = as.integer(0),
            usegr = as.integer(gradient),
            G = as.double(rep(0,npar)),
            H = as.double(rep(0,npar^2)),
            PACKAGE = "SharedHT2")
#	V.theta <- matrix(0, npar, npar)
  theta <- fit$estimate
  V.theta <- solve(matrix(fit$H,npar,npar))
  nu <- (2*d + 2)*(exp(theta[1])+1)
  if(d > 1) {
    Lambda.5 <- matrix(0, d, d)
    diag(Lambda.5) <- exp(theta[2:(1 + d)])
    Lambda.5[outer(1:d, 1:d, FUN = "<")] <- theta[ - (1:(1 + d))]
    Lambda <- t(Lambda.5) %*% Lambda.5
    Lambda.u <- c(Lambda[outer(1:d, 1:d, FUN = "<=")])
  }
  if(d == 1) Lambda <- Lambda.u <- matrix(exp(theta[-1]), 1, 1)
  if(!fit.only){
    Lambda.N <- t(matrix(c(Lambda), d^2, N))
    V <- (n - 1) * MVM
    V.i <- matrix(rowInvc(V), N, d^2)
    t1 <- matrix(rowProdc(Ybar, V.i, 1, d), N, d)
    HT2.stat <- (n - d + 1)/d * n * rowProdc(t1, Ybar, 1, d)
    HT2.pval <- 1 - pf(HT2.stat, d, n - d + 1)
    V <- Lambda.N + (n - 1) * MVM
    V.i <- matrix(rowInvc(V), N, d^2)
    t1 <- matrix(rowProdc(Ybar, V.i, 1, d), N, d)
    ShHT2.stat <- (nu + n - 2 * d - 1)/d * n * rowProdc(t1, Ybar, 1, d)
    ShHT2.pval <- 1 - pf(ShHT2.stat, d, nu + n - 2 * d - 1)
    ans <- list()
    ans$data <- as.data.frame(list(GeneId = id, ShHT2.stat = ShHT2.stat,
                              ShHT2.pval = ShHT2.pval, HT2.stat = HT2.stat,
                              HT2.pval = HT2.pval))
    .call.$theta0 <- .call.$gradient <- NULL
    stat.names <- c("ShHT2","HT2")
    names(stat.names) <- c("EB", "naive")
    fit <- list(coefficients = theta, variance = V.theta, gradient = fit$G,
                shape = nu, rate = Lambda, log.likelihood = -1*fit$objval,
                fail = fit$fail, fncount = fit$fncount, grcount = fit$grcount,
                mask = fit$mask, usegr = fit$usegr, call=.call.,
                stat.names=stat.names)
    class(fit) <- "EBfit"
    ans$EBfit <- fit
    class(ans) <- "fit.n.data"
  }
  else {
    fit <- list(coefficients = theta, variance = V.theta, gradient = fit$G,
                shape = nu, rate = Lambda, log.likelihood = -1*fit$objval,
                fail = fit$fail, fncount = fit$fncount, grcount = fit$grcount,
                mask = fit$mask, usegr = fit$usegr, call=.call.)
    class(fit) <- "EBfit"
    ans <- fit
  }
  ans
}

"SharedVar" <- 
  function(data, labels, H0="zero.means", M = NULL, verbose = TRUE, subset, theta0 = NULL,
           gradient = FALSE, fit.only = FALSE)
{
  nms <- names(data)
  m <- .call. <- match.call()
  m[[1]] <- as.name("model.frame")
  idx <- c(sapply(labels,FUN=grep, nms))
  form <- make.form("", nms[idx])
  m$formula <- form
  m$labels <- m$H0 <- m$M <- m$verbose <- m$theta0 <- m$gradient <- m$fit.only <- NULL
  m <- eval(m, sys.parent())
  id <- dimnames(m)[[1]]
  gnum <- dim(m)[1]
  clnum <- length(labels)
  N <- gnum
  d <- clnum
  n <- length(grep(labels[1], nms))
  mu.g <- matrix(NA, N, d)
  var.g <- matrix(NA, N, d^2)
  #find within group mean variance and sample size
  all.cols <- NULL
  Terms <- attr(m, "terms")
  Y <- model.matrix(Terms, m)[, -1, drop = F]
  S <- rep(0, N)
  for(i in 1:d) {
    mu.g[, i] <- rowMeans(Y[, n * (i - 1) + (1:n)], na.rm = T)
    Z.i <- Y[, n * (i - 1) + (1:n)] - mu.g[, i]
    S <- S + n * rowMeans(Z.i^2, na.rm = T)
  }

  if(missing(theta0)) theta0 <- rep(0, 2)
  nreps <- rep(n, N)
  
  fit <- .C("FitEqualVar",
            ptheta0 = as.double(theta0),
            S = as.double(S),
            pN = as.integer(N),
            pd = as.integer(d),
            pnreps = as.integer(nreps),
            pverbose = as.integer(verbose),
            objval = as.double(0),
            estimate = as.double(rep(0,2)),
            fail = as.integer(0),
            fncnt = as.integer(0),
            grcnt = as.integer(0),
            mask = as.integer(0),
            usegr = as.integer(gradient),
            G = as.double(rep(0,2)),
            H = as.double(rep(0,4)),
            PACKAGE = "SharedHT2")
  theta <- fit$estimate
  V.theta <- solve(matrix(fit$H,2,2))
  eth <- exp(theta)
  s <- eth[1]
  r <- eth[2]
  J.n.n <- matrix(1, n, n)
  J.d.d <- matrix(1, d, d)
  J.d.1 <- matrix(1, d, 1)
  if(!fit.only){
    if(H0 == "zero.means") M <- diag(d)
    if(H0 == "equal.means") {
      b <- diag(d) - 1/d * matrix(1, d, d)
      b.5 <- chol(b)
      b.5 <- b.5[ - d,  ]
      M <- rbind(b.5)
    }
    if(H0 == "no.trend") {
      x <- cbind(rep(1, d), (1:d))
      PI <- solve(t(x) %*% x) %*% t(x)
      M <- (J.d.1/d) %K% (PI[2,  , drop = F])
    }
  
    rnk <- qr(M)$rank
    degTop <- ifelse(rnk == d, n*d, rnk)
    M <- M %K% J.n.n
    tr.inds <- c(outer(n*(0:(d-1)), 1:n, FUN="+"))
    Top <- ((Y[,tr.inds] %*% M)*Y) %*% rep(1, n*d)
    ShUT2.stat <- Top/(2*r + S)*(2*s+d*(n-1))/degTop
    ShUT2.pval <- 1-pf(ShUT2.stat, degTop, 2*s+d*(n-1))
    UT2.stat <- Top/S*(d*(n-1))/degTop
    UT2.pval <- 1-pf(UT2.stat, degTop, d*(n-1))
    ans <- list()
    ans$data <- as.data.frame(list(GeneId = id, ShUT2.stat = ShUT2.stat,
                              ShUT2.pval = ShUT2.pval, UT2.stat = UT2.stat,
                              UT2.pval = UT2.pval))
    .call.$theta0 <- .call.$gradient <- NULL
    stat.names <- c("ShUT2","UT2")
    names(stat.names) <- c("EB", "naive")
    fit <- list(coefficients = theta, variance = V.theta, gradient = fit$G,
                shape=s, rate=r, log.likelihood = -1*fit$objval, fail = fit$fail,
                fncount = fit$fncount, grcount = fit$grcount, mask = fit$mask,
                usegr = fit$usegr, call=.call., stat.names = stat.names)
    class(fit) <- "EBfit"
    ans$EBfit <- fit
    class(ans) <- "fit.n.data"
  }
  else{
    fit <- list(coefficients = theta, variance = V.theta, gradient = fit$G,
                shape=s, rate=r, log.likelihood = -1*fit$objval, fail = fit$fail,
                fncount = fit$fncount, grcount = fit$grcount, mask = fit$mask,
                usegr = fit$usegr, call=.call.)
    class(fit) <- "EBfit"
    ans <- fit
  }
  ans
}

"TopGenes" <- 
function (obj, by = "EB", ref = 1, FDR = 0.05, allsig = FALSE, n.g = 20, 
          browse = FALSE, search.url = genecards, path = "", file = "") 
{
  if(class(obj)!="fit.n.data")
    stop("TopGenes expects an argument of class 'fit.n.data' which is returned by EB.Anova\n" %,%
         "make sure the argument 'fit.only' is set to FALSE (default)")
  if(by != "EB" && by != "naive")
    stop("Argument 'by' must be set to either \"EB\" (Empirical Bayes) or \"naive\"")
  fit <- obj$EBfit
  m <- obj.call <- fit$call
  this.call <- match.call()
  obj <- obj$data
  nosiggenes <- F
  N <- dim(obj)[1]
  stepwise <- (FDR * (1:N))/N
  m$id <- NULL
  labs <- eval(m$labels)
  labs <- substring(labs, 5, nchar(labs))
  m[[1]] <- as.name("means.var")
  m$H0 <- m$Var.Struct <- NULL
  m <- eval(m, sys.parent())
  mu <- m$mean
  id <- as.character(obj$GeneId)
  stat.names <- fit$stat.names
  stat.idx <- which(names(obj) == stat.names[by] %,% ".stat")
  p.idx <- which(names(obj) == stat.names[by] %,% ".pval")
  indx <- order(-obj[, stat.idx])
  position <- (1:N)[indx]
  if (allsig) {
    sig.inds <- obj[, p.idx][indx] <= stepwise
    if(sum(sig.inds)==0) n.g <- 0
    else n.g <- max(which(obj[, p.idx][indx] <= stepwise))
  }
  if (n.g == 0) 
    nosiggenes <- T
  if (!nosiggenes) {
    position <- position[1:n.g]
    g <- id[indx[1:n.g]]
    mu <- mu[indx[1:n.g], ]
    T2 <- obj[, stat.idx][indx[1:n.g]]
    pval <- obj[, p.idx][indx[1:n.g]]
  }
  p <- dim(rbind(mu))[2]
  p. <- switch(missing(ref) + 1, 2 * p + 3, 4 + p)
  if (!nosiggenes) {
    stepwise <- stepwise[1:n.g]
    ans <- cbind(1:length(g), 1:length(g), rbind(signif(2^mu, 3)))
    if (!missing(ref)) {
      mu <- rbind(mu)
      ratio <- signif(2^(mu[, -ref, drop = F] - mu[, ref]), 3)
      ans  <- cbind(ans, signif(ratio, 3), signif(T2, 3), signif(pval, 4),
                    signif(stepwise, 4))
    }
    else 
      ans <- cbind(ans, signif(T2, 3), signif(pval, 4), signif(stepwise, 4))
    msg <- ifelse(allsig,
                  "B-H FDR procedure detects " %,% n.g %,% " " %,%
                  "significant genes at FDR=" %,% FDR %,% ":\n",
                  "Top " %,% n.g %,% " Genes: \n")
  }
  else {
    ans <- rbind(rep(NA, p.))
    g <- NA
    position <- "<NA>"
    msg <- "B-H FDR procedure detects no significant genes at FDR=" %,% FDR %,% ".\n"
  }
  ans <- as.data.frame(ans)
  ans[[1]] <- position
  ans[[2]] <- g
  if (!missing(ref)) {
    labs.ind <- (1:2) * (ref == 1) + (2:1) * (ref == 2)
    if (ref == 1) 
      labs.ind <- 1:p
    if (ref == 2) 
      labs.ind <- c(2, 1, 3:p)
    if (ref > 2 && ref < p) 
      labs.ind <- c(ref, unique(1:(ref - 1)), unique((ref + 1):p))
    if (ref == p) 
      labs.ind <- c(p, unique(1:(p - 1)))
    nms <- c("RowNum", "GeneId", labs, labs[labs.ind[-1]] %,% 
             "/" %,% labs[labs.ind[1]], stat.names[by] %,% ".stat",
             stat.names[by] %,% ".p-val", "FDR.stepdown=" %,% FDR)
  }
  else nms <- c("RowNum", "GeneId", labs, stat.names[by] %,% ".stat",
                stat.names[by] %,% ".p-val", "FDR.stepdown=" %,% FDR)
  names(ans) <- nms
  if (browse) {
    cmd <- as.call(expression(htmtbl))
    cmd$x <- as.name("ans")
    if (!missing(file)) 
      cmd$file <- file
    if (!missing(path)) 
      cmd$path <- path
    cmd$search.url <- search.url
    cmd$id <- as.name("GeneId")
    eval(as.call(cmd))
    cat("Genelist linked to the \"GeneCards\" database loading in " %,%
        getOption("browser") %,% "...\n")
  }
  cat(msg)
  ans
}

"shape.rate" <- 
function(object)
{
  theta <- NULL
  if(class(object)=="fit.n.data")
    theta <- EBfit(object)$coefficients
  if(class(object)=="EBfit")
    theta <- object$coefficients
  if(is.null(theta)) theta <- object
  if(!is.numeric(theta))
    stop("argument must be a numeric vector, of class fit.n.data or of class EBfit")
  p <- length(theta)
  d <- ((8*(p-1)+1)^0.5-1)/2
  if(floor(d)!=d)
    stop("numeric argument not of appropriate length")
  nu <- (2*d + 2)*(exp(theta[1])+1)
  L.5 <- matrix(0, d, d)
  diag(L.5) <- exp(theta[2:(d+1)])
  L.5[outer(1:d,1:d,FUN="<")] <- theta[-(1:(d+1))]
  L <- t(L.5)%*%L.5
  list(shape = nu, rate = L)
}

"tlmgamma" <- 
function(x, d)
{
  ans <- .C("tlmgamma",
            px = as.double(x),
            pd = as.integer(d),
            pans = as.double(0),
            PACKAGE="SharedHT2")
  ans$pans
} 

"tdet" <- 
function(x)
{
  d <- dim(x)[1]
  ans <- .C("tdet",
            x = as.double(x),
            xd2buf = as.double(rep(0,d^2)),
            pd = as.integer(d),
            pans = as.double(0),
            PACKAGE="SharedHT2")
  ans$pans
}

"tChol" <- 
function(x){
  d <- dim(x)[1]
  ans <- .C("chol",
            s = as.double(x),
            t = as.double(rep(0, d^2)),
            pd = as.integer(d),
            PACKAGE="SharedHT2")
  matrix(ans$t, d, d)
}

"tloglik" <- 
function(theta, MVM, nreps)
{
	dimMVM <- dim(MVM)
	N <- dimMVM[1]
	d <- dimMVM[2]^0.5
	npar <- d*(d+1)/2 + 1
	if(npar!=length(theta)) stop("'theta' should be of length " %,% npar)

	ans <- .C("tloglik",
		ptheta=as.double(theta),
		MVM=as.double(t(MVM)),
		pN = as.integer(N),
		pd = as.integer(d),
		pnreps = as.integer(nreps),
		pans = as.double(0),
                PACKAGE="SharedHT2")
	ans$pans
}

"tGloglik" <-
function(theta, MVM, nreps)
{
	dimMVM <- dim(MVM)
	N <- dimMVM[1]
	d <- dimMVM[2]^0.5
	npar <- d*(d+1)/2 + 1
	if(npar!=length(theta)) stop("'theta' should be of length " %,% npar)
	ans <- .C("tGloglik",
                  ptheta = as.double(theta),
                  MVM = as.double(t(MVM)),
                  pN = as.integer(N),
                  pd = as.integer(d),
                  pnreps = as.integer(nreps),
                  pG = as.double(rep(0, npar)),
                  PACKAGE="SharedHT2")
	ans$pG
}

"genecards" <- "http://thr.cit.nih.gov/cgi-bin/cards/cardsearch.pl?search="

"htmtbl" <- 
function(x, path = "", file = "", main = "", id = NULL, search.url = NULL,
         browser = getOption("browser"))
{
	call. <- match.call()
	id.nm <- as.character(call.$id)
	id.idx <- 0
	if(!is.null(id.nm))
		id.idx <- grep(id.nm, names(x))
	if(file == "") {
		file <- tempfile("temp")
	}
	else {
		if(path == "")
			path <- system("pwd", intern=T)
		file <- path %,% "/" %,% file
	}
	table.attributes <- "BORDER"
	row.header.column.attributes <- ""
	data.column.attributes <- ""
	column.attributes <- c(row.header.column.attributes, rep(
		data.column.attributes, length = ncol(x)))
	row.header.attributes <- "ALIGN=RIGHT"
	cell.attributes <- ""
	if(!missing(id))
		id.vals <- as.character(x[, id.idx])
	out <- c()
	# format it as a character matrix
	if(inherits(x, "data.frame")) {
		dimnames.x <- dimnames(x)
		dim.x <- dim(x)
		x <- sapply(x, format)
		if(dim.x[1] <= 1) {
			x <- matrix(x, nrow = dim.x[1], ncol = dim.x[2], byrow
				 = T)
			dimnames(x) <- dimnames.x
		}
	}
	if(!missing(search.url) && !missing(id)) {
		fff <- function(x, search.url)
		{
			"<a href=" %,% search.url %,% x %,% "> " %,% x
		}
		x[, id.idx] <- sapply(id.vals, FUN = fff, search.url = 
			search.url)
	}
	out <- c(out, paste(collapse = " ", "<TABLE", paste(collapse = " ",
		table.attributes), ">"))
	if(!is.null(main) && main != "")
		out <- c(out, paste("<CAPTION> <H3>", main, "</H3> </CAPTION>",
			collapse = " "))
	out <- c(out, paste("<TR>", paste(paste("<TH", column.attributes, ">"),
		c(" ", dimnames.x[[2]]), "</TH>", collapse = " "), "</TR>",
		collapse = " "))
	# data rows
	row.header.attributes <- rep(row.header.attributes, len = nrow(x))
	for(i in seq(length = nrow(x))) {
		out <- c(out, paste(paste("\n   <TR", row.header.attributes[
			i], ">"), "<TH>", dimnames.x[[1]][i], "</TH>", paste(
			sep = "", paste("\n      <TD", cell.attributes, ">",
			sep = " "), x[i,  ], "</TD>", collapse = " "), "</TR>",
			collapse = " "))
	}
	out <- c(out, "</TABLE>")
	# return vector of strings or write to file and return file name 
	write(out, file = file, append = F)
        if (is.null(browser)) 
          shell.exec("file://" %,% file)
	else
          browseURL("file://" %,% file)
        invisible(file)
}

"means.var" <- 
function(data, labels, subset)
{
	nms <- names(data)
	m <- call. <- match.call()
	m[[1]] <- as.name("model.frame")
	idx <- c(sapply(labels,FUN=grep, nms))
	form <- make.form("", nms[idx])
	m$formula <- form
	m$labels <- NULL
	m <- eval(m, sys.parent())
	id <- dimnames(m)[[1]]
	gnum <- dim(m)[1]
	clnum <- length(labels)
	N <- gnum
	p <- clnum
	n <- length(grep(labels[1], nms))
	mu.g <- matrix(NA, N, p)
	Y <- matrix(NA, N, n * p)
	var.g <- matrix(NA, N, p^2)
	#find within group mean variance and sample size
	all.cols <- NULL
	Terms <- attr(m, "terms")
	Y <- model.matrix(Terms, m)[, -1, drop = F]
	for(i in 1:p) {
		mu.g[, i] <- rowMeans(Y[, n * (i - 1) + (1:n)], na.rm = T)
		for(j in unique(1:i)) {
			Z.i <- Y[, n * (i - 1) + (1:n)] - mu.g[, i]
			Z.j <- Y[, n * (j - 1) + (1:n)] - mu.g[, j]
			var.g[, p * (i - 1) + j] <- var.g[, p * (j - 1) + i] <-
				n/(n - 1) * rowMeans(Z.i * Z.j, na.rm = T)
		}
	}
	list(mean = mu.g, var = var.g)
}

"as.data.frame.fit.n.data" <- 
function(x, ...)
{
  ans <- x$data
  attr(ans, "EBfit") <- x$EBfit
  ans
}

"print.EBfit" <- 
function(x, ...){
  Var.Struct <- x$call$Var.Struct
  th <- x$coefficients
  V.th <- x$variance
  s <- x$shape
  if(Var.Struct!="simple"){
    d <- dim(x$rate)[1]
    J <- c((2*d+1)*exp(th[1]), exp(th[2:(d+1)]))
    upper <- outer(1:d, 1:d, FUN="<=")
    r <- x$rate[upper]
    if(d>1) J <- c(J, rep(1, d*(d-1)/2))
    ij <- outer(1:d, 1:d, FUN="%,%")[upper]
    r.nms <- "rate" %,% ij      
    J <- diag(J)
    V.s.r <- J %*% V.th %*% J
    se <- diag(V.s.r)^0.5
  }
  if(Var.Struct=="simple"){
    r <- x$rate
    J <- diag(exp(th))
    V.s.r <- J %*% V.th %*% J
    se <- diag(V.s.r)^0.5
    r.nms <- "rate"
  }
  coefs <- c(s, r)
  tbl <- cbind(coefs, se, 2*(1-pnorm(abs(coefs/se))))
  dimnames(tbl) <- list(c("shape", r.nms), c("est.","std.err","p-val"))
  print(x$call)
  print(tbl)
  invisible(x)
}

"print.fit.n.data" <- 
function(x, ...)
{
  top <- TopGenes(x, ...)
  n.x <- length(x)
  print(top)
  Var.Struct <- x$EBfit$call$Var.Struct
  m.type <- grep(Var.Struct, c("simple", "general"))
  m.nm <- c("Multivariate Normal/Inverse Wishart", "Normal/Inverse Gamma")[m.type]
  fit <- x$EBfit  
  cat("\n\n" %,% m.nm %,% " model fit: \n")
  print(fit)
  invisible(x)
}

"summary.fit.n.data" <-
function(object, ...)
{
  top <- TopGenes(object, ...)
  n.object <- length(object)
  n.top <- dim(top)[1]
  print(top)
  object[[n.object+1]] <- top
  names(object)[n.object+1] <- "Top" %,% n.top
  Var.Struct <- object$EBfit$call$Var.Struct
  m.type <- grep(Var.Struct, c("general", "simple"))
  m.nm <- c("Wishart/Inverse Wishart", "Chi-Squared/Inverse Gamma")[m.type]
  fit <- object$EBfit
  cat("\n\n" %,% m.nm %,% " model fit: \n")
  print(fit)
  invisible(object)
}

"EBfit" <-
function(object)
{
  fit <- object$EBfit
  fit
}

"update.fit.n.data" <-
function(object, ...)
{
  fit <- object$EBfit
  update(fit, ...)
}

"rowInvc" <- 
function(x)
{
	x <- t(x)
	dx <- dim(x)
	if(floor(dx[1]^0.5) != dx[1]^0.5)
		stop("x must be square")
	n <- dx[2]
	nx <- dx[1]^0.5
	ans <- .C("rowinv",
		x= as.double(x),
		z= as.double(rep(0, n * nx^2)),
		n= as.integer(n),
		nx= as.integer(nx),
                PACKAGE="SharedHT2")[[2]]
	t(matrix(ans, nx^2, n))
}

"rowCholc" <- 
function(x)
{
	x <- t(x)
	dx <- dim(x)
	if(floor(dx[1]^0.5) != dx[1]^0.5)
		stop("x must be square")
	n <- dx[2]
	nx <- dx[1]^0.5
	ans <- .C("rowchol",
		x = as.double(x),
		z = as.double(rep(0, n * nx^2)),
		n = as.integer(n),  
		nx = as.integer(nx),
                PACKAGE="SharedHT2")[[2]]
	t(matrix(ans, nx^2, n))
}

"rowProdc" <- 
function(x, y, na1, nb1)
{
	dx <- dim(x)
	dy <- dim(y)
	if(dx[1] != dy[1])
		stop("x and y must have same #rows")
	n <- dx[1]
	nb2 <- dy[2]/nb1
	na2 <- dx[2]/na1
	if(na2 != nb1)
		stop("Implied Matrices Not of Compatable Dimensions")
	ans <- .C("rowprod",
		x= as.double(x),
		y= as.double(y),
		z= as.double(rep(0, n * na1 * nb2)),
		n= as.integer(n),
		na1= as.integer(na1),
		nb1= as.integer(nb1),
		nb2= as.integer(nb2),
                PACKAGE="SharedHT2")[[3]]
	matrix(ans, n, na1 * nb2)
}

"SimW.IW" <- 
function(nsim, nu=NULL, Lambda=NULL, theta=NULL, nreps, Ngenes, effect.size,
         FDRlist = 0.05*(1:5), verbose=F, gradient=F)
{
	.call. <- match.call()
	is.nuL <- !(missing(nu)||missing(Lambda))
	is.theta <- !missing(theta)
	if(!is.nuL&&!is.theta) 
	   stop("Must provide _either_ (i) both 'nu' and 'Lambda' _or_ (ii) 'theta'")
	if(is.theta) {
		nuL <- shape.rate(theta)
		nu <- nuL$shape
		Lambda <- nuL$rate
	}
	d.L <- dim(Lambda)
	issym <- sum(abs(t(Lambda) - Lambda)==0)
	if(d.L[1]!=d.L[2]||!issym) stop("'Lambda' must be a square symmetric matrix")
	Lbdinvhlf <- chol(solve(Lambda))
	d <- dim(Lambda)[1]
	d2 <- d^2
        npar <- d*(d+1)/2 + 1
	if(length(nreps) == Ngenes) sumnreps <- sum(nreps)
	else {
		sumnreps <- Ngenes * nreps[1]
		nreps <- rep(nreps[1], Ngenes)
	}
        
        nFDRlist <- length(FDRlist)

        ans <- .C("SimW_IW",
                  verb = as.integer(verbose),
                  fail = as.integer(0),
                  fncnt = as.integer(0),
                  grcnt = as.integer(0),
                  mask = as.integer(0),
                  usegr = as.integer(gradient),
                  pnsim = as.integer(nsim),
                  nu = as.double(nu),
                  Lbdinvhlf = as.double(Lbdinvhlf),
                  pd = as.integer(d),
                  pnreps = as.integer(nreps),
                  pN = as.integer(Ngenes),
                  es = as.double(effect.size),
                  coef = as.double(rep(0,npar*nsim)),
                  coefEV = as.double(rep(0,2*nsim)),
                  FDRlist = as.double(FDRlist),
                  pnFDRlist = as.integer(nFDRlist),
                  fdrtbl = as.double(rep(0,8*nFDRlist)),
                  roctbl = as.double(rep(0, 8*Ngenes)),
                  PACKAGE = "SharedHT2")
        
        coef <- t(matrix(ans$coef,npar,nsim))
        dimnames(coef) <- list(1:nsim, "th" %,% (1:npar))

        coefEV <- t(matrix(ans$coefEV,2,nsim))
        dimnames(coefEV) <- list(1:nsim, c("s","r"))

        fdrtbl <- matrix(ans$fdrtbl, nFDRlist, 8)/nsim
        dimnames(fdrtbl) <- list(signif(FDRlist, 3),
                              c(t(outer(c("ShHT2","HT2","ShUT2","UT2"),
                                        c("-TPR","-FPR"), FUN="%,%"))))
        roctbl <- matrix(ans$roctbl/nsim, Ngenes, 8)
        dimnames(roctbl) <- list(rep("",Ngenes),
                             c(t(outer(c("ShHT2","HT2","ShUT2","UT2"),
                                       c("-TPR","-FPR"), FUN="%,%"))))       
        
	list(fdrtbl=fdrtbl, roctbl=roctbl, coef=coef,coefEV=coefEV, call = .call.)
}

"Simnu.mix" <- 
function(nsim, nu=NULL, Lambda=NULL, theta=NULL, nreps, Ngenes, effect.size,
         FDRlist = 0.05*(1:5), f1f2 = c(1/4, 1/2), verbose=F, gradient=F)
{
	.call. <- match.call()
	is.nuL <- !(missing(nu)||missing(Lambda))
	is.theta <- !missing(theta)
	if(!is.nuL&&!is.theta) 
	   stop("Must provide _either_ (i) both 'nu' and 'Lambda' _or_ (ii) 'theta'")
        if(f1f2[1]/f1f2[2] - 1 >=  -1e-10)
           stop("Ratio f1/f2 must be less than 1.  Try something else !!!")
	if(is.theta) {
		nuL <- shape.rate(theta)
		nu <- nuL$shape
		Lambda <- nuL$rate
	}
	d.L <- dim(Lambda)
	issym <- sum(abs(t(Lambda) - Lambda)==0)
	if(d.L[1]!=d.L[2]||!issym) stop("'Lambda' must be a square symmetric matrix")
	Lbdinvhlf <- chol(solve(Lambda))
	d <- dim(Lambda)[1]
	d2 <- d^2
        npar <- d*(d+1)/2 + 1
	if(length(nreps) == Ngenes) sumnreps <- sum(nreps)
	else {
		sumnreps <- Ngenes * nreps[1]
		nreps <- rep(nreps[1], Ngenes)
	}
        
        nFDRlist <- length(FDRlist)

        ans <- .C("Simnu_mix",
                  verb      = as.integer(verbose),
                  fail      = as.integer(0),
                  fncnt     = as.integer(0),
                  grcnt     = as.integer(0),
                  mask      = as.integer(0),
                  usegr     = as.integer(gradient),
                  pnsim     = as.integer(nsim),
                  nu        = as.double(nu),
                  Lbdinvhlf = as.double(Lbdinvhlf),
                  pd        = as.integer(d),
                  pnreps    = as.integer(nreps),
                  pN        = as.integer(Ngenes),
                  es        = as.double(effect.size),
                  f1f2      = as.double(f1f2),
                  coef      = as.double(rep(0,npar*nsim)),
                  coefEV    = as.double(rep(0,2*nsim)),
                  FDRlist   = as.double(FDRlist),
                  pnFDRlist = as.integer(nFDRlist),
                  fdrtbl    = as.double(rep(0,8*nFDRlist)),
                  roctbl = as.double(rep(0,8*Ngenes)),
                  PACKAGE   = "SharedHT2")
        
        coef <- t(matrix(ans$coef,npar,nsim))
        dimnames(coef) <- list(1:nsim, "th" %,% (1:npar))

        coefEV <- t(matrix(ans$coefEV,2,nsim))
        dimnames(coefEV) <- list(1:nsim, c("s","r"))

        fdrtbl <- matrix(ans$fdrtbl, nFDRlist, 8)/nsim
        dimnames(fdrtbl) <- list(signif(FDRlist, 3),
                              c(t(outer(c("ShHT2","HT2","ShUT2","UT2"),
                                        c("-TPR","-FPR"), FUN="%,%"))))

        roctbl <- matrix(ans$roctbl/nsim, Ngenes, 8)
        dimnames(roctbl) <- list(rep("",Ngenes),
                              c(t(outer(c("ShHT2","HT2","ShUT2","UT2"),
                                        c("-TPR","-FPR"), FUN="%,%"))))

	list(fdrtbl=fdrtbl, roctbl=roctbl, coef=coef,coefEV=coefEV, call = .call.)
}

"SimOneW.IW" <- 
function(nu=NULL, Lambda=NULL, theta=NULL, nreps, Ngenes, effect.size)
{
  .call. <- match.call()
  is.nuL <- !(missing(nu)||missing(Lambda))
  is.theta <- !missing(theta)
  if(!is.nuL&&!is.theta)
    stop("Must provide _either_ (i) both 'nu' and 'Lambda' _or_ (ii) 'theta'")
  if(is.theta) {
    nuL <- shape.rate(theta)
    nu <- nuL$shape
    Lambda <- nuL$rate
  }
  d.L <- dim(Lambda)
  issym <- sum(abs(t(Lambda) - Lambda)==0)
  if(d.L[1]!=d.L[2]||!issym) stop("'Lambda' must be a square symmetric matrix")
  Lbdinvhlf <- chol(solve(Lambda))
  d <- dim(Lambda)[1]
  d2 <- d^2

  ans <- .C("SimOneW_IW",
          nu = as.double(nu),
          Lbdinvhlf = as.double(Lbdinvhlf),
          pd = as.integer(d),
          pnreps = as.integer(nreps),
          pN = as.integer(Ngenes),
          es = as.double(effect.size),
          YY = as.double(rep(0,Ngenes*nreps*d)),
          PACKAGE="SharedHT2")
  nms <- c(outer("log2.grp" %,% (1:d), ".rep" %,% (1:nreps), FUN="%,%"))
  DAT <- data.frame(id=1:Ngenes)
  DAT[,nms] <- t(matrix(ans$YY,nreps*d, Ngenes))
  DAT
}

"SimOnenu.mix" <- 
function(nu=NULL, Lambda=NULL, theta=NULL, f1f2=c(1/4, 1/2), nreps, 
         Ngenes, effect.size)
{
  .call. <- match.call()
  is.nuL <- !(missing(nu)||missing(Lambda))
  is.theta <- !missing(theta)
  if(!is.nuL&&!is.theta)
    stop("Must provide _either_ (i) both 'nu' and 'Lambda' _or_ (ii) 'theta'")
  if(is.theta) {
    nuL <- shape.rate(theta)
    nu <- nuL$shape
    Lambda <- nuL$rate
  }
  d.L <- dim(Lambda)
  issym <- sum(abs(t(Lambda) - Lambda)==0)
  if(d.L[1]!=d.L[2]||!issym) stop("'Lambda' must be a square symmetric matrix")
  Lbdinvhlf <- chol(solve(Lambda))
  d <- dim(Lambda)[1]
  d2 <- d^2

  ans <- .C("SimOnenu_mix",
          nu = as.double(nu),
          Lbdinvhlf = as.double(Lbdinvhlf),
          f1f2 = as.double(f1f2),
          pd = as.integer(d),
          pnreps = as.integer(nreps),
          pN = as.integer(Ngenes),
          es = as.double(effect.size),
          YY = as.double(rep(0,Ngenes*nreps*d)),
          PACKAGE="SharedHT2")
  nms <- c(outer("log2.grp" %,% (1:d), ".rep" %,% (1:nreps), FUN="%,%"))
  DAT <- data.frame(id=1:Ngenes)
  DAT[,nms] <- t(matrix(ans$YY,nreps*d, Ngenes))
  DAT
}

"find.ncp" <- 
function(e.I, e.II, nreps, d)
{
  "%,%" <- function(x,y)paste(x,y,sep="")
  obj <- function(logtheta, nreps, d)
         {
           df1 <- nreps*d
           df2 <- d*(nreps-1)
           theta <- exp(logtheta)
           ncp <- nreps*d*exp(2*logtheta)
           power <- 1-pf(qf(1-e.I, df1, df2), df1, df2, ncp)
           cat(sprintf("power=%f, theta=%f, ncp=%f\n",power,theta,ncp))
           ans <- (1-e.II - power)
           attr(theta, "power") <- power
           ans
         }
  ans <- uniroot(f=obj, interval=c(-5,5), nreps=nreps, d=d)
  df1 <- nreps*d
  df2 <- d*(nreps-1)
  logtheta <- ans$root
  theta <- exp(logtheta)
  ncp <- nreps*d*exp(2*logtheta)
  power <- 1-pf(qf(1-e.I, df1, df2), df1, df2, ncp)
  c(theta=theta, power=power)
}

"trgamma" <- 
function(n, shape, scale)
{
	ans <- .C("rgammans",
	pn = as.integer(n),
	shape = as.double(shape),
	scale = as.double(scale),
	ans = as.double(rep(0,n)),
        PACKAGE="SharedHT2")
	ans$ans
}

"trnorm" <- 
function(n)
{
	ans <- .C("rnormns",
	pn = as.integer(n),
	ans = as.double(rep(0,n)),
        PACKAGE="SharedHT2")
	ans$ans
}

"rwishart" <- 
function(n, df, Lambda)
{
	Lambdahlf <- chol(Lambda)
	d <- dim(Lambdahlf)[1]
	ans <- .C("rwishartns",
	pn = as.integer(n),
	pdf = as.double(df),
	pd = as.integer(d),
	pSqrtSigma = as.double(Lambdahlf),
	prows = as.integer(1),
	pW = as.double(rep(0, n*d^2)),
	PACKAGE="SharedHT2")
	matrix(ans$pW, n, d^2)
}

"%,%" <- 
function (x, y)paste(x, y, sep = "")

"if.bad.0" <- 
function (x)
{
    ans <- x
    if (is.vector(x)) {
        ans[is.na(x) | (abs(x) == Inf)] <- 0
    }
    if (is.array(x)) {
        ans <- c(ans)
        ans[is.na(x) | (abs(x) == Inf)] <- 0
        ans <- array(ans, dim(x))
    }
    ans
}

"%K%" <-
function (x, y)
{
    kronecker(x, y)
}

"make.form" <- 
function (y, x)
{
    string <- paste(y, " ~ ", sep = "")
    p <- length(x)
    string <- paste(string, x[1], sep = "")
    if (p > 1)
        for (k in 2:p) string <- paste(string, paste(" + ", x[k],
            sep = ""), sep = "")
    form <- formula(string)
    form
}

ShHT2News <- function() {
    newsfile <- file.path(system.file(package="SharedHT2"), "NEWS")
    file.show(newsfile)
}

.onAttach <- function(libname, pkgname) {
    ShHT2ver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                      fields="Version")
    cat(paste(pkgname, ShHT2ver, "\n"))
    cat("Type ShHT2News() to see new features/changes/bug fixes.\n")
}

# Script to plot marginal ancestral probabilities.
# 
#
# Ignacio Quintero
#  t(-_-t)
# 
# 30 11 20
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-


library(ggtree)
library(data.table)


hmean = function(x) 1.0/mean(1.0/x)


linpred = function(t, t1, t2, y1, y2) {
  (y1 + (t - t1)*(y2 - y1)/(t2 - t1))
}


# prepare lineages to start simulations at time T
prepare_lineage_fs = function(tree_file, states_file, pars_file, T) {

  tr = read.tree(tree_file)
  nf = fread(states_file, header = TRUE)

  ## summarize ancestral states
  # get node names 
  cnams = colnames(nf)
  cnams = strsplit(cnams, '_')[3:length(cnams)]
  nd = as.numeric(sapply(cnams, '[', 2))
  st = as.numeric(sapply(cnams, '[', 4))

  # check if at least 80% is non NaN
  nans = nf[,lapply(.SD, is.nan)]
  rnan = apply(nans, 1, any)

  if (sum(rnan) < 0.5*nrow(nf)) {
    nf = nf[!rnan,]
    if (sum(rnan) > 0) {
      cat("Removed", sum(rnan), "rows with NaN \n")
    }
  } else {
    cat("File with too many NaNs \n")
  }

  ns = length(unique(st))
  nn = length(unique(nd))

  pr = nf[,3:ncol(nf)]
  if (ns == 3) {
    mpr = data.table(Iteration = numeric(0),
                 Posterior = numeric(0),
                 node      = numeric(0),
                 V1        = numeric(0),
                 V2        = numeric(0),
                 V3        = numeric(0))
  }
  if (ns == 6) {
    mpr = data.table(Iteration = numeric(0),
                     Posterior = numeric(0),
                     node      = numeric(0),
                     V1        = numeric(0),
                     V2        = numeric(0),
                     V3        = numeric(0),
                     V4        = numeric(0),
                     V5        = numeric(0),
                     V6        = numeric(0))
  }

  for (i in 1:nrow(pr)) {
    mpri = as.data.frame(matrix(as.numeric(pr[i,]), nn, ns, byrow = TRUE))
    setDT(mpri)
    mpr = rbind(mpr, cbind(nf[i,1], nf[i,2], node = unique(nd), mpri))
  } 

  tb = fortify(tr)
  setDT(tb)

  # add state probabilities
  tb = merge(tb, mpr, by = 'node')

  # make tb have start and end probabilities
  ui = unique(tb[,Iteration])
  for (i in ui) {
    tbi = tb[Iteration == i]
    pr_p = tbi[match(tbi[,parent],tbi[,node]), 12:(11+ns)]
    setnames(pr_p, paste0('Vi', 1:ns))
    tb[Iteration == i, paste0('Vi', 1:ns) := pr_p]
  }

  # add start time
  tb[, xi := x - branch.length]

  # get lineages that intersect time T
  rT = max(tb[,x]) - T
  al = tb[xi <= rT & x >= rT]

  # get average states at time rT
  al[, Vp1 := linpred(rT, xi, x, Vi1, V1)]
  al[, Vp2 := linpred(rT, xi, x, Vi2, V2)]
  al[, Vp3 := linpred(rT, xi, x, Vi3, V3)]
  if (ns == 6) {
    al[, Vp4 := linpred(rT, xi, x, Vi4, V4)]
    al[, Vp5 := linpred(rT, xi, x, Vi5, V5)]
    al[, Vp6 := linpred(rT, xi, x, Vi6, V6)]
  }

  # match with parameters
  pf = fread(pars_file)

  pf = pf[Iteration %in% ui]
  pf[,Posterior := NULL]
  al = merge(al, pf, by = 'Iteration')

  return(al)
}


# to cartesian coordinates
xf = function(a, r, j = 0.) r*cos((a*pi/180)) + j
yf = function(a, r, k = 0.) r*sin((a*pi/180)) + k


# estimate probabilities and prepare tree to plot
prepare_nmsp = function(tree_file, states_file, pars_file, n_slices = 200) {

  tr = read.tree(tree_file)
  nf = fread(states_file, header = TRUE)

  # get node names 
  cnams = colnames(nf)
  cnams = strsplit(cnams, '_')[3:length(cnams)]
  nd = as.numeric(sapply(cnams, '[', 2))
  st = as.numeric(sapply(cnams, '[', 4))

  # check if at least 80% is non NaN
  nans = nf[,lapply(.SD, is.nan)]
  rnan = apply(nans, 1, any)

  if (sum(rnan) < 0.5*nrow(nf)) {
    nf = nf[!rnan,]
    if (sum(rnan) > 0) {
      cat("Removed", sum(rnan), "rows with NaN \n")
    }
  } else {
    cat("File with too many NaNs \n")
  }

  # estimate marginal probabilities
  ss = sum(nf[,2])
  mp = nf[ , lapply(.SD, "*", nf[,2])]
  mp = apply(mp, 2, sum, na.rm = TRUE)
  pr = mp[3:length(mp)]/ss

  ns = length(unique(st))
  nn = length(unique(nd))

  mpr = as.data.frame(matrix(pr, nn, ns, byrow = TRUE))
  mpr = cbind(node = unique(nd), mpr)
  setDT(mpr)

  # make plotting table from ggplot
  tb = fortify(tr)
  setDT(tb)

  # add state probabilities
  tb = merge(tb, mpr, by = 'node')

  # make tb have start and end probabilities
  pr_p = tb[match(tb[,parent],tb[,node]), 10:(9+ns)]
  setnames(pr_p, paste0('Vi', 1:ns))
  tb[, paste0('Vi', 1:ns) := pr_p]

  # make circular coordinates
  tb = tb[order(parent, angle)]

  tb[, xi := x - branch.length]

  tb[, `:=`(xci = xf(angle, xi),
            xcf = xf(angle, x),
            yci = yf(angle, xi),
            ycf = yf(angle, x)
            )]

  # make circular segments joining straight lines
  ypair = tb[,.(xi, angle, node), parent]
  dups  = which(duplicated(ypair, by = 'parent'))

  ypair = ypair[order(parent, angle)]

  # initial and end angle
  ai = ypair[dups-1,angle]
  af = ypair[dups,  angle]
  xi = ypair[dups-1, xi]

  # make x and y circle coordinates
  pairl = list()
  for (i in 1:length(ai)) {
    yseq = seq(ai[i], af[i], 0.01)
    xc = xf(yseq, xi[i])
    yc = yf(yseq, xi[i])
    pairl[[i]] = list(x = xc, y = yc)
  }

  # prepare pars densities
  prf = fread(pars_file, header = TRUE)

  sdat = state_slice(tb, ns, nn, by = 1.0)

  return(list(tb = tb, prf = prf, sdat = sdat, pairl = pairl, ns = ns, nn = nn))
}



# plot node marginal state probabilities
plot_nmsp = function(d,
  coltree   = 'black', 
  colstates = c('purple', 'orange', 'grey', 'blue', 'red', 'black'),
  w         = 2,
  tiplabels = FALSE,
  tlper     = 0.03,
  wcircle   = 0.15,
  lrect     = 0.2) {

  tb    = .subset2(d, 'tb')
  pairl = .subset2(d, 'pairl')
  ns    = .subset2(d, 'ns')
  nn    = .subset2(d, 'nn')
  nt = which(tb$isTip == TRUE)
  ni = which(tb$isTip == FALSE)

  # reorder to have hidden states of observed states together
  if (ns == 6) {
    tb = tb[, c(1:9,10,13,11,14,12,15,(10+ns):ncol(tb)), with = FALSE]
  }

  r  = max(tb[,x])

  # make pie graphs centered at 0 for internal nodes
  pies = list()
  for (j in 1:length(ni)) {

    xx = as.numeric(tb[ni[j],10:(9+ns)])
    cs = cumsum(xx)
    ax = 360*cs/cs[length(cs)]
    ax = c(0.0, ax)

    pslice = list()
    for (i in 1:ns) {
      axi = ax[i]
      axf = ax[i+1]
      if (axi == axf) {
        pslice[[i]] = list(x = NULL, y = NULL)
        next
      }
      yseq = seq(axi, axf, 0.5)
      xc = xf(yseq, wcircle*r)
      yc = yf(yseq, wcircle*r)
      pslice[[i]] = list(x = c(0., xc, 0.), y = c(0., yc, 0.))
    }

    pies[[j]] = pslice
  }

  # determine rectangle width
  rwd = 360/length(nt)/2

  rects = list()
  # make rectangles centered for tips
  for (j in 1:length(nt)) {
    tbi = tb[nt[j],]
    xx  = as.numeric(tbi[,10:(9+ns)])
    ang = tbi[1,angle]
    cs  = cumsum(xx)
    ax  = c(0.0, cs)

    prect = list()
    for (i in 1:ns) {
      axi = ax[i]
      axf = ax[i+1]
      if (axi == axf) {
        prect[[i]] = list(x = NULL, y = NULL)
        next
      }
      # make semi-circles
      yseq = seq(ang-rwd, ang+rwd, 0.01)
      xc0 = xf(yseq,r+lrect*r*axi)
      yc0 = yf(yseq,r+lrect*r*axi)
      xc1 = xf(yseq,r+lrect*r*axf)
      yc1 = yf(yseq,r+lrect*r*axf)
      prect[[i]] = 
        list(x = c(xc0, rev(xc1), xc0[1]), 
             y = c(yc0, rev(yc1), yc0[1]))
    }
    rects[[j]] = prect
  }


  # make tree plot
  plot(1, xlim = c(-r-w,r+w), ylim = c(-r-w,r+w), type = 'n', bty = 'n',
    xaxt = 'n', yaxt = 'n', ylab = '', xlab = '')

  # plot tree
  segments(tb[,xci], tb[,yci], tb[,xcf], tb[,ycf], col = coltree)
  for (i in 1:length(pairl)) {
    sl = .subset2(pairl,i)
    lines(sl$x, sl$y, col = coltree)
  }

  # plot tip labels
  if (tiplabels == TRUE) {
    tb[isTip == TRUE, 
      `:=`(tx = xf(angle, x * (1+tlper*r)),
           ty = yf(angle, x * (1+tlper*r)))]
    for (i in which(tb[,isTip])) {
      text(tb[i, tx], tb[i, ty], 
        labels = tb[i, label], 
        srt = tb[i, angle])
    }
  }

  # plot pies
  for (j in 1:length(ni)) {
    psli = pies[[j]]
    for (i in 1:ns) {
        polygon(psli[[i]]$x + tb[ni[j],xcf], psli[[i]]$y + tb[ni[j],ycf], 
          col = colstates[i], border = colstates[i])
    }
  }

  # plot rects
  for (j in 1:length(rects)) {
    recti = rects[[j]]
    for (i in 1:ns) {
        polygon(recti[[i]]$x, recti[[i]]$y, 
          col = colstates[i], border = colstates[i])
    }
  }
}





# plot node marginal state probabilities
plot_nmsp_hs = function(d,
  coltree   = 'black', 
  colstates = c('purple', 'orange', 'grey', 'blue', 'red', 'black'),
  w         = 2,
  tiplabels = FALSE,
  tlper     = 0.03,
  wcircle   = 0.15,
  lrect     = 0.2) {

  tb    = .subset2(d, 'tb')
  pairl = .subset2(d, 'pairl')
  ns    = .subset2(d, 'ns')
  nn    = .subset2(d, 'nn')
  nt = which(tb$isTip == TRUE)
  ni = which(tb$isTip == FALSE)

  # reorder to have hidden states of observed states together
  if (ns == 6) {
    tb = tb[, c(1:9,10,13,11,14,12,15,(10+ns):ncol(tb)), with = FALSE]
  }

  # sum over observed states to leave only hidden states
  tb[,V1 := V1 + V2 + V3]
  tb[,V2 := V4 + V5 + V6]
  tb[,V3 := NULL]
  tb[,V4 := NULL]
  tb[,V5 := NULL]
  tb[,V6 := NULL]

  tb[,Vi1 := Vi1 + Vi2 + Vi3]
  tb[,Vi2 := Vi4 + Vi5 + Vi6]
  tb[,Vi3 := NULL]
  tb[,Vi4 := NULL]
  tb[,Vi5 := NULL]
  tb[,Vi6 := NULL]

  r  = max(tb[,x])

  # make pie graphs centered at 0 for internal nodes
  pies = list()
  for (j in 1:length(ni)) {

    xx = as.numeric(tb[ni[j],10:11])
    cs = cumsum(xx)
    ax = 360*cs/cs[length(cs)]
    ax = c(0.0, ax)

    pslice = list()
    for (i in 1:2) {
      axi = ax[i]
      axf = ax[i+1]
      if (axi == axf) {
        pslice[[i]] = list(x = NULL, y = NULL)
        next
      }
      yseq = seq(axi, axf, 0.5)
      xc = xf(yseq, wcircle*r)
      yc = yf(yseq, wcircle*r)
      pslice[[i]] = list(x = c(0., xc, 0.), y = c(0., yc, 0.))
    }

    pies[[j]] = pslice
  }

  # determine rectangle width
  rwd = 360/length(nt)/2

  rects = list()
  # make rectangles centered for tips
  for (j in 1:length(nt)) {
    tbi = tb[nt[j],]
    xx  = as.numeric(tbi[,10:11])
    ang = tbi[1,angle]
    cs  = cumsum(xx)
    ax  = c(0.0, cs)

    prect = list()
    for (i in 1:2) {
      axi = ax[i]
      axf = ax[i+1]
      if (axi == axf) {
        prect[[i]] = list(x = NULL, y = NULL)
        next
      }
      # make semi-circles
      yseq = seq(ang-rwd, ang+rwd, 0.01)
      xc0 = xf(yseq,r+lrect*r*axi)
      yc0 = yf(yseq,r+lrect*r*axi)
      xc1 = xf(yseq,r+lrect*r*axf)
      yc1 = yf(yseq,r+lrect*r*axf)
      prect[[i]] = 
        list(x = c(xc0, rev(xc1), xc0[1]), 
             y = c(yc0, rev(yc1), yc0[1]))
    }
    rects[[j]] = prect
  }


  # make tree plot
  plot(1, xlim = c(-r-w,r+w), ylim = c(-r-w,r+w), type = 'n', bty = 'n',
    xaxt = 'n', yaxt = 'n', ylab = '', xlab = '')

  # plot tree
  segments(tb[,xci], tb[,yci], tb[,xcf], tb[,ycf], col = coltree)
  for (i in 1:length(pairl)) {
    sl = .subset2(pairl,i)
    lines(sl$x, sl$y, col = coltree)
  }

  # plot tip labels
  if (tiplabels == TRUE) {
    tb[isTip == TRUE, 
      `:=`(tx = xf(angle, x * (1+tlper*r)),
           ty = yf(angle, x * (1+tlper*r)))]
    for (i in which(tb[,isTip])) {
      text(tb[i, tx], tb[i, ty], 
        labels = tb[i, label], 
        srt = tb[i, angle])
    }
  }

  # plot pies
  for (j in 1:length(ni)) {
    psli = pies[[j]]
    for (i in 1:2) {
        polygon(psli[[i]]$x + tb[ni[j],xcf], psli[[i]]$y + tb[ni[j],ycf], 
          col = colstates[i], border = colstates[i])
    }
  }

  # plot rects
  for (j in 1:length(rects)) {
    recti = rects[[j]]
    for (i in 1:2) {
        polygon(recti[[i]]$x, recti[[i]]$y, 
          col = colstates[i], border = colstates[i])
    }
  }
}





# linear interpolation function
linpred = function(t, t0, t1, x0, x1) {
  (x0 + (t - t0)*(x1 - x0)/(t1 - t0))
}



# through time preparation
state_slice = function(tb, ns, nn, by = 1.0) {

  slices  = seq(0, max(tb[,x]), by = by)
  n_slices = length(slices)

  tb[, xirev := abs(xi - max(x))]
  tb[, xrev := abs(x - max(x))]

  # preallocate data.table
  sdat = data.table(slice = slices)

  # make empty columns
  for (k in 1:ns) {
    sdat[, paste0('V',k) := NA_real_]
  }

  # which are start prob and which are end
  pri = grep('Vi', colnames(tb))
  prf = grep('V[0-9]', colnames(tb))

  # Pr at end of branches at sli
  for (j in 1:nrow(sdat)) {
    sli = slices[j]

    bi = tb[xirev >= sli & sli >= xrev,]

    if (sli == 0.0) {
      lp = tb[isTip == TRUE,prf,with=FALSE]
        sdat[slice == sli, (2:(1+ns)) := lp[,lapply(.SD,sum)]]
    } else { 
      lp = linpred(sli, bi[,xirev], bi[, xrev], 
        bi[,pri,with=FALSE], bi[,prf,with=FALSE])
      sdat[slice == sli, (2:(1+ns)) := lp[,lapply(.SD,sum)]]
    }
  }

  return(sdat)
}




plot_nmsp_prop = function(d,
    colstates = c('purple', 'orange', 'grey', 'blue', 'red', 'black')) {

  sdat  = .subset2(d, 'sdat')
  ns    = .subset2(d, 'ns')

  # reorder to have hidden states of observed states together
  if (ns == 6) {
    sdat = sdat[, c(1,2,5,3,6,4,7), with = FALSE]
  }

  # cumsum
  csum = t(apply(sdat[,2:(1+ns), with = FALSE], 1, cumsum))

  # standardized curve
  cssd = csum/csum[,ns]

  # make polygons
  pols = list()
  for (i in 1:ns) {
    x = c(sdat[,slice], rev(sdat[,slice]))
    
    if (i == 1) {
      y = c(cssd[,1], rep(0.0, nrow(cssd)))
    } else {
      if (i == ns) {
        y = c(cssd[,i-1], rep(1.0, nrow(cssd)))
      } else {
        y = c(cssd[,i],rev(cssd[,i-1]))
      }
    }
    pols[[i]] = list(x = x, y = y)
  }

  mxx = max(sdat[,slice])
  xax = seq(0, mxx, 30)

  plot(1, type = 'n', xlim = c(-mxx, 0), 
    ylim = c(0,1), 
    bty = 'n', las = 1, ylab = 'Proportion of lineages in each state',
    xlab = 'Time (Mya)', xaxt = 'n')
  axis(1, at = -xax, labels = xax)

  nr = nrow(sdat)

  for (i in 1:ns) {
    p = pols[[i]]
    polygon(-p$x, p$y, col = colstates[i], border = colstates[i])
  }

}





plot_nmsp_ltt = function(d,
    colstates = c('purple', 'orange', 'grey', 'blue', 'red', 'black'),
    ...) {

  sdat  = .subset2(d, 'sdat')
  ns    = .subset2(d, 'ns')

  # reorder to have hidden states of observed states together
  if (ns == 6) {
    sdat = sdat[, c(1,2,5,3,6,4,7), with = FALSE]
  }

  mxx = max(sdat[,slice])
  xax = seq(0, mxx, 30)

  plot(1, type = 'n', xlim = c(-mxx, 0), 
    ylim = c(1, max(sdat[,2:ncol(sdat), with = FALSE])), 
    bty = 'n', las = 1, ylab = 'ln(lineages) in each state',
    xlab = 'Time(Mya)', xaxt = 'n', log = 'y')
  axis(1, at = -xax, labels = xax)

  cnam = colnames(sdat)
  for (i in 1:ns) {
    lines(-sdat[,slice], sdat[,get(cnam[i+1])], 
      col = colstates[i], ...)
  }
}




plot_densities_ns3 = function(d, 
  colstates = c('purple', 'orange', 'grey'),
  logax = FALSE,
  viowd = 0.3) {

  prf = .subset2(d, 'prf')

  # separate beta and other pars
  #organize columns 
  cs = setdiff(3L:ncol(prf), 10:13)

  # make density pols
  dpp = lapply(prf[,cs,with=FALSE], vpol, viowd, log = logax)

  # x axis labels
  pcn = c(expression(lambda[T]),
          expression(lambda[nT]),
          expression(lambda[W]),
          expression(mu[T]),
          expression(mu[nT]),
          'T_nT',
          'nT_T')

  if (ncol(prf) == 13) {
    bs = ncol(prf) - 1:0
    dpb = lapply(prf[,bs,with=FALSE], vpol, viowd, log = logax)
    pcn = c(pcn, c(expression(beta[T]),
                   expression(beta[nT])))
  }

  if (logax) {
    plot(1, type = 'n', xlim = c(0, length(cs) + 3), 
      ylim = c(-10, max(log(prf[,cs,with = FALSE]))), 
      bty = 'n', las = 1, ylab = expression(paste('Pr(',theta,')')),
      xlab = '', xaxt = 'n')
    axis(1, at = 1:(length(cs)), labels = pcn)
  } else {
    plot(1, type = 'n', xlim = c(0, length(cs) + 3), 
      ylim = c(0, max(prf[,cs,with = FALSE])), 
      bty = 'n', las = 1, ylab = expression(paste('Pr(',theta,')')),
      xlab = '', xaxt = 'n')
    axis(1, at = 1:(length(pcn)), labels = pcn)
  }

  # non beta pars
  colo = colstates[c(1,2,3,1,2,1,2)]
  for (i in 1:length(dpp)) {
    dppi = dpp[[i]]
    polygon(.subset2(dppi,'x') + i, .subset2(dppi,'y'), 
      col = colo[i],
      border = colo[i])
  }

  if (ncol(prf) == 13) {
    par(new = TRUE)
    plot(1, type = 'n', xlim = c(0, length(cs) + 3), 
      ylim = range(prf[,bs,with = FALSE]), 
      bty = 'n', las = 1, axes = FALSE, xlab = "", ylab = "")

    # beta pars
    colo = colstates[c(1,2)]
    for (i in 1:length(dpb)) {
      dpbi = dpb[[i]]
      polygon(.subset2(dpbi,'x') + length(cs) + i, .subset2(dpbi,'y'), 
        col = colo[i],
        border = colo[i])
    }
    axis(side=4, at = pretty(range(prf[,bs,with = FALSE])), line = -1.5)
    segments(length(cs)+1, 0.0, length(cs)+3, 0.0, lty = 3)
  }

}




plot_densities_ns6 = function(d, 
  colstates = c('purple', 'orange', 'grey', 'blue', 'red', 'black'),
  logax = FALSE,
  viowd = 0.3) {

  prf = .subset2(d, 'prf')

  #organize columns 
  cs = setdiff(3:ncol(prf),17:20)
  cswob = setdiff(cs, 23:26)

  # make density pols
  dps = lapply(prf[,cs,with=FALSE], vpol, viowd, log = logax)

  # x axis labels
  pcn = c(expression(lambda[T]),
          expression(lambda[nT]),
          expression(lambda[W]),
          expression(mu[T]),
          expression(mu[nT]),
          'T_nT',
          'nT_T',
          '0_1',
          '1_0')

  if (ncol(prf) == 26) {
    bs = ncol(prf) - 3:0
    dpb = lapply(prf[,bs,with=FALSE], vpol, viowd, log = logax)
    pcn = c(pcn, c(expression(beta[T]),
                   expression(beta[nT])))
  }

  if (logax == TRUE) {
    plot(1, type = 'n', xlim = c(0, length(cs)/2 + 2), 
      ylim = c(-10, max(log(prf[,cswob,with = FALSE]))), 
      bty = 'n', las = 1, ylab = expression(paste('Pr(',theta,')')),
      xlab = '', xaxt = 'n')
    axis(1, at = 1:(length(cs)/2 +1), labels = pcn)
  } else {
    plot(1, type = 'n', xlim = c(0, length(cs)/2 + 2), 
      ylim = c(0, max(prf[,cswob,with = FALSE])), 
      bty = 'n', las = 1, ylab = expression(paste('Pr(',theta,')')),
      xlab = '', xaxt = 'n')
    axis(1, at = 1:(length(cs)/2 + 1), labels = pcn)
  }

  col0 = colstates[c(1,3,5,1,3,1,3)]
  col1 = colstates[c(2,4,6,2,4,2,4)]

  i0 = c(1:3,7:8,11:12)
  i1 = c(4:6,9:10,13:14)

  # 0 state
  for (i in 1:length(i0)) {
    dpsi = dps[[i0[i]]]
    polygon(-.subset2(dpsi,'x') + i, .subset2(dpsi,'y'), 
      col = col0[i],
      border = col0[i])
  }

  # 1 state
  for (i in 1:length(i1)) {
    dpsi = dps[[i1[i]]]
    polygon(.subset2(dpsi,'x') + i, .subset2(dpsi,'y'), 
      col = col1[i],
      border = col1[i])
  }

  # hidden transitions
  dpsi = dps[[15]]
  polygon(.subset2(dpsi,'x') + 8, .subset2(dpsi,'y'), 
    col = 'darkgrey',
    border = 'darkgrey')

  dpsi = dps[[16]]
  polygon(.subset2(dpsi,'x') + 9, .subset2(dpsi,'y'), 
    col = 'grey',
    border = 'grey')

  #betas
  if (ncol(prf) == 26) {
    par(new = TRUE)
    plot(1, type = 'n', xlim = c(0, length(cs)/2 + 2), 
      ylim = range(prf[,bs,with = FALSE]), 
      bty = 'n', las = 1, axes = FALSE, xlab = "", ylab = "")

    # beta pars 0
    colo = colstates[c(1,3)]
    for (i in 1:2) {
      dpbi = dpb[[i]]
      polygon(-.subset2(dpbi,'x') + length(cs)/2 - 1 + i, .subset2(dpbi,'y'), 
        col = colo[i],
        border = colo[i])
    }

    # beta pars 1
    colo = colstates[c(2,4)]
    for (i in 3:4) {
      dpbi = dpb[[i]]
      polygon(.subset2(dpbi,'x') + length(cs)/2 - 3 + i, .subset2(dpbi,'y'), 
        col = colo[i-2],
        border = colo[i-2])
    }


    axis(side=4, at = pretty(range(prf[,bs,with = FALSE])), line = -1.5)
    segments(length(cs)/2, 0.0, length(cs)/2+2, 0.0, lty = 3)
  }


}






plot_densities_ns9 = function(d, 
  colstates = c('purple', 'orange', 'grey', 'blue', 'red', 'black'),
  logax = FALSE,
  viowd = 0.3) {

  prf = .subset2(d, 'prf')

  #organize columns 
  cs = setdiff(3:ncol(prf),24:29)

  # make density pols
  dps = lapply(prf[,cs,with=FALSE], vpol, viowd, log = logax)

  # x axis labels
  pcn = c(expression(lambda[T]),
          expression(lambda[nT]),
          expression(lambda[W]),
          expression(mu[T]),
          expression(mu[nT]),
          'T_nT',
          'nT_T',
          '0_1',
          '0_2',
          '1_0',
          '1_0')

  if (logax == TRUE) {
    plot(1, type = 'n', xlim = c(0, length(cs)/2 + 2), 
      ylim = c(-10, max(log(prf[,cs,with = FALSE]))), 
      bty = 'n', las = 1, ylab = expression(paste('Pr(',theta,')')),
      xlab = '', xaxt = 'n')
    axis(1, at = 1:(length(cs)/2 +1), labels = pcn)
  } else {
    plot(1, type = 'n', xlim = c(0, length(cs)/2 + 2), 
      ylim = c(0, max(prf[,cs,with = FALSE])), 
      bty = 'n', las = 1, ylab = expression(paste('Pr(',theta,')')),
      xlab = '', xaxt = 'n')
    axis(1, at = 1:(length(cs)/2 + 1), labels = pcn)
  }

  col0 = colstates[c(1,3,5,1,3,1,3)]
  col1 = colstates[c(2,4,6,2,4,2,4)]

  i0 = c(1:3,7:8,11:12)
  i1 = c(4:6,9:10,13:14)

  # 0 state
  for (i in 1:length(i0)) {
    dpsi = dps[[i0[i]]]
    polygon(-.subset2(dpsi,'x') + i, .subset2(dpsi,'y'), 
      col = col0[i],
      border = col0[i])
  }

  # 1 state
  for (i in 1:length(i1)) {
    dpsi = dps[[i1[i]]]
    polygon(.subset2(dpsi,'x') + i, .subset2(dpsi,'y'), 
      col = col1[i],
      border = col1[i])
  }

  # hidden transitions
  dpsi = dps[[15]]
  polygon(.subset2(dpsi,'x') + 8, .subset2(dpsi,'y'), 
    col = 'darkgrey',
    border = 'darkgrey')

  dpsi = dps[[16]]
  polygon(.subset2(dpsi,'x') + 9, .subset2(dpsi,'y'), 
    col = 'grey',
    border = 'grey')

}







# make a density polygon function
vpol = function(ns, vw = 0.5, log = FALSE) {
  if (log) {
    d  = density(log(ns))
  } else {
    d  = density(ns)
  }
  p = list(x = c(d$y, d$y[1]), y = c(d$x, d$x[1])) 
  p$x = p$x/max(p$x) * vw
  return(p)
}

is.odd = function(x) x %% 2 != 0 






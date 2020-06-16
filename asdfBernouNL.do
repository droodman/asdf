*
* After Benoulli diffusion ML fitting, compute nonlinear derived quantities including s, B, delta, sigma, no-take-off probability, median take-off year
* Uses Stata's nlcom command and *replaces* current estimation results
*

cap program drop asdfBernouNL
program define asdfBernouNL
  scalar s = exp([/lna]) * [/gamma] * ([/nu] + [/gamma])
  local nlcomcmd (lna  : [/lna]) (b:[/b]) (nu:[/nu]) (gamma:[/gamma]) ///
                 (s    : `=cond(e(rank)<e(k), "0", "exp([/lna]) * [/gamma] * ([/nu] + [/gamma]) / s")') ///  // under constrained model, s=0 by fiat and nlcom complains
                 (B    : -1 / [/gamma]) ///
                 (delta: [/b] * [/gamma]) ///
                 (sigma: sqrt(2) * exp([/lna]/2) * abs([/gamma]))
  if ([/nu]+[/gamma]) / [/b] < 0 local nlcomcmd `nlcomcmd' (Y_b  : (-exp([/lna]) * ([/nu]+[/gamma]) / [/b]) ^ [/gamma])   // zero-drift level

  foreach if in i f {  // initial, final
    if [/nu]<=0 & [/gamma]<0 {
      local notakeoffprob gammap(-[/nu], [/b] / exp([/lna]) * Y`if'^(1/[/gamma]))
      scalar notakeoffprob`if' = `notakeoffprob'
      if notakeoffprob`if' local nlcomcmd `nlcomcmd' (notakeoffprob`if': `notakeoffprob' / `=notakeoffprob`if'')  // trick to avoid crash stata.com/statalist/archive/2009-03/msg01244.html
                      else local nlcomcmd `nlcomcmd' (notakeoffprob`if': 0                                     )

      local takeoffb0 Y`if'^(1/[/gamma]) / exp([/lna]) / invgammap(-[/nu], .5)
      if notakeoffprob`if'<.5 {
        local mediantakeoff`if' cond([/b], -ln1m([/b]*`takeoffb0')/[/b], `takeoffb0')
        scalar mediantakeoff`if' = `mediantakeoff`if''
        local nlcomcmd `nlcomcmd' (mediantakeoff`if': `mediantakeoff`if'' / mediantakeoff`if')
      }
    }
    else {
      scalar notakeoffprobi = 1
      scalar notakeoffprobf = 1
    }
  }
  nlcom `nlcomcmd', post

  local nlcomcmd
  foreach var in `:colnames e(b)' {
         if inlist(substr("`var'",1,9), "s", "notakeoff") local nlcomcmd `nlcomcmd' (`var': _b[`var'] * `var'                           )  // finish trick to avoid crash
    else if        substr("`var'",1,13)=="mediantakeoff"  local nlcomcmd `nlcomcmd' (`var': _b[`var'] * `var' + t`=substr("`var'",14,1)')
    else                                                  local nlcomcmd `nlcomcmd' (`var': _b[`var']                                   )
  }
  nlcom `nlcomcmd', post
end

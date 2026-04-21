# Simulation Wrapper

simulate.stochastic <- function(init, transitions, 
                                ratefunc, params, tf = 100, tl.params = list(epsilon = 0.05)) {
  # params must include: t1 (start of intervention), t2 (stop of intervention), N, 
  #beta, gamma, lambda (influx rate), lev (intensity of intervention),  u_func= function(t, params){return(0)}, influx = FALSE
  
  #influx (TRUE /FALSE)
  # define a maxtau, wie finde ich raus, was ein sinnvolles maxtau ist. Was sind die Werte von tau? 
  # if influx = TRUE (default: FALSE) for one infectious import one susceptible less in population 
  # under default parameters, see transitions_influx. 
  
  if(is.null(params$u_func)){
    params$u_func = function(t, params){return(0)}
  }
  if(is.null(params$influx)){
    params$influx = FALSE 
  }
  
  influx <- params$influx
  if (!influx){
  transitions <- rbind(
    
    S = c(-1,  0),
    
    I = c(+1, -1),
    
    R = c( 0, +1)
    
  )
  
  colnames(transitions) <- c("infection", "recovery")
  }else if (influx){
  transitions <- transitions_influx
  
  
  colnames(transitions_influx) <- c("infection", "recovery", "influx")
  } else {
    stop("Influx has to be TRUE or FALSE")
  }
  
  result <- ssa.adaptivetau(
    
    init.values  = init,
    
    transitions  = transitions,
    
    rateFunc     = ratefunc,
    
    params       = params,
    
    tf           = tf,
    
    tl.params    = tl.params
    
  )
  
  # Convert to data.frame with proportions
  
  df <- as.data.frame(result)
  
  df$S_prop <- df$S / params$N
  
  df$I_prop <- df$I / params$N
  
  df$R_prop <- df$R / params$N
  
  
  
  return(df)
  
}


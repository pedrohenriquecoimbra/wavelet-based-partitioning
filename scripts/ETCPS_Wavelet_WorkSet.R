


wavelet <- function(
    # arguments for the wavelet code
  PY_env_path,
  WL_handler_path,
  site_ID,
  output_dir_EP,
  output_dir_WL,
  time_range, 
  raw.freq,
  # arguments for the fetching the needed files  
  RF.cs.dir,                                      # directory containing the RFlux cleanset file
  p.mode,
  BM.dir,
  ST.dir,
  AUX.dir,
  FFP.dir,
  path_ecmd
){
  
  
  # Run the wavelet handler ---------------------------------------------------
  # called by cmd prompt
  WL.cmd <- paste(WL_handler_path, #,
                  site_ID,
                  output_dir_EP,
                  output_dir_WL,
                  time_range,
                  '-sr', raw.freq
  )
  system2(PY_env_path, args = WL.cmd) 
  
  
  # . Import the RFlux WorkSet file ------------------------------------------
  rfcs.0 <- data.frame()
  rfcs.0 <- fread(grep('WorkSet', list.files(RF.cs.dir, full.names = TRUE), value = TRUE), data.table = F, na.strings = '-9999')
  rfcs.0$'TIMESTAMP' <- as.character(rfcs.0$'TIMESTAMP')
  # match TIMESTAMP in WorkSet (START) with Eddypro and Wavelet (END) 
  rfcs.0$'TIMESTAMP' <- format(strptime(rfcs$'TIMESTAMP', '%Y%m%d%H%M', tz='GMT') + 1800, '%Y%m%d%H%M')
  
  # Build a continuous time line starting from the first half-hour of the cleanset
  # This will be used to verify/force the alignment of the output base files
  TIMESTAMP_ref <- data.frame(
    TIMESTAMP = format(seq(strptime(rfcs$'TIMESTAMP'[1], '%Y%m%d%H%M', tz='GMT'), 
                               strptime(tail(rfcs$'TIMESTAMP',1), '%Y%m%d%H%M', tz='GMT'),
                               1800),
                           '%Y%m%d%H%M')
  )
  
  # align the cleanset file
  rfcs <- data.frame()
  rfcs <- base::merge(rfcs.0, TIMESTAMP_ref, by = 'TIMESTAMP', all.y = TRUE)
  # reorder the timestamps
  rfcs <- rfcs[, c('TIMESTAMP',names(rfcs)[-c(1:2)])]
  rm(rfcs.0)
  
  # . Import the EddyPro full_output file ------------------------------------------
  epfo.0 <- data.frame()
  epfo.0 <- fread(grep('full_output', list.files(output_dir_EP, full.names = TRUE), value = TRUE), data.table = F, na.strings = '-9999.0', skip = 3)
  names(epfo.0) <- scan(grep('full_output', list.files(output_dir_EP, full.names = TRUE), value = TRUE), what = 'character', skip = 1, nlines = 1, sep = ',', quiet = T)
  epfo.0$'TIMESTAMP' <- as.character(format(strptime(paste(epfo.0$'date', epfo.0$'time'), '%Y-%m-%d %H:%M', tz = 'GMT'), '%Y%m%d%H%M'))
  
  # align the EP file
  epfo <- data.frame()
  epfo <- base::merge(epfo.0, TIMESTAMP_ref, by = 'TIMESTAMP', all.y = TRUE)
  rm(epfo.0)
  
  # . Import the wavelet combined file ------------------------------------------
  wlfo.0 <- data.frame()
  wlfo.0 <- fread(grep('CDWT_full_cospectra', list.files(output_dir_WL, full.names = TRUE), value = TRUE), data.table = F, na.strings = '-9999')
  wlfo.0$'TIMESTAMP' <- format(strptime(wlfo.0$'TIMESTAMP', '%Y-%m-%d %H:%M', tz = 'GMT'), '%Y%m%d%H%M')
  
  # align the WL file
  wlfo <- data.frame()
  wlfo <- base::merge(wlfo.0, TIMESTAMP_ref, by = 'TIMESTAMP', all.y = TRUE)
  rm(wlfo.0)
  
  
  # . Replace FC --------------------------------------------------------------
  # correct FCwl (== covariance) with air molar volume
  # FC = 1/vd * w'r' , vd= va*Pa/Pd, Pd=Pa-e [https://www.licor.com/env/support/EddyPro/topics/calculate-flux-level-0.html#Frictionvelocityums1]
  rfcs$'CO2flux' <- NA
  rfcs$'CO2flux' <- wlfo$'wco2'/epfo$'air_molar_volume'
  rfcs$'CO2flux' <- rfcs$'CO2flux'*rfcs$'CO2scf'
  # windows()
  # plot(epfo$'co2_flux', type='o')
  # points(rfcs$'CO2flux', type='o', col=4)
  
  # . Replace LE (correct LEwl with lambda) -----------------------------------
  rfcs$'LE' <- NA
  rfcs$'LE' <- wlfo$'wh2o'/epfo$'air_molar_volume'*(10^-3)*0.01802*(10^3*(3147.5-2.37*epfo$'air_temperature'))
  rfcs$'LE' <- rfcs$'LE'*rfcs$'H2Oscf'
  # windows()
  # plot(epfo$LE, type='o')
  # points(rfcs$'LE', type='o', col=4)
  
  # . Replace H (correct Hwl with ?) ------------------------------------------
  rfcs$'H' <- NA
  rfcs$'H' <- wlfo$'wts'*epfo$'air_heat_capacity'*epfo$'air_density'
  rfcs$'H' <- rfcs$'H'*rfcs$'Hscf'
  # windows()
  # plot(epfo$H, type='o')
  # points(rfcs$'H', type='o', col=4)
  
  # . Replace Tau (correct Tauwl with ?) ------------------------------------------
  rfcs$'Tau' <- NA
  rfcs$'Tau' <- epfo$'air_density'*(wlfo$'uw'^2 * wlfo$'vw'^2)^(1/2)
  # windows()
  # plot(epfo$Tau, type='o') # what is Tau in EP output?
  # points(rfcs$'Tau', type='o', col=4)
  
  # . Save WorkSet_WL file (within the WL folder) ----------------------------
  data.table::fwrite(rfcs, paste0(output_dir_WL, site_ID, '_WorkSet_WL.csv'), quote = F, sep=',', na = '-9999')
  
}



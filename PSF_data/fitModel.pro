pro fitModel
   ;This will fit a kinematic model based on OII MUSE observations

    fito=2
   version='1.0'
   model='slp'
   options='_mclean5.0'
   option1='_ssmooth'
   file='input_fit_o2.txt'
   path='o2/'
   line='_OII3729'
   fit_massiv, fito, /ifix, pfix=0, /xfix, /yfix, plot=0, file=file, options=options, option1=option1, path=path, line=line
end

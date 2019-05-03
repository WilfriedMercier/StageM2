pro batch_kin_analysis

    ; correspondance de noms
    ; corrfile='/media/bepinat/WD2To/data/Instruments/MUSE/analyse/Commissioning/tables/names_correspondances_tc_jb.txt'
    
;     On commence par le nettoyage.
;
; xxxxxxxxxxxxxx
;     - cleanning fait avec le programme python 'clean_data.py' + PyQubeVis (clean.fits) + 2eme passe de clean_data.py pour les mclean
;     - puis selection basée sur le SNR a partir du programme 'sample_selection.py'
; xxxxxxxxxxxxxx
    
    ; Et on s'occupe du fit
    
    ; Cas de la raie OII
    
    ; il faut créer le fichier input_fit_o2.txt qui doit contenir
    ; gal        x     y          pa         i             vs         vm    d  sig  psfx psfz  smooth
;     
;     on utilise 'kinemorpho_catalogs.py' (convert_morphocat_kincat) pour déterminer x, y, pa, et i et ainsi faciliter la création de ce fichier
    
    fito=2
    version='1.0'
    model='slp'
    options='_mclean5.0'
    option1='_ssmooth'
    file='input_fit_o2.txt'
    path='o2/'
    line='_OII3729'
    fit_massiv, fito, /ifix, pfix=0, /xfix, /yfix, plot=0, file=file, options=options, option1=option1, path=path, line=line
    
    ; Cas des raies OIII + Hb
    
    ; il faut créer le fichier input_fit_o3hb.txt qui doit contenir
    ; gal        x     y          pa         i             vs         vm    d  sig  psfx psfz  smooth
    
    fito=2
    version='1.0'
    model='slp'
    options='_mclean5.0'
    option1='_ssmooth'
    file='input_fit_o3hb.txt'
    path='o3hb/'
    line='_OIII5007'
    fit_massiv, fito, /ifix, pfix=0, /xfix, /yfix, plot=0, file=file, options=options, option1=option1, path=path, line=line
    
;     pour faire les maps
;     il faut créer le fichier maps_input_o2.txt qui doit contenir (mettre à jour avec les résultats du fit, en particulier pa et vs)
; galaxy_name         name     morpho_image   opt1   line opt2 x    y  pa  i  vs   zpsf   z   r12  psf
    model='slp'
    option1='_ssmooth'
    ;pathgf = '/media/bepinat/WD2To/data/Instruments/MUSE/groups/morpho/acs_mosaic_2.0_gr28/'
    pathgf = '/home/wilfried/hst/x'
    pushd,'o2'
    list='../maps_input_o2.txt'
     create_article_maps_muse,list=list,model=model, option1=option1, pos_name=1
;    create_article_maps_musehdfs_galfit,list=list,model=model, option1=option1, pathgf=pathgf, pos_name=2, field='_CGr30B'
;     pathgp = '/media/bepinat/WD2To/data/Instruments/MUSE/analyse/Commissioning/v031c/galpak/'
;     create_article_maps_musehdfs_galpak_galfit,list=list,model=model, option1=option1, pathgp=pathgp, pathgf=pathgf, pos_name=1
    popd
    pathgf = '/home/vabril/MUSE/Groups_Deep/morphology/HST_CGr30/CGr30_z0.7/'
    pushd,'o3hb'
    list='../maps_input_o3hb.txt'
    ; create_article_maps_muse,list=list,model=model, option1=option1, pos_name=1
    create_article_maps_musehdfs_galfit,list=list,model=model, option1=option1, pathgf=pathgf, pos_name=2, field='_CGr30B'
;     pathgp = '/media/bepinat/WD2To/data/Instruments/MUSE/analyse/Commissioning/v031c/galpak/'
;     create_article_maps_musehdfs_galpak_galfit,list=list,model=model, option1=option1, pathgp=pathgp, pathgf=pathgf, pos_name=1
    popd

end


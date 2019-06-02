pro generateMaps
   ;This creates the final plots comparing the input maps and the kinematical model

   model   ='slp'
   option1 ='_ssmooth'   
   pathgf  = '/home/wilfried/ST2/data/hst/CGr'
   pushd,'o2'
   list    ='../maps_input_o2.txt'
   create_article_maps_muse,list=list,model=model, option1=option1, pos_name=1
   popd
   
end

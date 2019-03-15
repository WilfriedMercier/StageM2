#       #######################################################################
# 
#            			Morphological Catalogue v1.0 
#
#	  Lidia Tasca, Laboratoire d'Astrophysique de Marseille, 23 July 2008
# 
#       ####################################################################### 
#
#
#				~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#			 	  Selection of the catalogue:
#			              F814W(auto)<24.5
#
#			          Total number of objects:
#			              237914
#
#				~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#
#   i)   The catalogue contains morphological information for 237912 objects 
#	 in the ACS catalogue (Leauthaud et al. 2007, apjs, 172, 219)
# 
#   ii)  Morphological parametric quantities have been computed using Morpheus 2005
#        (Abraham , R.~G. et al. 2007, apj, 669, 184)
# 
#   iii) SExtractor (Bertin & Arnouts) is used in a hot/cold configuration to
#        optimimse the deblending
#   
#   iv)  Morphological types are estimated in three different ways.
#        Please do acknowledge the proper paper according to the classification you use:
#        - dist_int    --> Tasca et al. in preparation ; Cassata et al. in preparation
#        - class_linee --> Abraham, R.~G et al. (1996), mnras, 279, L47 ; Cassata et al. in preparation
#        - class_SVM   --> Huertas-Company, M., Rouan, D., Tasca, L., et al., 2008, A&A, 478, 971
# 
#   v)   The catalogue contains  237912 lines of data in 12 columns
# 
# 
#	 MORE IN DETAILS:
#	 ----------------
#	 The morphological study is performed on HST-ACS images 0.03"/pix.
#	 Particular attention has been put in the tuning of the SExtractor
#	 parameters used to produce segmentation images. 
#	 The reason is that Morpheus (Bob's code) computes morphological 
#	 measurements starting from SExtractor segmentation images. 
#	 It is therefore fundamental to reduce cases where galaxies are not deblended.  
#	 Since this morphological catalogue is HST-ACS selected we used the
#	 same hot/cold method as Leauthaud et al. (2007) in order to provide
#	 a one-one match with the official HST/ACS catalogue.
#	 
#	 To transform morphological parameters into morphological classes is
#	 not a trivial job.
#	 To help the user we provide three different morphological classification
#	 all obtained starting from structural parameters measured with Morpheus.
#	 Our preference is for the class named "dist_int". This method has already been
#	 used in Cassata et al. (2007) and optimised for this release. 
#	 The method is explained in Tasca et al. (in prep) and Cassata et al. (in prep).
#	 Another classification based on a support vector machine is the
#	 parameter named "class_SVM". We refer to Huertas-Company et al.(2008)
#	 for further details.
#	 The parameter "class_linee" is the morphological class obtained with an
#	 optimisation of the standard and widely used Abraham et al. (1996) technique 
#	 of subdividing the galaxy population into classes on the basis of 
#	 position on an asymmetry versus concentration diagram. 
#
# 
# Please, contact 
#   lidia.tasca@oamp.fr
#   for any questions related to this catalogue.
# 
# 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#
#
# 
# 1 ID 				my ID, useless for others (sorry!)
# 2 ACS_ALPHA_J2000
# 3 ACS_DELTA_J2000
# 4 ACS_MAG_AUTO (F814W)
# 5 ba				axis ratio
# 6 hlr				half light radius
# 7 g				gini coefficient
# 8 c				concentration
# 9 a				asymmetry
# 10 class_int  	        type= 1         Early type
#			        type= 2         Spirals
#			        type= 3         Irregulars
# 11 class_linee		type= 1         Early type
#                               type= 2         Spirals
#                               type= 3         Irregulars
# 12 class_SVM			type= 1         Early type
#                               type= 2         Spirals
#                               type= 3         Irregulars
#
#
#
#~~~~~~~~~~~~~~~~~~
# 
#   Have fun!
#
##################################################################################

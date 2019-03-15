#
# Cassata's morphological catalog (paolo@astro.umass.edu for any question!)
#
#	232022 galaxies up to mag_i(ACS)=25
#
# Concentration, Asymmetry, Gini and M20 are measured within petrosian
# apertures with a self-written code. The petrosian radius is the
# radius rp at which the ratio of the surface brightness at rp to the mean
# surface brightness within rp is equal to 0.2, and it is measured in 
# elliptical apertures.
# The concentration, asymmetry, gini and m20 are defined as in Cassata 
# et al. 2007 and Cassata et al. 2008 (in prep.).
# B/A is the axial ratio measured by SExtractor.
# 
# The morphological parameters are combined to classify galaxies in
# early-types, disks and irregulars using the following technique (see
# Cassata et al. 2008 in prep. for details and Tasca et al. 2008 in
# prep., using the same technique):
# 
# 1. A reference sample of 500 galaxies with known morphological
# parameters is visually classified, and galaxies are assigned to 
# the three classes. This reference sample allows to characterize the 
# parameters space.
# 
# 2. For each new galaxy that must be classified, the distances in the
# parameters space to the 500 reference galaxies are measured. Then, the 
# 11 closest neighbors are picked out.
#
# 3. Each galaxy is assigned to the most frequent class among
# these 11 neighbors, and the frequency of the chosen class is used to
# assign a reliability to this morphological classification (weight=1: all 
# of the 11 closest neighbors are visually classified in the same class; 
# weight~0.3: the three classes are equally represented among the 11 
# closest neighbors).
# 
# columns:
#  1 ID		 
#  2 X_WORLD (right ascension J2000)
#  3 Y_WORLD (declination J2000)
#  4 MAG_AUTO_ACS (magnitude in I_ACS band)
#  5 PETROSIAN RADIUS
#  6 HALF_LIGHT_RADIUS
#  7 CONCENTRATION
#  8 ASYMMETRY
#  9 GINI COEFFICIENT
#  10 M20
#  11 B/A
#  12 AUTO_CLASS (1=early-type; 2=spiral; 3=irregular)
#  13 AUTO_CLASS_WEIGHT (weight for AUTO CLASS)

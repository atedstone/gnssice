
 obs_file
  <S07> rinex/<S07>/<S07><day>0.<S09>o F L1C L2W C1C C2W L1C L2L C1C C2L
  <S08> rinex/<S08>/<S08><day>0.<S09>o K L1C L2L C1C C2L
 
 nav_file sp3/com<S09><day>.sp3 sp3

 tr_gnss GRE

 interval 10

 mode long

 ion_stats <S06>

 out_type NEU+GEOD
 
 back_type smooth

 SEARCH_TYPE None
 float_type 1 1 LC 0.25 0.5 <S04> <S05> 25

 site_pos
  <S07> 1586032.319 -1932258.392 5848546.971
  <S08> <S01> -<S02> <S03>  
   
 site_stats
 <S08> 10 10 10 1 1 1
# <S08> 100 100 100 10 10 10 
 
 bf_set 2 40

 IONEX_FILE ionex/igsg<day>0.<S09>i

 dcb_file config/dcb.dat

 ante_off 
  <S07> .000 .000 .000 LEIAR25.R4 C
  <S08> .000 .000 .000 SFETOP106 C

 antmod_file config/igs14_2215_plus.atx

 atm_modelc VMF3

 vmf_dir vmf/


 obs_file
  <S07> rinex/<S07>/<S07><day>0.<S09>o F
  <S08> rinex/<S08>/<S08><day>0.<S09>o K
 
 nav_file sp3/igs<S09><day>0.sp3 sp3

 interval 10

 mode long

 ion_stats <S06>

 out_type NEU+GEOD
 
 back_type smooth

 SEARCH_TYPE None
 float_type 1 1 LC 0.25 0.5 <S04> <S05> 25

 site_pos
  <S07> 1594640.5090 -1905874.6780 5855003.8940
  <S08> <S01> -<S02> <S03>  
   
 site_stats
 <S08> 10 10 10 1 1 1
# <S08> 100 100 100 10 10 10 
 
 bf_set 2 40

 IONEX_FILE ionex/igsg<day>0.<S09>i

 dcb_file /work/atedstone/gnss_config/dcb.dat

 ante_off 
  <S07> .000 .000 .000 TWIVC6150 C
  <S08> .000 .000 .000 LEIAX1202GG C

 antmod_file /work/atedstone/gnss_config/igs14_2215_plus.atx

 atm_modelc VMF3

 vmf_dir vmf/

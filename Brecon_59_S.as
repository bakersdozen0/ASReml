!ARGS 3 Dre_28 1.2 !XML #Arguments
Brecon_59 SpatialR
 Genotype_id  Family_id !I Family_name !A Sca_id * Mum_type !A Dad_type !A Crosstype * Trial_id !A Stem_id !I 1680
 Ploc !A Prow * Ppos *
 Block * Av_20  Br_11  Dm_11  Dm_20  Ht_06  Pil_11  St_11  Sur_03  Sur_06  Sur_11  Sur_20  Sur_28  Dbhob_28  Dbhub_28  Drs_28  Dre_28  Bte_28  Cull_28  
 Origin !G 8
 #Ro_north_hg  Ro_north_unk  Ro_north_wc  Ro_south_unk  Ro_hg  Ro_ledmore  Ro_wc  Ro_filler
Brecon_59_S.csv !SKIP 1 !CSV !MAXIT 40 !DOPART $A !EXTRA 5 !MVINCLUDE !NOREORDER !DDF 1

!PART 1
# Design model
$B ~  Origin mv !r Block Family_id uni(Crosstype,2) !GU
1 2 # There are 1 Dummys, each of 2 dimensions
Prow Prow IDEN # Rows
Ppos Ppos IDEN # Columns

!PART 2
# Design+ model
$B ~  Origin mv !r Block Block.Prow Block.Ppos Family_id uni(Crosstype,2) !GU
1 2 # There are 1 Dummys, each of 2 dimensions
Prow Prow IDEN # Rows
Ppos Ppos IDEN # Columns

!PART 3
# Design+ + Spatial AR1 model
$B ~  Origin mv !r Block Block.Prow Block.Ppos Family_id uni(Crosstype,2) !GU units 10
1 2 # There are 1 Dummys, each of 2 dimensions
Prow Prow AR 0.8 !S2=$C # Rows. $C=Ve/10
Ppos Ppos AR 0.8 # Columns

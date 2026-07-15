# datasets.R
#
# Central registry for all circadian datasets.
# To add a new dataset, add ONE entry to CIRC_DATASETS below.
# All other files (server.R, ui.R, comparison_functions.R, explorer_module.R,
# precompute_positions.R) derive their behaviour from this registry automatically.
#
# Field reference:
#   tissue          — tissue id used as grouping key in dropdowns
#   tissue_label    — human-readable tissue name shown in dropdown
#   dataset_label   — human-readable dataset name shown in dropdown
#   folder          — path to the data folder (relative to app root)
#   prefix          — file prefix for network files ({prefix}_nodes_with_tooltips.rds, etc.)
#   period          — oscillation period in hours
#   t               — time-point indices passed to fit_and_generate_curve(); NULL for atlas
#   columnsToFormat — column indices to format in the sidebar DT
#   col_type        — "single"     → BH.Q / phase columns in display_dt
#                     "wt_primary" → BH.Q_wt / phase_wt columns
#                     "ko_primary" → BH.Q_ko / phase_ko columns
#   plot_fn         — "osmo"    → fitted_plot_osmo()
#                     "atlas"   → fitted_plot_atlas()   (also triggers atlas legend)
#                     "cart_wt" → fitted_plot_cart(), is_primary_wt = TRUE
#                     "cart_ko" → fitted_plot_cart(), is_primary_wt = FALSE
#                     "gobs"    → fitted_plot_gobs()
#                     "zt24"    → fitted_plot_zt24()  (6 tp × 4h, triplicates, 24h period)
#   wt_data_file,
#   ko_data_file,
#   display_file    — explicit paths for paired WT/KO datasets (col_type != "single")
#   species         — key into CIRC_SPECIES (e.g. "mouse", "baboon")

# ============================================================
# Species registry
# Add a new species here; set has_data = FALSE for a placeholder.
# ============================================================
CIRC_SPECIES <- list(
  mouse  = list(label = "Mouse (Mus musculus)",  has_data = TRUE),
  baboon = list(label = "Baboon (Papio anubis)", has_data = TRUE)
)

CIRC_DATASETS <- list(

  # ---- Cartilage -----------------------------------------------------------
  cartilage = list(
    tissue          = "cartilage",
    tissue_label    = "Cartilage",
    dataset_label   = "Hip articular WT",
    folder          = "data/mouse/cartilage_circ",
    prefix          = "cartilage",
    period          = 6,
    t               = 1:12,
    columnsToFormat = c(2, 3, 4, 5, 6, 7),
    col_type        = "wt_primary",
    plot_fn         = "cart_wt",
    wt_data_file    = "data/mouse/cartilage_circ/cartilage_wt_data.rds",
    ko_data_file    = "data/mouse/cartilage_circ/cartilage_ko_data.rds",
    display_file    = "data/mouse/cartilage_circ/cartilage_display_dt.rds",
    species         = "mouse"
  ),

  ko_cartilage = list(
    tissue          = "cartilage",
    tissue_label    = "Cartilage",
    dataset_label   = "Hip articular Bmal1 KO",
    folder          = "data/mouse/cartilage_circ",
    prefix          = "ko_cart",
    period          = 6,
    t               = 1:12,
    columnsToFormat = c(2, 3, 4, 5, 6, 7),
    col_type        = "ko_primary",
    plot_fn         = "cart_ko",
    wt_data_file    = "data/mouse/cartilage_circ/cartilage_wt_data.rds",
    ko_data_file    = "data/mouse/cartilage_circ/cartilage_ko_data.rds",
    display_file    = "data/mouse/cartilage_circ/cartilage_display_dt.rds",
    species         = "mouse"
  ),
  # ---- Xiphoid cartilage -------------------------------------------------------
  xiphoid = list(
    tissue          = "cartilage",
    tissue_label    = "cartilage",
    dataset_label   = "Xiphoid cartilage",
    folder          = "data/mouse/xiphoid_cartilage",

    prefix          = "xiphoid",
    period          = 6,
    t               = 1:11,
    columnsToFormat = c(2, 3, 4),
    col_type        = "single",
    plot_fn         = "cart_wt",
    species         = "mouse"
  ),

  # ---- Chondrocytes (Added under cartilage tissue)-------------------------------
  osmo = list(
    tissue          = "cartilage",
    tissue_label    = "cartilage",
    dataset_label   = "Primary Chondrocytes - Osmotic stress",
    folder          = "data/mouse/osmo",
    prefix          = "osmo",
    period          = 6,
    t               = 1:12,
    columnsToFormat = c(2, 3, 4),
    col_type        = "single",
    plot_fn         = "osmo",
    species         = "mouse"
  ),

  dex = list(
    tissue          = "cartilage",
    tissue_label    = "cartilage",
    dataset_label   = "Primary Chondrocytes - Dexamethasone",
    folder          = "data/mouse/dex",
    prefix          = "dex",
    period          = 6,
    t               = 1:12,
    columnsToFormat = c(2, 3, 4),
    col_type        = "single",
    plot_fn         = "osmo",
    species         = "mouse"
  ),

  hs = list(
    tissue          = "cartilage",
    tissue_label    = "cartilage",
    dataset_label   = "Primary Chondrocytes - Heat shock",
    folder          = "data/mouse/hs",
    prefix          = "hs",
    period          = 6,
    t               = 1:12,
    columnsToFormat = c(2, 3, 4),
    col_type        = "single",
    plot_fn         = "osmo",
    species         = "mouse"
  ),

  # ---- Intervertebral Disc -------------------------------------------------
  # ivd = list(
  #   tissue          = "ivd",
  #   tissue_label    = "Intervertebral Disc",
  #   dataset_label   = "WT",
  #   folder          = "data/mouse/ivd_circ",
  #   prefix          = "ivd",
  #   period          = 6,
  #   t               = 1:12,
  #   columnsToFormat = c(2, 3, 4, 5, 6, 7),
  #   col_type        = "wt_primary",
  #   plot_fn         = "cart_wt",
  #   wt_data_file    = "data/mouse/ivd_circ/ivd_wt_data.rds",
  #   ko_data_file    = "data/mouse/ivd_circ/ivd_ko_data.rds",
  #   display_file    = "data/mouse/ivd_circ/ivd_display_dt.rds",
  #   species         = "mouse"
  # ),

  # ko_ivd = list(
  #   tissue          = "ivd",
  #   tissue_label    = "Intervertebral Disc",
  #   dataset_label   = "Bmal1 KO",
  #   folder          = "data/mouse/ivd_circ",
  #   prefix          = "ko_ivd",
  #   period          = 6,
  #   t               = 1:12,
  #   columnsToFormat = c(2, 3, 4, 5, 6, 7),
  #   col_type        = "ko_primary",
  #   plot_fn         = "cart_ko",
  #   wt_data_file    = "data/mouse/ivd_circ/ivd_wt_data.rds",
  #   ko_data_file    = "data/mouse/ivd_circ/ivd_ko_data.rds",
  #   display_file    = "data/mouse/ivd_circ/ivd_display_dt.rds",
  #   species         = "mouse"
  # ),

  # "12h_ivd" = list(
  #   tissue          = "ivd",
  #   tissue_label    = "Intervertebral Disc",
  #   dataset_label   = "Bmal1 KO 12h rhythm",
  #   folder          = "data/mouse/12h_ivd",
  #   prefix          = "12h_ko_ivd",
  #   period          = 3,
  #   t               = 1:12,
  #   columnsToFormat = c(2, 3, 4, 5, 6, 7),
  #   col_type        = "ko_primary",
  #   plot_fn         = "cart_ko",
  #   wt_data_file    = "data/mouse/12h_ivd/12h_ivd_wt_data.rds",
  #   ko_data_file    = "data/mouse/12h_ivd/12h_ivd_ko_data.rds",
  #   display_file    = "data/mouse/12h_ivd/12h_ivd_display_dt.rds",
  #   species         = "mouse"
  # ),
  # ---- NIH3T3 Forskolin -------------------------------------------------------
  nih3t3 = list(
    tissue          = "nih3t3",
    tissue_label    = "NIH3T3 Forskolin",
    dataset_label   = "NIH3T3 Forskolin",
    folder          = "data/mouse/nih3t3",
    prefix          = "nih3t3",
    period          = 24,
    t               = NULL,
    columnsToFormat = c(2, 3, 4),
    col_type        = "single",
    plot_fn         = "atlas",
    species         = "mouse"
  ),
  # ---- Mammary gland -------------------------------------------------------
  mammary = list(
    tissue          = "mammary_gland",
    tissue_label    = "Mammary gland",
    dataset_label   = "Mammary gland",
    folder          = "data/mouse/mammary_gland",
    prefix          = "mammary_gland",
    period          = 6,
    t               = 1:12,
    columnsToFormat = c(2, 3, 4),
    col_type        = "single",
    plot_fn         = "cart_wt",
    species         = "mouse"
  ),
  # ---- Tendon --------------------------------------------------------------
  tendon = list(
    tissue          = "tendon",
    tissue_label    = "Tendon",
    dataset_label   = "Tail tendon",
    folder          = "data/mouse/tendon",
    prefix          = "tendon",
    period          = 6,
    t               = 1:6,
    columnsToFormat = c(2, 3, 4),
    col_type        = "single",
    plot_fn         = "osmo",
    species         = "mouse"
  ),

  # ---- Goblet cells --------------------------------------------------------
  # gobs = list(
  #   tissue          = "gobs",
  #   tissue_label    = "Goblet cells",
  #   dataset_label   = "Goblet cells",
  #   folder          = "data/mouse/gobs",
  #   prefix          = "gobs",
  #   period          = 6,
  #   t               = 1:12,
  #   columnsToFormat = c(2, 3, 4),
  #   col_type        = "single",
  #   plot_fn         = "gobs",
  #   species         = "mouse"
  # ),

  # ---- Kidney (grouped: podocytes  + glomeruli + whole kidney) -------------
  pods = list(
    tissue          = "kidney",
    tissue_label    = "Kidney",
    dataset_label   = "Podocytes dexamethasone",
    folder          = "data/mouse/podocytes_dexamethasone",
    prefix          = "pods",
    period          = 6,
    t               = 1:12,
    columnsToFormat = c(2, 3, 4),
    col_type        = "single",
    plot_fn         = "osmo",
    species         = "mouse"
  ),

  # ---- Kidney (grouped: podocytes + glomeruli + whole kidney) --------------
  gloms = list(
    tissue          = "kidney",
    tissue_label    = "Kidney",
    dataset_label   = "Glomeruli",
    folder          = "data/mouse/gloms",
    prefix          = "gloms",
    period          = 6,
    t               = 1:12,
    columnsToFormat = c(2, 3, 4),
    col_type        = "single",
    plot_fn         = "cart_wt",
    species         = "mouse"
  ),
  # ---- Liver RNAseq Ad Libidum and Night Restricted Feeding GSE158600-----------------
  liver_AL_WT = list(
    tissue          = "liver",
    tissue_label    = "Liver",
    dataset_label   = "Liver WT Ad Libidum",
    folder          = "data/mouse/liver_RNAseq/AL_WT",
    prefix          = "AL_WT",
    period          = 24,
    t               = 1:6,
    columnsToFormat = c(2, 3, 4),
    col_type        = "single",
    plot_fn         = "zt24",
    species         = "mouse"
  ),

  liver_AL_KO = list(
    tissue          = "liver",
    tissue_label    = "Liver",
    dataset_label   = "Liver KO Ad Libidum",
    folder          = "data/mouse/liver_RNAseq/AL_KO",
    prefix          = "AL_KO",
    period          = 24,
    t               = 1:6,
    columnsToFormat = c(2, 3, 4),
    col_type        = "single",
    plot_fn         = "zt24",
    species         = "mouse"
  ),

  liver_AL_RE = list(
    tissue          = "liver",
    tissue_label    = "Liver",
    dataset_label   = "Liver RE Ad Libidum",
    folder          = "data/mouse/liver_RNAseq/AL_RE",
    prefix          = "AL_RE",
    period          = 24,
    t               = 1:6,
    columnsToFormat = c(2, 3, 4),
    col_type        = "single",
    plot_fn         = "zt24",
    species         = "mouse"
  ),

  liver_TRF_WT = list(
    tissue          = "liver",
    tissue_label    = "Liver",
    dataset_label   = "Liver TRF WT",
    folder          = "data/mouse/liver_RNAseq/TRF_WT",
    prefix          = "TRF_WT",
    period          = 24,
    t               = 1:6,
    columnsToFormat = c(2, 3, 4),
    col_type        = "single",
    plot_fn         = "zt24",
    species         = "mouse"
  ),

  liver_TRF_KO = list(
    tissue          = "liver",
    tissue_label    = "Liver",
    dataset_label   = "Liver TRF KO",
    folder          = "data/mouse/liver_RNAseq/TRF_KO",
    prefix          = "TRF_KO",
    period          = 24,
    t               = 1:6,
    columnsToFormat = c(2, 3, 4),
    col_type        = "single",
    plot_fn         = "zt24",
    species         = "mouse"
  ),

  liver_TRF_RE = list(
    tissue          = "liver",
    tissue_label    = "Liver",
    dataset_label   = "Liver TRF RE",
    folder          = "data/mouse/liver_RNAseq/TRF_RE",
    prefix          = "TRF_RE",
    period          = 24,
    t               = 1:6,
    columnsToFormat = c(2, 3, 4),
    col_type        = "single",
    plot_fn         = "zt24",
    species         = "mouse"
  ),

  # ---- Atlas tissues (period = 24, t = NULL) --------------------------------
  adrenal = list(
    tissue          = "adrenal_gland",
    tissue_label    = "Adrenal gland",
    dataset_label   = "Adrenal gland",
    folder          = "data/mouse/adrenal_gland",
    prefix          = "adrenal_gland",
    period          = 24,
    t               = NULL,
    columnsToFormat = c(2, 3, 4),
    col_type        = "single",
    plot_fn         = "atlas",
    species         = "mouse"
  ),

  aorta = list(
    tissue          = "aorta",
    tissue_label    = "Aorta",
    dataset_label   = "Aorta",
    folder          = "data/mouse/aorta",
    prefix          = "aorta",
    period          = 24,
    t               = NULL,
    columnsToFormat = c(2, 3, 4),
    col_type        = "single",
    plot_fn         = "atlas",
    species         = "mouse"
  ),

  brown_adipose = list(
    tissue          = "brown_adipose",
    tissue_label    = "Brown adipose",
    dataset_label   = "Brown adipose",
    folder          = "data/mouse/brown_adipose",
    prefix          = "brown_adipose",
    period          = 24,
    t               = NULL,
    columnsToFormat = c(2, 3, 4),
    col_type        = "single",
    plot_fn         = "atlas",
    species         = "mouse"
  ),

  brain_stem = list(
    tissue          = "brain_stem",
    tissue_label    = "Brain stem",
    dataset_label   = "Brain stem",
    folder          = "data/mouse/brain_stem",
    prefix          = "brain_stem",
    period          = 24,
    t               = NULL,
    columnsToFormat = c(2, 3, 4),
    col_type        = "single",
    plot_fn         = "atlas",
    species         = "mouse"
  ),

  cerebellum = list(
    tissue          = "cerebellum",
    tissue_label    = "Cerebellum",
    dataset_label   = "Cerebellum",
    folder          = "data/mouse/cerebellum",
    prefix          = "cerebellum",
    period          = 24,
    t               = NULL,
    columnsToFormat = c(2, 3, 4),
    col_type        = "single",
    plot_fn         = "atlas",
    species         = "mouse"
  ),

  heart = list(
    tissue          = "heart",
    tissue_label    = "Heart",
    dataset_label   = "Heart",
    folder          = "data/mouse/heart",
    prefix          = "heart",
    period          = 24,
    t               = NULL,
    columnsToFormat = c(2, 3, 4),
    col_type        = "single",
    plot_fn         = "atlas",
    species         = "mouse"
  ),

  hypothalamus = list(
    tissue          = "hypothalamus",
    tissue_label    = "Hypothalamus",
    dataset_label   = "Hypothalamus",
    folder          = "data/mouse/hypothalamus",
    prefix          = "hypothalamus",
    period          = 24,
    t               = NULL,
    columnsToFormat = c(2, 3, 4),
    col_type        = "single",
    plot_fn         = "atlas",
    species         = "mouse"
  ),

  kidney = list(
    tissue          = "kidney",
    tissue_label    = "Kidney",
    dataset_label   = "Kidney",
    folder          = "data/mouse/kidney",
    prefix          = "kidney",
    period          = 24,
    t               = NULL,
    columnsToFormat = c(2, 3, 4),
    col_type        = "single",
    plot_fn         = "atlas",
    species         = "mouse"
  ),

  liver = list(
    tissue          = "liver",
    tissue_label    = "Liver",
    dataset_label   = "Liver",
    folder          = "data/mouse/liver",
    prefix          = "liver",
    period          = 24,
    t               = NULL,
    columnsToFormat = c(2, 3, 4),
    col_type        = "single",
    plot_fn         = "atlas",
    species         = "mouse"
  ),

  lung = list(
    tissue          = "lung",
    tissue_label    = "Lung",
    dataset_label   = "Lung",
    folder          = "data/mouse/lung",
    prefix          = "lung",
    period          = 24,
    t               = NULL,
    columnsToFormat = c(2, 3, 4),
    col_type        = "single",
    plot_fn         = "atlas",
    species         = "mouse"
  ),

  muscle = list(
    tissue          = "muscle",
    tissue_label    = "Muscle",
    dataset_label   = "Muscle",
    folder          = "data/mouse/muscle",
    prefix          = "muscle",
    period          = 24,
    t               = NULL,
    columnsToFormat = c(2, 3, 4),
    col_type        = "single",
    plot_fn         = "atlas",
    species         = "mouse"
  ),

  white_adipose = list(
    tissue          = "white_adipose",
    tissue_label    = "White adipose",
    dataset_label   = "White adipose",
    folder          = "data/mouse/white_adipose",
    prefix          = "white_adipose",
    period          = 24,
    t               = NULL,
    columnsToFormat = c(2, 3, 4),
    col_type        = "single",
    plot_fn         = "atlas",
    species         = "mouse"
  ),

  # ---- Baboon --------------------------------------------------------------
  baboon_adrenal_cortex = list(
    tissue = "adrenal_cortex", tissue_label = "Adrenal cortex", dataset_label = "Adrenal cortex",
    folder = "data/baboon/adrenal_cortex", prefix = "adrenal_cortex",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_adrenal_medulla = list(
    tissue = "adrenal_medulla", tissue_label = "Adrenal medulla", dataset_label = "Adrenal medulla",
    folder = "data/baboon/adrenal_medulla", prefix = "adrenal_medulla",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_amygdala = list(
    tissue = "amygdala", tissue_label = "Amygdala", dataset_label = "Amygdala",
    folder = "data/baboon/amygdala", prefix = "amygdala",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_antrum = list(
    tissue = "antrum", tissue_label = "Antrum", dataset_label = "Antrum",
    folder = "data/baboon/antrum", prefix = "antrum",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_aorta_endothelium = list(
    tissue = "aorta_endothelium", tissue_label = "Aorta endothelium", dataset_label = "Aorta endothelium",
    folder = "data/baboon/aorta_endothelium", prefix = "aorta_endothelium",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_arcuate_nucleus = list(
    tissue = "arcuate_nucleus", tissue_label = "Arcuate nucleus", dataset_label = "Arcuate nucleus",
    folder = "data/baboon/arcuate_nucleus", prefix = "arcuate_nucleus",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_ascending_colon = list(
    tissue = "ascending_colon", tissue_label = "Ascending colon", dataset_label = "Ascending colon",
    folder = "data/baboon/ascending_colon", prefix = "ascending_colon",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_axillary_lymphonodes = list(
    tissue = "axillary_lymphonodes", tissue_label = "Axillary lymph nodes", dataset_label = "Axillary lymph nodes",
    folder = "data/baboon/axillary_lymphonodes", prefix = "axillary_lymphonodes",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_bladder = list(
    tissue = "bladder", tissue_label = "Bladder", dataset_label = "Bladder",
    folder = "data/baboon/bladder", prefix = "bladder",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_bone_marrow = list(
    tissue = "bone_marrow", tissue_label = "Bone marrow", dataset_label = "Bone marrow",
    folder = "data/baboon/bone_marrow", prefix = "bone_marrow",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_cecum = list(
    tissue = "cecum", tissue_label = "Cecum", dataset_label = "Cecum",
    folder = "data/baboon/cecum", prefix = "cecum",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_cerebellum = list(
    tissue = "cerebellum", tissue_label = "Cerebellum", dataset_label = "Cerebellum",
    folder = "data/baboon/cerebellum", prefix = "cerebellum",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_cornea = list(
    tissue = "cornea", tissue_label = "Cornea", dataset_label = "Cornea",
    folder = "data/baboon/cornea", prefix = "cornea",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_descending_colon = list(
    tissue = "descending_colon", tissue_label = "Descending colon", dataset_label = "Descending colon",
    folder = "data/baboon/descending_colon", prefix = "descending_colon",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_dorsomedial_hypothalamus = list(
    tissue = "dorsomedial_hypothalamus", tissue_label = "Dorsomedial hypothalamus", dataset_label = "Dorsomedial hypothalamus",
    folder = "data/baboon/dorsomedial_hypothalamus", prefix = "dorsomedial_hypothalamus",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_duodenum = list(
    tissue = "duodenum_opposite_of_ampula", tissue_label = "Duodenum (opp. ampulla)", dataset_label = "Duodenum (opp. ampulla)",
    folder = "data/baboon/duodenum_opposite_of_ampula", prefix = "duodenum_opposite_of_ampula",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_habenula = list(
    tissue = "habenula", tissue_label = "Habenula", dataset_label = "Habenula",
    folder = "data/baboon/habenula", prefix = "habenula",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_heart = list(
    tissue = "heart", tissue_label = "Heart", dataset_label = "Heart",
    folder = "data/baboon/heart", prefix = "heart",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_hippocampus = list(
    tissue = "hippocampus", tissue_label = "Hippocampus", dataset_label = "Hippocampus",
    folder = "data/baboon/hippocampus", prefix = "hippocampus",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_ileum = list(
    tissue = "ileum", tissue_label = "Ileum", dataset_label = "Ileum",
    folder = "data/baboon/ileum", prefix = "ileum",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_iris = list(
    tissue = "iris", tissue_label = "Iris", dataset_label = "Iris",
    folder = "data/baboon/iris", prefix = "iris",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_kidney_cortex = list(
    tissue = "kidney_cortex", tissue_label = "Kidney cortex", dataset_label = "Kidney cortex",
    folder = "data/baboon/kidney_cortex", prefix = "kidney_cortex",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_kidney_medulla = list(
    tissue = "kidney_medulla", tissue_label = "Kidney medulla", dataset_label = "Kidney medulla",
    folder = "data/baboon/kidney_medulla", prefix = "kidney_medulla",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_lateral_globus_pallidus = list(
    tissue = "lateral_globus_pallidus", tissue_label = "Lateral globus pallidus", dataset_label = "Lateral globus pallidus",
    folder = "data/baboon/lateral_globus_pallidus", prefix = "lateral_globus_pallidus",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_lateral_hypothalamus = list(
    tissue = "lateral_hypothalamus", tissue_label = "Lateral hypothalamus", dataset_label = "Lateral hypothalamus",
    folder = "data/baboon/lateral_hypothalamus", prefix = "lateral_hypothalamus",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_liver = list(
    tissue = "liver", tissue_label = "Liver", dataset_label = "Liver",
    folder = "data/baboon/liver", prefix = "liver",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_lung = list(
    tissue = "lung", tissue_label = "Lung", dataset_label = "Lung",
    folder = "data/baboon/lung", prefix = "lung",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_mammilary_bodies = list(
    tissue = "mammilary_bobies", tissue_label = "Mammillary bodies", dataset_label = "Mammillary bodies",
    folder = "data/baboon/mammilary_bobies", prefix = "mammilary_bobies",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_medial_globus_pallidus = list(
    tissue = "medial_globus_pallidus", tissue_label = "Medial globus pallidus", dataset_label = "Medial globus pallidus",
    folder = "data/baboon/medial_globus_pallidus", prefix = "medial_globus_pallidus",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_mesenteric_lymphonodes = list(
    tissue = "mesenteric_lymphonodes", tissue_label = "Mesenteric lymph nodes", dataset_label = "Mesenteric lymph nodes",
    folder = "data/baboon/mesenteric_lymphonodes", prefix = "mesenteric_lymphonodes",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_muscle_abdominal = list(
    tissue = "muscle_abdominal", tissue_label = "Muscle (abdominal)", dataset_label = "Muscle (abdominal)",
    folder = "data/baboon/muscle_abdominal", prefix = "muscle_abdominal",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_muscle_gastrocnemian = list(
    tissue = "muscle_gastrocnemian", tissue_label = "Muscle (gastrocnemius)", dataset_label = "Muscle (gastrocnemius)",
    folder = "data/baboon/muscle_gastrocnemian", prefix = "muscle_gastrocnemian",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_oesophagus = list(
    tissue = "oesophagus", tissue_label = "Oesophagus", dataset_label = "Oesophagus",
    folder = "data/baboon/oesophagus", prefix = "oesophagus",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_olfactory_bulb = list(
    tissue = "olfactory_bulb", tissue_label = "Olfactory bulb", dataset_label = "Olfactory bulb",
    folder = "data/baboon/olfactory_bulb", prefix = "olfactory_bulb",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_omental_fat = list(
    tissue = "omental_fat", tissue_label = "Omental fat", dataset_label = "Omental fat",
    folder = "data/baboon/omental_fat", prefix = "omental_fat",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_optic_nerve_head = list(
    tissue = "optic_nerve_head", tissue_label = "Optic nerve head", dataset_label = "Optic nerve head",
    folder = "data/baboon/optic_nerve_head", prefix = "optic_nerve_head",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_pancreas = list(
    tissue = "pancreas", tissue_label = "Pancreas", dataset_label = "Pancreas",
    folder = "data/baboon/pancreas", prefix = "pancreas",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_paraventricular_nuclei = list(
    tissue = "paraventricular_nuclei", tissue_label = "Paraventricular nuclei", dataset_label = "Paraventricular nuclei",
    folder = "data/baboon/paraventricular_nuclei", prefix = "paraventricular_nuclei",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_pineal = list(
    tissue = "pineal", tissue_label = "Pineal gland", dataset_label = "Pineal gland",
    folder = "data/baboon/pineal", prefix = "pineal",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_pituitary = list(
    tissue = "pituitary", tissue_label = "Pituitary", dataset_label = "Pituitary",
    folder = "data/baboon/pituitary", prefix = "pituitary",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_pons = list(
    tissue = "pons", tissue_label = "Pons", dataset_label = "Pons",
    folder = "data/baboon/pons", prefix = "pons",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_prefrontal_cortex = list(
    tissue = "prefrontal_cortex", tissue_label = "Prefrontal cortex", dataset_label = "Prefrontal cortex",
    folder = "data/baboon/prefrontal_cortex", prefix = "prefrontal_cortex",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_preoptic_area = list(
    tissue = "preoptic_area", tissue_label = "Preoptic area", dataset_label = "Preoptic area",
    folder = "data/baboon/preoptic_area", prefix = "preoptic_area",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_prostate = list(
    tissue = "prostate", tissue_label = "Prostate", dataset_label = "Prostate",
    folder = "data/baboon/prostate", prefix = "prostate",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_putamen = list(
    tissue = "putamen", tissue_label = "Putamen", dataset_label = "Putamen",
    folder = "data/baboon/putamen", prefix = "putamen",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_retina = list(
    tissue = "retina", tissue_label = "Retina", dataset_label = "Retina",
    folder = "data/baboon/retina", prefix = "retina",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_retinal_pigment_epithelium = list(
    tissue = "retinal_pigment_epithelium", tissue_label = "Retinal pigment epithelium", dataset_label = "Retinal pigment epithelium",
    folder = "data/baboon/retinal_pigment_epithelium", prefix = "retinal_pigment_epithelium",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_skin = list(
    tissue = "skin_from_underbelly", tissue_label = "Skin (underbelly)", dataset_label = "Skin (underbelly)",
    folder = "data/baboon/skin_from_underbelly", prefix = "skin_from_underbelly",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_smooth_muscle = list(
    tissue = "smooth_muscle_opposite_to_antrum", tissue_label = "Smooth muscle (opp. antrum)", dataset_label = "Smooth muscle (opp. antrum)",
    folder = "data/baboon/smooth_muscle_opposite_to_antrum", prefix = "smooth_muscle_opposite_to_antrum",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_spleen = list(
    tissue = "spleen", tissue_label = "Spleen", dataset_label = "Spleen",
    folder = "data/baboon/spleen", prefix = "spleen",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_stomach_fundus = list(
    tissue = "stomach_fundus", tissue_label = "Stomach fundus", dataset_label = "Stomach fundus",
    folder = "data/baboon/stomach_fundus", prefix = "stomach_fundus",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_substantia_nigra = list(
    tissue = "substantia_nigra", tissue_label = "Substantia nigra", dataset_label = "Substantia nigra",
    folder = "data/baboon/substantia_nigra", prefix = "substantia_nigra",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_suprachiasmatic_nuclei = list(
    tissue = "suprachiasmatic_nuclei", tissue_label = "Suprachiasmatic nuclei", dataset_label = "Suprachiasmatic nuclei",
    folder = "data/baboon/suprachiasmatic_nuclei", prefix = "suprachiasmatic_nuclei",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_supraoptic_nucleus = list(
    tissue = "supraoptic_nucleus", tissue_label = "Supraoptic nucleus", dataset_label = "Supraoptic nucleus",
    folder = "data/baboon/supraoptic_nucleus", prefix = "supraoptic_nucleus",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_testicles = list(
    tissue = "testicles", tissue_label = "Testicles", dataset_label = "Testicles",
    folder = "data/baboon/testicles", prefix = "testicles",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_thalamus = list(
    tissue = "thalamus", tissue_label = "Thalamus", dataset_label = "Thalamus",
    folder = "data/baboon/thalamus", prefix = "thalamus",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_thyroid = list(
    tissue = "thyroid", tissue_label = "Thyroid", dataset_label = "Thyroid",
    folder = "data/baboon/thyroid", prefix = "thyroid",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_ventromedial_hypothalamus = list(
    tissue = "ventro_medial_hypothalamus", tissue_label = "Ventromedial hypothalamus", dataset_label = "Ventromedial hypothalamus",
    folder = "data/baboon/ventro_medial_hypothalamus", prefix = "ventro_medial_hypothalamus",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_visual_cortex = list(
    tissue = "visual_cortex", tissue_label = "Visual cortex", dataset_label = "Visual cortex",
    folder = "data/baboon/visual_cortex", prefix = "visual_cortex",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_white_adipose_mesenteric = list(
    tissue = "white_adipose_mesenteric", tissue_label = "White adipose (mesenteric)", dataset_label = "White adipose (mesenteric)",
    folder = "data/baboon/white_adipose_mesenteric", prefix = "white_adipose_mesenteric",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_white_adipose_pericardial = list(
    tissue = "white_adipose_pericardial", tissue_label = "White adipose (pericardial)", dataset_label = "White adipose (pericardial)",
    folder = "data/baboon/white_adipose_pericardial", prefix = "white_adipose_pericardial",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_white_adipose_perirenal = list(
    tissue = "white_adipose_perirenal", tissue_label = "White adipose (perirenal)", dataset_label = "White adipose (perirenal)",
    folder = "data/baboon/white_adipose_perirenal", prefix = "white_adipose_perirenal",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_white_adipose_retroperitoneal = list(
    tissue = "white_adipose_retroperitoneal", tissue_label = "White adipose (retroperitoneal)", dataset_label = "White adipose (retroperitoneal)",
    folder = "data/baboon/white_adipose_retroperitoneal", prefix = "white_adipose_retroperitoneal",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  ),
  baboon_white_adipose_subcutaneous = list(
    tissue = "white_adipose_subcutaneous", tissue_label = "White adipose (subcutaneous)", dataset_label = "White adipose (subcutaneous)",
    folder = "data/baboon/white_adipose_subcutaneous", prefix = "white_adipose_subcutaneous",
    period = 24, t = NULL, columnsToFormat = c(2, 3, 4), col_type = "single", plot_fn = "baboon_atlas", species = "baboon"
  )
)

# ============================================================
# Derive dropdown structures from the registry.
# Insertion order in CIRC_DATASETS determines dropdown order.
# ============================================================

build_tissue_choices <- function(registry) {
  seen   <- character(0)
  result <- character(0)
  for (meta in registry) {
    key <- meta$tissue
    if (!key %in% seen) {
      seen   <- c(seen, key)
      result <- c(result, setNames(key, meta$tissue_label))
    }
  }
  result
}

build_dataset_choices_map <- function(registry) {
  result <- list()
  for (ds_id in names(registry)) {
    meta   <- registry[[ds_id]]
    tissue <- meta$tissue
    result[[tissue]] <- c(result[[tissue]],
                          setNames(ds_id, meta$dataset_label))
  }
  result
}

CIRC_TISSUE_CHOICES      <- build_tissue_choices(CIRC_DATASETS)
CIRC_DATASET_CHOICES_MAP <- build_dataset_choices_map(CIRC_DATASETS)

# ---- Species-aware helpers (for the 3-level Species → Tissue → Dataset cascade) ----

build_species_choices <- function(species_registry) {
  setNames(
    names(species_registry),
    vapply(species_registry, `[[`, character(1), "label")
  )
}

build_tissue_choices_by_species <- function(dataset_registry) {
  result <- list()
  for (meta in dataset_registry) {
    sp  <- meta$species
    key <- meta$tissue
    if (is.null(result[[sp]])) result[[sp]] <- character(0)
    if (!key %in% result[[sp]]) {
      result[[sp]] <- c(result[[sp]], setNames(key, meta$tissue_label))
    }
  }
  result
}

build_dataset_choices_by_species <- function(dataset_registry) {
  result <- list()
  for (ds_id in names(dataset_registry)) {
    meta   <- dataset_registry[[ds_id]]
    sp     <- meta$species
    tissue <- meta$tissue
    if (is.null(result[[sp]])) result[[sp]] <- list()
    result[[sp]][[tissue]] <- c(result[[sp]][[tissue]],
                                setNames(ds_id, meta$dataset_label))
  }
  result
}

CIRC_SPECIES_CHOICES           <- build_species_choices(CIRC_SPECIES)
CIRC_TISSUE_CHOICES_BY_SPECIES <- build_tissue_choices_by_species(CIRC_DATASETS)
CIRC_DATASET_CHOICES_BY_SPECIES <- build_dataset_choices_by_species(CIRC_DATASETS)

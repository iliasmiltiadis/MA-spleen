# ==============================================================================
# Meta-analysis of Splenic Dimensions and Vasculature
# ==============================================================================

library(meta)
library(grid)

# 1. Ensure output directory exists
dir.create("figures", showWarnings = FALSE)

# 2. Load Data
dat <- read.csv("data/spleen-database.csv", stringsAsFactors = FALSE)

# 3. Define Custom Color Palette
col_overall <- "#0072b2"  # Blue (Base models)
col_mod     <- "#d55e00"  # Vermillion (Modality subgroups)
col_qa      <- "#009e73"  # Bluish Green (QA subgroups)
col_branch2 <- "#845da7"  # Purple (Specific for 2-branch baseline)

# ==============================================================================
# HELPER FUNCTIONS FOR CLEAN PLOTTING
# ==============================================================================

# Function to draw the 3-panel Forest Plot
export_3panel_forest <- function(m_base, m_mod, m_qa, filename, x_label, sm_label, base_color) {
  pdf(filename, width = 18, height = 17)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(3, 1)))
  
  plot_row <- function(model, row_num, label, color) {
    pushViewport(viewport(layout.pos.row = row_num, layout.pos.col = 1))
    
    if (inherits(model, "metaprop")) {
      forest(model, sortvar = TE, prediction = TRUE, col.square = color, col.square.lines = "black", 
             col.diamond = color, col.diamond.lines = "black", leftlabs = c("Study", "Events", "Total"), 
             xlab = x_label, smlab = sm_label, digits = 2, new = FALSE)
             
    } else if (inherits(model, "metacont")) {
      forest(model, sortvar = TE, prediction = TRUE, col.square = color, col.square.lines = "black", 
             col.diamond = color, col.diamond.lines = "black", label.e = "Male", label.c = "Female", 
             xlab = x_label, smlab = sm_label, digits = 2, new = FALSE)
             
    } else {
      forest(model, sortvar = TE, prediction = TRUE, col.square = color, col.square.lines = "black", 
             col.diamond = color, col.diamond.lines = "black", leftlabs = c("Study", "Mean", "SE"), 
             xlab = x_label, smlab = sm_label, digits = 2, digits.se = 2, new = FALSE)
    }
    
    grid.text(label, x = unit(0.02, "npc"), y = unit(0.92, "npc"), 
              gp = gpar(fontsize = 30, fontfamily = "sans", fontface = "bold"))
    popViewport()
  }
  
  plot_row(m_base, 1, "a", base_color)
  plot_row(m_mod,  2, "b", col_mod)
  plot_row(m_qa,   3, "c", col_qa)
  
  dev.off()
}

# Function to run and export Trim-and-Fill Funnel Plot
export_funnel <- function(model, filename, title) {
  tf_model <- trimfill(model)
  pdf(filename, width = 8, height = 8)
  funnel(tf_model, pch = ifelse(tf_model$trimfill, 1, 16), level = 0.95, 
         contour = c(0.9, 0.95, 0.99), col.contour = c("gray75", "gray85", "gray95"), 
         main = title)
  legend("topright", c("Observed", "Imputed"), pch = c(16, 1), bg = "white")
  dev.off()
}


# ==============================================================================
# SECTION A: SINGLE-ARM CONTINUOUS VARIABLES
# ==============================================================================

# 1. Splenic Length
df_len <- subset(dat, !is.na(Len_M) & !is.na(Len_SE))
m_len <- metagen(TE = Len_M, seTE = Len_SE, studlab = Study, data = df_len, sm = "MD", random = TRUE, common = FALSE, method.tau = "REML")
m_len_mod <- update(m_len, subgroup = Modality, subgroup.name = "Modality")
m_len_qa  <- update(m_len, subgroup = QA, subgroup.name = "Quality")

export_3panel_forest(m_len, m_len_mod, m_len_qa, "figures/Length_Panel.pdf", "Pooled mean", "Pooled mean", col_overall)
export_funnel(m_len, "figures/Length_Funnel.pdf", "Funnel Plot: Splenic Length")


# 2. Splenic Width
df_wid <- subset(dat, !is.na(Wid_M) & !is.na(Wid_SE))
m_wid <- metagen(TE = Wid_M, seTE = Wid_SE, studlab = Study, data = df_wid, sm = "MD", random = TRUE, common = FALSE, method.tau = "REML")
m_wid_mod <- update(m_wid, subgroup = Modality, subgroup.name = "Modality")
m_wid_qa  <- update(m_wid, subgroup = QA, subgroup.name = "Quality")

export_3panel_forest(m_wid, m_wid_mod, m_wid_qa, "figures/Width_Panel.pdf", "Pooled mean", "Pooled mean", col_overall)
export_funnel(m_wid, "figures/Width_Funnel.pdf", "Funnel Plot: Splenic Width")


# 3. Splenic Thickness
df_thick <- subset(dat, !is.na(Thick_M) & !is.na(Thick_SE))
m_thick <- metagen(TE = Thick_M, seTE = Thick_SE, studlab = Study, data = df_thick, sm = "MD", random = TRUE, common = FALSE, method.tau = "REML")
m_thick_mod <- update(m_thick, subgroup = Modality, subgroup.name = "Modality")
m_thick_qa  <- update(m_thick, subgroup = QA, subgroup.name = "Quality")

export_3panel_forest(m_thick, m_thick_mod, m_thick_qa, "figures/Thickness_Panel.pdf", "Pooled mean", "Pooled mean", col_overall)
export_funnel(m_thick, "figures/Thickness_Funnel.pdf", "Funnel Plot: Splenic Thickness")


# 4. Splenic Volume
df_vol <- subset(dat, !is.na(Vol_M) & !is.na(Vol_SE))
m_vol <- metagen(TE = Vol_M, seTE = Vol_SE, studlab = Study, data = df_vol, sm = "MD", random = TRUE, common = FALSE, method.tau = "REML")
m_vol_mod <- update(m_vol, subgroup = Modality, subgroup.name = "Modality")
m_vol_qa  <- update(m_vol, subgroup = QA, subgroup.name = "Quality")

export_3panel_forest(m_vol, m_vol_mod, m_vol_qa, "figures/Volume_Panel.pdf", "Pooled mean", "Pooled mean", col_overall)
export_funnel(m_vol, "figures/Volume_Funnel.pdf", "Funnel Plot: Splenic Volume")


# 5. Splenic Artery Diameter
df_art <- subset(dat, !is.na(ArtDiam_M) & !is.na(ArtDiam_SE))
m_art <- metagen(TE = ArtDiam_M, seTE = ArtDiam_SE, studlab = Study, data = df_art, sm = "MD", random = TRUE, common = FALSE, method.tau = "REML")
m_art_mod <- update(m_art, subgroup = Modality, subgroup.name = "Modality")
m_art_qa  <- update(m_art, subgroup = QA, subgroup.name = "Quality")

export_3panel_forest(m_art, m_art_mod, m_art_qa, "figures/Artery_Diam_Panel.pdf", "Pooled mean", "Pooled mean", col_overall)
export_funnel(m_art, "figures/Artery_Diam_Funnel.pdf", "Funnel Plot: Artery Diameter")


# ==============================================================================
# SECTION B: TWO-GROUP CONTINUOUS VARIABLE
# ==============================================================================

# 6. Volume Differences (Male vs Female)
df_sex <- subset(dat, !is.na(VolSex_Male.N) & !is.na(VolSex_Female.N))
m_sex <- metacont(n.e = VolSex_Male.N, mean.e = VolSex_Male.M, sd.e = VolSex_Male.SD,
                  n.c = VolSex_Female.N, mean.c = VolSex_Female.M, sd.c = VolSex_Female.SD,
                  studlab = Study, data = df_sex, sm = "MD", random = TRUE, common = FALSE, method.tau = "REML")
m_sex_mod <- update(m_sex, subgroup = Modality, subgroup.name = "Modality")
m_sex_qa  <- update(m_sex, subgroup = QA, subgroup.name = "Quality")

export_3panel_forest(m_sex, m_sex_mod, m_sex_qa, "figures/Volume_Sex_Panel.pdf", "Mean Difference", "MD", col_overall)
export_funnel(m_sex, "figures/Volume_Sex_Funnel.pdf", "Funnel Plot: Volume (Male vs Female)")


# ==============================================================================
# SECTION C: PROPORTIONS (FREEMAN-TUKEY TRANSFORMATION)
# ==============================================================================

# 7. Two (2) Primary Branches
df_b2 <- subset(dat, !is.na(Branches_N) & !is.na(Branches_2))
m_b2 <- metaprop(event = Branches_2, n = Branches_N, studlab = Study, data = df_b2, sm = "PFT", random = TRUE, common = FALSE, method.tau = "REML")
m_b2_mod <- update(m_b2, subgroup = Modality, subgroup.name = "Modality")
m_b2_qa  <- update(m_b2, subgroup = QA, subgroup.name = "Quality")

export_3panel_forest(m_b2, m_b2_mod, m_b2_qa, "figures/Branches_2_Panel.pdf", "Proportion", "Proportion", col_branch2)
export_funnel(m_b2, "figures/Branches_2_Funnel.pdf", "Funnel Plot: Two Branches")


# 8. Three (3) Primary Branches
df_b3 <- subset(dat, !is.na(Branches_N) & !is.na(Branches_3))
m_b3 <- metaprop(event = Branches_3, n = Branches_N, studlab = Study, data = df_b3, sm = "PFT", random = TRUE, common = FALSE, method.tau = "REML")
m_b3_mod <- update(m_b3, subgroup = Modality, subgroup.name = "Modality")
m_b3_qa  <- update(m_b3, subgroup = QA, subgroup.name = "Quality")

export_3panel_forest(m_b3, m_b3_mod, m_b3_qa, "figures/Branches_3_Panel.pdf", "Proportion", "Proportion", col_overall)
export_funnel(m_b3, "figures/Branches_3_Funnel.pdf", "Funnel Plot: Three Branches")
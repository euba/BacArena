model <- readSBMLmod("data/ecoli_core.xml", bndCond = FALSE)
save(model, file="data/ecore_model.R")

model <- readSBMLmod("data/ecoli_iAF1260.xml", bndCond = FALSE)
save(model, file="data/Bcoli_model.R")

model <- readSBMLmod("data/barkeri_iAF692.xml", bndCond = FALSE)
save(model, file="data/barkeri_model.R")

model <- readSBMLmod("data/Cbeijerinckii_iCM925.xml", bndCond = FALSE)
save(model, file="data/clos_model.R")
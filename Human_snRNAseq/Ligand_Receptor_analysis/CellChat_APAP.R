library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

###-------------------------------------------------
# Cell Chat Analysis - Without projectData
###-------------------------------------------------



# Here we load a scRNA-seq data matrix and its associated cell meta data

apap <- readRDS("/Human/apap.RDS")
data.input = apap@assays$RNA@data # normalized data matrix
meta = apap@meta.data # a dataframe with rownames containing cell mata data
meta$cluster_anno_number <- factor(meta$cluster_anno_number, levels = c("Hep_0", "Hep_1", "Hep_2", "Hep_3", "Hep_4", "Hep_5", "Hep_6","Hep_7", "Hep_8",
                                                          "End_0", "End_1","End_2", "End_3", "End_4", "End_5", "End_6", 
                                                          "End_7", "End_8", "End_9", "End_10", "Cho_0", "Cho_1", "Cho_2", "Cho_3",
                                                          "Cho_4", "Cho_5", "Mes_0", "Mes_1", "Mes_2", "Mes_3", "Mes_4", "MP_0", "MP_1", 
                                                          "MP_2", "MP_3", "MP_4", "MP_5", "Plasma", "T cells", "ILCs", "B cells"))
# Create CellChat Object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cluster_anno_number")


# Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human

# set the used database in the object
cellchat@DB <- CellChatDB


# Preprocessing the expression data for cell-cell communication analysis

## subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)


# Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat, population.size = FALSE, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Extract the inferred cellular communication network as a data frame
df.net <- subsetCommunication(cellchat)

# Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

# Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

saveRDS(cellchat, "/Human/CellChat_apap_Results/CellChat.RDS")
write.csv(df.net, "/Human/CellChat_apap_Results/interactions.csv")

cp <- SeuratPipe:::discrete_col_pal[c(1, 2, 5, 7, 4, 6, 3, 9, 10, 8)]
groupSize <- as.numeric(table(cellchat@idents))


target.use <- "Hep_2"
top = 0.2

png(paste0("/Human/CellChat_apap_Results/Migratory_target_circle_plot_top_0.2.png"), width = 15, height = 20, units = "in", res = 500)
 netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, targets.use =  target.use, top = top, color.use = c(cp[1],cp[1],cp[1],cp[1],cp[1],cp[1],cp[1],cp[1],cp[1],
                                                                                                                                              cp[3],cp[3],cp[3],cp[3],cp[3],cp[3],cp[3],cp[3],cp[3],cp[3],cp[3],
                                                                                                                                              cp[2],cp[2],cp[2],cp[2],cp[2],cp[2],
                                                                                                                                              cp[4],cp[4],cp[4],cp[4],cp[4],
                                                                                                                                              cp[5],cp[5],cp[5],cp[5],cp[5],cp[5],
                                                                                                                                              cp[9],cp[6],cp[7],cp[8]))
dev.off()


net <- cellchat@net$weight
thresh <- stats::quantile(net, probs = 1 - top)
net[net < thresh] <- 0
df.net2 <- reshape2::melt(net, value.name = "value")
colnames(df.net2)[1:2] <- c("source", "target")
df.net2 <- subset(df.net2, target %in% target.use)
df.net2 <- df.net2[df.net2$value >0,]
df.net2


png(paste0("/Human/CellChat_apap_Results/Migratory_specific_bubble_plot_top_0.2.png"), width = 7, height = 7, units = "in", res = 500)
print(netVisual_bubble(cellchat, sources.use = df.net2$source, targets.use = target.use, signaling = c("TGFb","BMP")))
dev.off()



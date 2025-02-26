
load("./Rdata/mydata_DMP_DMR.Rdata")
myDMP = myDMP$Control_to_GDM
rm(myDMR)


tiff(file = "./Figure/chr22.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
myDMP %>% 
  as_tibble() %>% 
  mutate(group = ifelse(logFC > 0, "hyper-DMPs", "hypo-DMPs")) %>% 
  mutate(MAPINFO = as.integer(MAPINFO)) %>% 
  mutate(ids = MAPINFO/10000000) %>% 
  filter(CHR == 22) %>% 
  ggplot(aes(x = ids, y = logFC, color = group)) + 
  geom_point() +
  scale_x_continuous(breaks = seq(0,6,1), limits = c(0,5.1))+
  scale_y_continuous(limits = c(-0.4,0.4)) +
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5,size = 22),
        axis.text.x = element_text(hjust = 0.5,size = 18), 
        axis.text.y = element_text(hjust = 0.5,size = 18),
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        legend.text = element_text(size = 18),
        legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_blank()) +
  ylab("Δβ") +
  xlab("Chr22 (Mb, x10^6 bp)") +
  scale_color_cosmic()
dev.off()


tiff(file = "./Figure/chr1.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
myDMP %>% 
  as_tibble() %>% 
  mutate(group = ifelse(logFC > 0, "hyper-DMPs", "hypo-DMPs")) %>% 
  mutate(MAPINFO = as.integer(MAPINFO)) %>% 
  mutate(ids = MAPINFO/10000000) %>% 
  filter(CHR == 1) %>% 
  ggplot(aes(x = ids, y = logFC, color = group)) + 
  geom_point() +
  scale_x_continuous(breaks = seq(0,25,1), limits = c(0,25))+
  scale_y_continuous(limits = c(-0.4,0.4)) +
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5,size = 22),
        axis.text.x = element_text(hjust = 0.5,size = 18), 
        axis.text.y = element_text(hjust = 0.5,size = 18),
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        legend.text = element_text(size = 18),
        legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_blank()) +
  ylab("Δβ") +
  xlab("Chr1 (Mb, x10^6 bp)") +
  scale_color_cosmic()
dev.off()

tiff(file = "./Figure/chr2.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
myDMP %>% 
  as_tibble() %>% 
  mutate(group = ifelse(logFC > 0, "hyper-DMPs", "hypo-DMPs")) %>% 
  mutate(MAPINFO = as.integer(MAPINFO)) %>% 
  mutate(ids = MAPINFO/10000000) %>% 
  filter(CHR == 2) %>% 
  ggplot(aes(x = ids, y = logFC, color = group)) + 
  geom_point() +
  scale_x_continuous(breaks = seq(0,25,1), limits = c(0,25))+
  scale_y_continuous(limits = c(-0.4,0.4)) +
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5,size = 22),
        axis.text.x = element_text(hjust = 0.5,size = 18), 
        axis.text.y = element_text(hjust = 0.5,size = 18),
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        legend.text = element_text(size = 18),
        legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_blank()) +
  ylab("Δβ") +
  xlab("Chr2 (Mb, x10^6 bp)") +
  scale_color_cosmic()
dev.off()

tiff(file = "./Figure/chr3.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
myDMP %>% 
  as_tibble() %>% 
  mutate(group = ifelse(logFC > 0, "hyper-DMPs", "hypo-DMPs")) %>% 
  mutate(MAPINFO = as.integer(MAPINFO)) %>% 
  mutate(ids = MAPINFO/10000000) %>% 
  filter(CHR == 3) %>% 
  ggplot(aes(x = ids, y = logFC, color = group)) + 
  geom_point() +
  scale_x_continuous(breaks = seq(0,20,1), limits = c(0,20))+
  scale_y_continuous(limits = c(-0.4,0.4)) +
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5,size = 22),
        axis.text.x = element_text(hjust = 0.5,size = 18), 
        axis.text.y = element_text(hjust = 0.5,size = 18),
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        legend.text = element_text(size = 18),
        legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_blank()) +
  ylab("Δβ") +
  xlab("Chr3 (Mb, x10^6 bp)") +
  scale_color_cosmic()
dev.off()

tiff(file = "./Figure/chr4.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
myDMP %>% 
  as_tibble() %>% 
  mutate(group = ifelse(logFC > 0, "hyper-DMPs", "hypo-DMPs")) %>% 
  mutate(MAPINFO = as.integer(MAPINFO)) %>% 
  mutate(ids = MAPINFO/10000000) %>% 
  filter(CHR == 4) %>% 
  ggplot(aes(x = ids, y = logFC, color = group)) + 
  geom_point() +
  scale_x_continuous(breaks = seq(0,19,1), limits = c(0,19))+
  scale_y_continuous(limits = c(-0.4,0.4)) +
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5,size = 22),
        axis.text.x = element_text(hjust = 0.5,size = 18), 
        axis.text.y = element_text(hjust = 0.5,size = 18),
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        legend.text = element_text(size = 18),
        legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_blank()) +
  ylab("Δβ") +
  xlab("Chr4 (Mb, x10^6 bp)") +
  scale_color_cosmic()
dev.off()

tiff(file = "./Figure/chr5.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
myDMP %>% 
  as_tibble() %>% 
  mutate(group = ifelse(logFC > 0, "hyper-DMPs", "hypo-DMPs")) %>% 
  mutate(MAPINFO = as.integer(MAPINFO)) %>% 
  mutate(ids = MAPINFO/10000000) %>% 
  filter(CHR == 5) %>% 
  ggplot(aes(x = ids, y = logFC, color = group)) + 
  geom_point() +
  scale_x_continuous(breaks = seq(0,18,1), limits = c(0,18))+
  scale_y_continuous(limits = c(-0.4,0.4)) +
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5,size = 22),
        axis.text.x = element_text(hjust = 0.5,size = 18), 
        axis.text.y = element_text(hjust = 0.5,size = 18),
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        legend.text = element_text(size = 18),
        legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_blank()) +
  ylab("Δβ") +
  xlab("Chr5 (Mb, x10^6 bp)") +
  scale_color_cosmic()
dev.off()

tiff(file = "./Figure/chr6.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
myDMP %>% 
  as_tibble() %>% 
  mutate(group = ifelse(logFC > 0, "hyper-DMPs", "hypo-DMPs")) %>% 
  mutate(MAPINFO = as.integer(MAPINFO)) %>% 
  mutate(ids = MAPINFO/10000000) %>% 
  filter(CHR == 6) %>% 
  ggplot(aes(x = ids, y = logFC, color = group)) + 
  geom_point() +
  scale_x_continuous(breaks = seq(0,17,1), limits = c(0,17))+
  scale_y_continuous(limits = c(-0.4,0.4)) +
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5,size = 22),
        axis.text.x = element_text(hjust = 0.5,size = 18), 
        axis.text.y = element_text(hjust = 0.5,size = 18),
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        legend.text = element_text(size = 18),
        legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_blank()) +
  ylab("Δβ") +
  xlab("Chr6 (Mb, x10^6 bp)") +
  scale_color_cosmic()
dev.off()

tiff(file = "./Figure/chr7.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
myDMP %>% 
  as_tibble() %>% 
  mutate(group = ifelse(logFC > 0, "hyper-DMPs", "hypo-DMPs")) %>% 
  mutate(MAPINFO = as.integer(MAPINFO)) %>% 
  mutate(ids = MAPINFO/10000000) %>% 
  filter(CHR == 7) %>% 
  ggplot(aes(x = ids, y = logFC, color = group)) + 
  geom_point() +
  scale_x_continuous(breaks = seq(0,16,1), limits = c(0,16))+
  scale_y_continuous(limits = c(-0.4,0.4)) +
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5,size = 22),
        axis.text.x = element_text(hjust = 0.5,size = 18), 
        axis.text.y = element_text(hjust = 0.5,size = 18),
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        legend.text = element_text(size = 18),
        legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_blank()) +
  ylab("Δβ") +
  xlab("Chr7 (Mb, x10^6 bp)") +
  scale_color_cosmic()
dev.off()

tiff(file = "./Figure/chr8.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
myDMP %>% 
  as_tibble() %>% 
  mutate(group = ifelse(logFC > 0, "hyper-DMPs", "hypo-DMPs")) %>% 
  mutate(MAPINFO = as.integer(MAPINFO)) %>% 
  mutate(ids = MAPINFO/10000000) %>% 
  filter(CHR == 8) %>% 
  ggplot(aes(x = ids, y = logFC, color = group)) + 
  geom_point() +
  scale_x_continuous(breaks = seq(0,14.5,1), limits = c(0,14.5))+
  scale_y_continuous(limits = c(-0.4,0.4)) +
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5,size = 22),
        axis.text.x = element_text(hjust = 0.5,size = 18), 
        axis.text.y = element_text(hjust = 0.5,size = 18),
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        legend.text = element_text(size = 18),
        legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_blank()) +
  ylab("Δβ") +
  xlab("Chr8 (Mb, x10^6 bp)") +
  scale_color_cosmic()
dev.off()

tiff(file = "./Figure/chr9.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
myDMP %>% 
  as_tibble() %>% 
  mutate(group = ifelse(logFC > 0, "hyper-DMPs", "hypo-DMPs")) %>% 
  mutate(MAPINFO = as.integer(MAPINFO)) %>% 
  mutate(ids = MAPINFO/10000000) %>% 
  filter(CHR == 9) %>% 
  ggplot(aes(x = ids, y = logFC, color = group)) + 
  geom_point() +
  scale_x_continuous(breaks = seq(0,14,1), limits = c(0,14))+
  scale_y_continuous(limits = c(-0.4,0.4)) +
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5,size = 22),
        axis.text.x = element_text(hjust = 0.5,size = 18), 
        axis.text.y = element_text(hjust = 0.5,size = 18),
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        legend.text = element_text(size = 18),
        legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_blank()) +
  ylab("Δβ") +
  xlab("Chr9 (Mb, x10^6 bp)") +
  scale_color_cosmic()
dev.off()

tiff(file = "./Figure/chr10.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
myDMP %>% 
  as_tibble() %>% 
  mutate(group = ifelse(logFC > 0, "hyper-DMPs", "hypo-DMPs")) %>% 
  mutate(MAPINFO = as.integer(MAPINFO)) %>% 
  mutate(ids = MAPINFO/10000000) %>% 
  filter(CHR == 10) %>% 
  ggplot(aes(x = ids, y = logFC, color = group)) + 
  geom_point() +
  scale_x_continuous(breaks = seq(0,13.5,1), limits = c(0,13.5))+
  scale_y_continuous(limits = c(-0.4,0.4)) +
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5,size = 22),
        axis.text.x = element_text(hjust = 0.5,size = 18), 
        axis.text.y = element_text(hjust = 0.5,size = 18),
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        legend.text = element_text(size = 18),
        legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_blank()) +
  ylab("Δβ") +
  xlab("Chr10 (Mb, x10^6 bp)") +
  scale_color_cosmic()
dev.off()

tiff(file = "./Figure/chr11.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
myDMP %>% 
  as_tibble() %>% 
  mutate(group = ifelse(logFC > 0, "hyper-DMPs", "hypo-DMPs")) %>% 
  mutate(MAPINFO = as.integer(MAPINFO)) %>% 
  mutate(ids = MAPINFO/10000000) %>% 
  filter(CHR == 11) %>% 
  ggplot(aes(x = ids, y = logFC, color = group)) + 
  geom_point() +
  scale_x_continuous(breaks = seq(0,13.5,1), limits = c(0,13.5))+
  scale_y_continuous(limits = c(-0.4,0.4)) +
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5,size = 22),
        axis.text.x = element_text(hjust = 0.5,size = 18), 
        axis.text.y = element_text(hjust = 0.5,size = 18),
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        legend.text = element_text(size = 18),
        legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_blank()) +
  ylab("Δβ") +
  xlab("Chr11 (Mb, x10^6 bp)") +
  scale_color_cosmic()
dev.off()

tiff(file = "./Figure/chr12.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
myDMP %>% 
  as_tibble() %>% 
  mutate(group = ifelse(logFC > 0, "hyper-DMPs", "hypo-DMPs")) %>% 
  mutate(MAPINFO = as.integer(MAPINFO)) %>% 
  mutate(ids = MAPINFO/10000000) %>% 
  filter(CHR == 12) %>% 
  ggplot(aes(x = ids, y = logFC, color = group)) + 
  geom_point() +
  scale_x_continuous(breaks = seq(0,13.5,1), limits = c(0,13.5))+
  scale_y_continuous(limits = c(-0.4,0.4)) +
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5,size = 22),
        axis.text.x = element_text(hjust = 0.5,size = 18), 
        axis.text.y = element_text(hjust = 0.5,size = 18),
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        legend.text = element_text(size = 18),
        legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_blank()) +
  ylab("Δβ") +
  xlab("Chr12 (Mb, x10^6 bp)") +
  scale_color_cosmic()
dev.off()

tiff(file = "./Figure/chr13.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
myDMP %>% 
  as_tibble() %>% 
  mutate(group = ifelse(logFC > 0, "hyper-DMPs", "hypo-DMPs")) %>% 
  mutate(MAPINFO = as.integer(MAPINFO)) %>% 
  mutate(ids = MAPINFO/10000000) %>% 
  filter(CHR == 13) %>% 
  ggplot(aes(x = ids, y = logFC, color = group)) + 
  geom_point() +
  scale_x_continuous(breaks = seq(0,11.5,1), limits = c(0,11.5))+
  scale_y_continuous(limits = c(-0.4,0.4)) +
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5,size = 22),
        axis.text.x = element_text(hjust = 0.5,size = 18), 
        axis.text.y = element_text(hjust = 0.5,size = 18),
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        legend.text = element_text(size = 18),
        legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_blank()) +
  ylab("Δβ") +
  xlab("Chr13 (Mb, x10^6 bp)") +
  scale_color_cosmic()
dev.off()

tiff(file = "./Figure/chr14.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
myDMP %>% 
  as_tibble() %>% 
  mutate(group = ifelse(logFC > 0, "hyper-DMPs", "hypo-DMPs")) %>% 
  mutate(MAPINFO = as.integer(MAPINFO)) %>% 
  mutate(ids = MAPINFO/10000000) %>% 
  filter(CHR == 14) %>% 
  ggplot(aes(x = ids, y = logFC, color = group)) + 
  geom_point() +
  scale_x_continuous(breaks = seq(0,11,1), limits = c(0,11))+
  scale_y_continuous(limits = c(-0.4,0.4)) +
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5,size = 22),
        axis.text.x = element_text(hjust = 0.5,size = 18), 
        axis.text.y = element_text(hjust = 0.5,size = 18),
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        legend.text = element_text(size = 18),
        legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_blank()) +
  ylab("Δβ") +
  xlab("Chr14 (Mb, x10^6 bp)") +
  scale_color_cosmic()
dev.off()

tiff(file = "./Figure/chr15.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
myDMP %>% 
  as_tibble() %>% 
  mutate(group = ifelse(logFC > 0, "hyper-DMPs", "hypo-DMPs")) %>% 
  mutate(MAPINFO = as.integer(MAPINFO)) %>% 
  mutate(ids = MAPINFO/10000000) %>% 
  filter(CHR == 15) %>% 
  ggplot(aes(x = ids, y = logFC, color = group)) + 
  geom_point() +
  scale_x_continuous(breaks = seq(0,10,1), limits = c(0,10))+
  scale_y_continuous(limits = c(-0.4,0.4)) +
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5,size = 22),
        axis.text.x = element_text(hjust = 0.5,size = 18), 
        axis.text.y = element_text(hjust = 0.5,size = 18),
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        legend.text = element_text(size = 18),
        legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_blank()) +
  ylab("Δβ") +
  xlab("Chr15 (Mb, x10^6 bp)") +
  scale_color_cosmic()
dev.off()

tiff(file = "./Figure/chr16.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
myDMP %>% 
  as_tibble() %>% 
  mutate(group = ifelse(logFC > 0, "hyper-DMPs", "hypo-DMPs")) %>% 
  mutate(MAPINFO = as.integer(MAPINFO)) %>% 
  mutate(ids = MAPINFO/10000000) %>% 
  filter(CHR == 16) %>% 
  ggplot(aes(x = ids, y = logFC, color = group)) + 
  geom_point() +
  scale_x_continuous(breaks = seq(0,9,1), limits = c(0,9))+
  scale_y_continuous(limits = c(-0.4,0.4)) +
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5,size = 22),
        axis.text.x = element_text(hjust = 0.5,size = 18), 
        axis.text.y = element_text(hjust = 0.5,size = 18),
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        legend.text = element_text(size = 18),
        legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_blank()) +
  ylab("Δβ") +
  xlab("Chr16 (Mb, x10^6 bp)") +
  scale_color_cosmic()
dev.off()

tiff(file = "./Figure/chr17.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
myDMP %>% 
  as_tibble() %>% 
  mutate(group = ifelse(logFC > 0, "hyper-DMPs", "hypo-DMPs")) %>% 
  mutate(MAPINFO = as.integer(MAPINFO)) %>% 
  mutate(ids = MAPINFO/10000000) %>% 
  filter(CHR == 17) %>% 
  ggplot(aes(x = ids, y = logFC, color = group)) + 
  geom_point() +
  scale_x_continuous(breaks = seq(0,8,1), limits = c(0,8))+
  scale_y_continuous(limits = c(-0.4,0.4)) +
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5,size = 22),
        axis.text.x = element_text(hjust = 0.5,size = 18), 
        axis.text.y = element_text(hjust = 0.5,size = 18),
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        legend.text = element_text(size = 18),
        legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_blank()) +
  ylab("Δβ") +
  xlab("Chr17 (Mb, x10^6 bp)") +
  scale_color_cosmic()
dev.off()

tiff(file = "./Figure/chr18.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
myDMP %>% 
  as_tibble() %>% 
  mutate(group = ifelse(logFC > 0, "hyper-DMPs", "hypo-DMPs")) %>% 
  mutate(MAPINFO = as.integer(MAPINFO)) %>% 
  mutate(ids = MAPINFO/10000000) %>% 
  filter(CHR == 18) %>% 
  ggplot(aes(x = ids, y = logFC, color = group)) + 
  geom_point() +
  scale_x_continuous(breaks = seq(0,8,1), limits = c(0,8))+
  scale_y_continuous(limits = c(-0.4,0.4)) +
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5,size = 22),
        axis.text.x = element_text(hjust = 0.5,size = 18), 
        axis.text.y = element_text(hjust = 0.5,size = 18),
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        legend.text = element_text(size = 18),
        legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_blank()) +
  ylab("Δβ") +
  xlab("Chr18 (Mb, x10^6 bp)") +
  scale_color_cosmic()
dev.off()

tiff(file = "./Figure/chr19.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
myDMP %>% 
  as_tibble() %>% 
  mutate(group = ifelse(logFC > 0, "hyper-DMPs", "hypo-DMPs")) %>% 
  mutate(MAPINFO = as.integer(MAPINFO)) %>% 
  mutate(ids = MAPINFO/10000000) %>% 
  filter(CHR == 19) %>% 
  ggplot(aes(x = ids, y = logFC, color = group)) + 
  geom_point() +
  scale_x_continuous(breaks = seq(0,6,1), limits = c(0,6))+
  scale_y_continuous(limits = c(-0.4,0.4)) +
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5,size = 22),
        axis.text.x = element_text(hjust = 0.5,size = 18), 
        axis.text.y = element_text(hjust = 0.5,size = 18),
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        legend.text = element_text(size = 18),
        legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_blank()) +
  ylab("Δβ") +
  xlab("Chr19 (Mb, x10^6 bp)") +
  scale_color_cosmic()
dev.off()

tiff(file = "./Figure/chr20.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
myDMP %>% 
  as_tibble() %>% 
  mutate(group = ifelse(logFC > 0, "hyper-DMPs", "hypo-DMPs")) %>% 
  mutate(MAPINFO = as.integer(MAPINFO)) %>% 
  mutate(ids = MAPINFO/10000000) %>% 
  filter(CHR == 20) %>% 
  ggplot(aes(x = ids, y = logFC, color = group)) + 
  geom_point() +
  scale_x_continuous(breaks = seq(0,6.5,1), limits = c(0,6.5))+
  scale_y_continuous(limits = c(-0.4,0.4)) +
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5,size = 22),
        axis.text.x = element_text(hjust = 0.5,size = 18), 
        axis.text.y = element_text(hjust = 0.5,size = 18),
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        legend.text = element_text(size = 18),
        legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_blank()) +
  ylab("Δβ") +
  xlab("Chr20 (Mb, x10^6 bp)") +
  scale_color_cosmic()
dev.off()

tiff(file = "./Figure/chr21.tiff", 
     width = 10 * 300,
     height = 8 * 300, 
     res = 300, 
     compression = "lzw")
myDMP %>% 
  as_tibble() %>% 
  mutate(group = ifelse(logFC > 0, "hyper-DMPs", "hypo-DMPs")) %>% 
  mutate(MAPINFO = as.integer(MAPINFO)) %>% 
  mutate(ids = MAPINFO/10000000) %>% 
  filter(CHR == 21) %>% 
  ggplot(aes(x = ids, y = logFC, color = group)) + 
  geom_point() +
  scale_x_continuous(breaks = seq(0,5,1), limits = c(0,5))+
  scale_y_continuous(limits = c(-0.4,0.4)) +
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5,size = 22),
        axis.text.x = element_text(hjust = 0.5,size = 18), 
        axis.text.y = element_text(hjust = 0.5,size = 18),
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        legend.text = element_text(size = 18),
        legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_blank()) +
  ylab("Δβ") +
  xlab("Chr21 (Mb, x10^6 bp)") +
  scale_color_cosmic()
dev.off()
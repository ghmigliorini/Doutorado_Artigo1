### Análises utilizadas para gerar os resultados do artigo sobre aquecimento e decomposição

### Artigo publicado na revista Scientific Reports - link para o artigo https://www.nature.com/articles/s41598-020-77382-7




# importar pacotes --------------------------------------------------------

library(FD)
library(picante)
library(vegan)
library(ggplot2)
library(dplyr)
library(FactoMineR)
library(corrplot)
library(cowplot)
library(nlme)
library(multcompView)
library(gridExtra)


# PCA com todos os traits e espécies de detritos (Figura 2 do artigo) --------------------------

# importar matriz de atributos dos detritos
trait_spp <- read.table("C:/Users/Gustavo/OneDrive/ECOLOGIA/Pós-graduação/Doutorado/Projeto/Atributos detritos/trait_matrix.txt", h=T)

ncol(comp_sp)
nrow (trait_spp)


# PCA de todos os traits
library(grid)
library(ggrepel)

#Corrigindo os nomes nas linhas e colunas
sp_names<-c("Jacaranda", "I. subnuda", "Alchornea", "Pera", "M. glabra", "M. racemosa", "Andira", "Abarema", "Cupania", "Miconia", "I. edulis", "Lacistema")
rownames(trait_spp)<-NULL
rownames(trait_spp)<-sp_names
trait_names<-c("C", "N", "P", "Lignin", "Phenolics", "Tannins", "C:N", "N:P")
colnames(trait_spp)<-NULL
colnames(trait_spp)<-trait_names
trait_spp

# gerar PCA
pca_atributos<- rda(trait_spp, scale=TRUE)
scores_atributos <- scores(pca_atributos, choice=c(1,2))$sites
summary(pca_atributos)

windowsFonts(Times=windowsFont("TT Times New Roman"))

# criar dataframe com dados da PCA para gerar gráfico
smry <- summary(pca_atributos)
df1  <- data.frame(smry$sites[,1:2])
df2  <- data.frame(smry$species[,1:2])

lt1<-subset(df1, rownames(df1) %in% c("Cupania", "Pera", "Alchornea", "Miconia"))
lt2<-subset(df1, rownames(df1) %in% c("Lacistema", "M. glabra", "Jacaranda", "M. racemosa"))
lt3<-subset(df1, rownames(df1) %in% c("Lacistema", "Jacaranda", "Abarema", "M. racemosa"))
lt4<-subset(df1, rownames(df1) %in% c("Cupania", "Jacaranda", "I. subnuda", "Abarema"))
lt5<-subset(df1, rownames(df1) %in% c("I. subnuda", "Andira", "Abarema", "I. edulis"))

lt1<-mutate(lt1, LT="LT1")
lt2<-mutate(lt2, LT="LT2")
lt3<-mutate(lt3, LT="LT3")
lt4<-mutate(lt4, LT="LT4")
lt5<-mutate(lt5, LT="LT5")

lt_ellipse<-bind_rows(lt1,lt2,lt3,lt4,lt5)

hulls <- lt_ellipse %>%
  group_by(LT) %>%
  slice(chull(PC1, PC2))

# gráfico PCA

ggplot(df1, aes(x=PC1, y=PC2)) +
  geom_text(aes(label=rownames(df1)),size=4.5, color="grey40", fontface=c("bold.italic"), family="Times") +
  geom_hline(yintercept=0, linetype="dotted", alpha=0.3) +
  geom_vline(xintercept=0, linetype="dotted", alpha=0.3) +
  scale_x_continuous(limits = c(-2, 2)) +
  scale_y_continuous(limits = c(-2, 2)) +
  geom_segment(data=df2, aes(x=0, xend=PC1, y=0, yend=PC2),
               color="red", arrow=arrow(length=unit(0.01,"npc")), size=0.2) +
  geom_text_repel(data=df2,
                  aes(x=PC1,y=PC2,label=rownames(df2)),
                  color="black", size=4, fontface="bold") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x="PC1 (53.5%)", y="PC2 (19.2%)") +
  theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size = 9),
        axis.title.y = element_text(size = 10), axis.text.y = element_text(size = 9)) +
  geom_polygon(data = hulls, aes(x =  PC1, y = PC2, fill = LT, color = LT), size=0.1, alpha=0.1) +
  geom_point(size=0.5) +
  theme(legend.title = element_blank(), legend.position = c(0.08, 0.87)) +
  scale_color_manual(values=c("green", "blue", "red", "yellow", "black")) +
  scale_fill_manual(values=c("green", "blue", "red", "yellow", "black"))

ggsave("pca_traits.svg", width = 12, height = 10, units = "cm", dpi=600)


### A partir da PCA acima foi possível gerar 5 tratamentos contendo 4 espécies em cada
### Esses tratamentos foram criados de acordo com a posição de cada espécie na PCA, ou seja a proximidade com cada trait
### Seguindo um gradiente de qualidade LT1 -> LT5 - Figura 2 do artigo


# Diversidade funcional ---------------------------------------------------

# Padronizando os traits com Z-score
trait_spp <- read.table("C:/Users/Gustavo/OneDrive/ECOLOGIA/Pós-graduação/Doutorado/Projeto/Atributos detritos/trait_matrix.txt", h=T)

traits_std<-decostand(trait_spp, "standardize")

comp_sp_new<-comp_sp[-c(2:5,7:10,12:15,17:20,22:25),]
lt_names<-as.vector(rbind(rep(c("LT1", "LT2", "LT3", "LT4", "LT5"))))
rownames(comp_sp_new)<-NULL
rownames(comp_sp_new)<-lt_names
com_cwm<-as.matrix(comp_sp_new)

rao_q <- dbFD(traits_std, com_cwm)$RaoQ

C <- dbFD(traits_std[,1, drop=FALSE], com_cwm)$RaoQ
N <- dbFD(traits_std[,2, drop=FALSE], com_cwm)$RaoQ
P <- dbFD(traits_std[,3, drop=FALSE], com_cwm)$RaoQ
Lignin <- dbFD(traits_std[,4, drop=FALSE], com_cwm)$RaoQ
Phenolics <- dbFD(traits_std[,5, drop=FALSE], com_cwm)$RaoQ
Tannins <- dbFD(traits_std[,6, drop=FALSE], com_cwm)$RaoQ
C.N <- dbFD(traits_std[,7, drop=FALSE], com_cwm)$RaoQ
N.P <- dbFD(traits_std[,8, drop=FALSE], com_cwm)$RaoQ

rao_df<-data.frame(C, N, P, Lignin, Phenolics, Tannins, C.N, N.P)

rao_pc<-PCA(rao_df, graph=FALSE, scale=TRUE)
rao_pc$eig
rao_pc$var$coord
rao_pc$ind$coord
plot(rao_pc, choix = "ind")
dimdesc(rao_pc, axes = 1:2)


##PCA FD
#rotation
rotation_PCA<-t(apply(rao_pc$var$coord, 1, function(x) {x/sqrt(rao_pc$eig[,1])}))
labels<-c("C", "N", "P", "Lignin", "Phenolics", "Tannins", "C:N", "N:P")
rownames(rotation_PCA)<-NULL
rownames(rotation_PCA)<-labels

# função para criar biplot
PCbiplot <- function(rao_pc, x="Dim.1", y="Dim.2") {
  data <- data.frame(obsnames=row.names(rao_pc$ind$coord), rao_pc$ind$coord)
  plot <- ggplot(data, aes_string(x=x, y=y)) +
    geom_text(color=1, size=4, aes(label=obsnames)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(axis.title.x=element_text(size=8), axis.text.x=element_text(size=6), axis.title.y=element_text(size=8),
          axis.text.y=element_text(size=6)) +
    labs(x="FD1 (50.95%)", y="FD2 (31.48%)") +
    geom_vline(xintercept=c(-0,0), linetype="dotted") +
    geom_hline(yintercept = c(-0,0), linetype="dotted")
  datapc <- data.frame(varnames=rownames(rotation_PCA), rotation_PCA)
  mult <- min(
    (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
    (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
  )
  datapc <- transform(datapc,
                      v1 = .7 * mult * (get(x)),
                      v2 = .7 * mult * (get(y))
  )
  plot <- plot + geom_text(data=datapc, aes(x=v1, y=v2, label=varnames), size = 3, vjust=-0.6,  color="grey20")
  plot <- plot + geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), arrow = arrow(length=unit(0.1,"cm")), alpha=1, color="black")
  plot
}

# aplicar a função aos dados
p1<-PCbiplot(rao_pc)
p1

# salvar o gráfico
ggsave("pca_fd.png", width=8, height=7, units="cm", dpi=600)
ggsave("pca_fd.svg", width=8, height=7, units="cm", dpi=600)


# Qualidade média ---------------------------------------------------------

mean_traits<-functcomp(traits_std, com_cwm)

mean_pc<-PCA(mean_traits, graph=FALSE, scale=TRUE)
mean_pc$eig
mean_pc$var$coord
mean_pc$ind$coord
plot(mean_pc, choix = "ind")
dimdesc(mean_pc, axes = 1:2)

##PCA quality
#rotation
rotation_PCA<-t(apply(mean_pc$var$coord, 1, function(x) {x/sqrt(mean_pc$eig[,1])}))
labels<-c("C", "N", "P", "Lignin", "Phenolics", "Tannins", "C:N", "N:P")
rownames(rotation_PCA)<-NULL
rownames(rotation_PCA)<-labels

# função para criar biplot
PCbiplot <- function(mean_pc, x="Dim.1", y="Dim.2") {
  data <- data.frame(obsnames=row.names(mean_pc$ind$coord), mean_pc$ind$coord)
  plot <- ggplot(data, aes_string(x=x, y=y)) +
    geom_text(color=1, size=4, aes(label=obsnames)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(axis.title.x=element_text(size=8), axis.text.x=element_text(size=6), axis.title.y=element_text(size=8),
          axis.text.y=element_text(size=6)) +
    labs(x="Quality1 (74.42%)", y="Quality2 (19.83%)") +
    geom_vline(xintercept=c(-0,0), linetype="dotted") +
    geom_hline(yintercept = c(-0,0), linetype="dotted")
  datapc <- data.frame(varnames=rownames(rotation_PCA), rotation_PCA)
  mult <- min(
    (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
    (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
  )
  datapc <- transform(datapc,
                      v1 = .7 * mult * (get(x)),
                      v2 = .7 * mult * (get(y))
  )
  plot <- plot + geom_text(data=datapc, aes(x=v1, y=v2, label=varnames), size = 3, vjust=-0.1,  color="grey20")
  plot <- plot + geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), arrow = arrow(length=unit(0.1,"cm")), alpha=1, color="black")
  plot
}

# aplicar a função aos dados
p2<-PCbiplot(mean_pc)
p2

# salvar o gráfico
ggsave("pca_quality.png", width=8, height=7, units="cm", dpi=600)
ggsave("pca_quality.svg", width=8, height=7, units="cm", dpi=600)


# juntar os dois plots PCA de FD e Quality
cowplot::plot_grid(p1, p2, labels = c("A", "B"), align = "h", vjust = 8)


# correlações dos traits com eixos para gerar uma figura; figura 2 artigo

fd1_cor<-as.data.frame(rao_pc$var$cor[,1])
colnames(fd1_cor)<-NULL
fd1_cor_name<-as.vector(rbind("FD1"))
colnames(fd1_cor)<-fd1_cor_name
fd1_cor

mean_cor<-as.data.frame(mean_pc$var$cor[,1])
colnames(mean_cor)<-NULL
mean_cor_name<-as.vector(rbind("Quality1"))
colnames(mean_cor)<-mean_cor_name
mean_cor

# Juntando os dados de correlações dos eixos de qualidade e diversidade funcional
cor_eixos<-cbind(fd1_cor, mean_cor)
rownames(cor_eixos)[rownames(cor_eixos) == "carbono"] <- "C"
rownames(cor_eixos)[rownames(cor_eixos) == "nitrogenio"] <- "N"
rownames(cor_eixos)[rownames(cor_eixos) == "fosforo"] <- "P"
rownames(cor_eixos)[rownames(cor_eixos) == "lignina"] <- "Lignin"
rownames(cor_eixos)[rownames(cor_eixos) == "fenois"] <- "Phenolics"
rownames(cor_eixos)[rownames(cor_eixos) == "taninos"] <- "Tannins"
rownames(cor_eixos)[rownames(cor_eixos) == "C.N"] <- "C:N"
rownames(cor_eixos)[rownames(cor_eixos) == "N.P"] <- "N:P"
cor_eixos

cor_eixos<-as.matrix(cor_eixos)

pdf(file="corrplot_traits_eixos.pdf")
tiff("corrplot_traits_eixos.tiff", width=10, height=8, units="cm", res=300)
eps("corrplot_traits_eixos.eps")
corrplot(cor_eixos, cl.pos = "b", cl.length = 3, tl.col= "black", tl.cex = 2, cl.cex = 1.5)

dev.off()


# correlação FD - CWM
library(dplyr)

fd_traits<-rename(rao_df, FD_C=C, FD_N=N, FD_P=P, FD_lignin=Lignin, FD_phenolics=Phenolics, FD_tannins=Tannins, "FD_C:N"=C.N, "FD_N:P"=N.P)
cwm_traits<-rename(mean_traits, CWM_C=carbono, CWM_N=nitrogenio, CWM_P=fosforo, CWM_lignin=lignina, CWM_phenolics=fenois, CWM_tannins=taninos, "CWM_C:N"=C.N, "CWM_N:P"=N.P)

fd_cwm_cor<-bind_cols(fd_traits, cwm_traits)
mt<-cor(fd_cwm_cor)
corrplot(mt, method = "number")

svg("corrplot_cwm_fd.svg", width=8, height=8)
corrplot(mt, method = "number")

dev.off()


# plots dos indices de qualidade e fd em relação aos tratamentos - fig.1 do artigo

fine<-read.table("C:/Users/Gustavo/OneDrive/ECOLOGIA/Pós-graduação/Doutorado/Projeto/Campo doutorado 2015/Planilhas R doc/decay.fine3.txt", h=T)

p1 <- ggplot(fine, aes(x = treatment, y = quality1, group = 1)) +
  geom_smooth(method = "loess", span = 0.8, se = FALSE, color = "grey50", size = 2) +
  geom_point(size = 5, shape = 21) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_discrete(breaks = c("trat1", "trat2", "trat3", "trat4", "trat5"), labels = c("LT1", "LT2", "LT3", "LT4", "LT5")) +
  labs(y = "Quality")+
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size=16),
        axis.title.y = element_text(size=18), axis.text.y = element_text(size = 16))

ggsave("quality_LT.png", width=8, height=7, units="cm", dpi=600)


p2 <- ggplot(fine, aes(x = treatment, y = FD1, group = 1)) +
  geom_smooth(method = "loess", span = 0.8, se = FALSE, color="grey50", size = 2) +
  geom_point(size = 5, shape = 21) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_discrete(breaks = c("trat1", "trat2", "trat3", "trat4", "trat5"), labels = c("LT1", "LT2", "LT3", "LT4", "LT5")) +
  labs(y = "Functional Diversity") +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 16))

ggsave("FD_LT.png", width=8, height=7, units="cm", dpi=600)


#--------------------------------------------------------------------------------------------------------------------------------------------------


# análises decomposição ---------------------------------------------------

library(nlme)


# coarse x fine
fine_coarse<-read.table("C:/Users/Gustavo/OneDrive/ECOLOGIA/Pós-graduação/Doutorado/Projeto/Campo doutorado 2015/Planilhas R doc/fine_coarse.txt", h=T)
head(fine_coarse)
fine_coarse$block<-as.factor(fine_coarse$block)

# check model
m<-lm(k~treatment*temp.value*mesh, data = fine.coarse)
plot(m1)

# model
model1<-lme(log(k)~treatment+temp.value+mesh, random=~1|block, method="REML", data=fine.coarse)
anova(model1)

#pairwise test
lsmeans::lsmeans(model1, list(pairwise~treatment), adjust="tukey")

# gerar letras das diferenças para inserir sobre as barras no gráfico
model = lm(log(k)~treatment, data = fine_coarse)
mod.ANOVA=aov(model)
TUKEY <- TukeyHSD(mod.ANOVA)

generate_label_df <- function(TUKEY, treatment){
  Tukey.levels <- TUKEY$treatment[,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  Tukey.labels$treatment=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
  return(Tukey.labels)
}

# Apply the function on my dataset
LABELS <- generate_label_df(TUKEY , "treatment")

names(LABELS) <- c('Letters','treatment')
yvalue <- fine_coarse |>
  group_by(treatment) |>
  summarise_all(mean)

final <- merge(LABELS, yvalue)


# gráfico - decomposição x malhas (material suplementar)
p1 <- ggplot(fine_coarse, aes(mesh, log(k))) +
  theme_bw() +
  geom_boxplot(outlier.size=0.5) +
  geom_jitter(aes(color=mesh), width = 0.2, alpha = 0.3, size=1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values = c("black", "black"))+
  scale_x_discrete(breaks=c("coarse", "fine"), labels=c("Coarse-mesh", "Fine-mesh")) +
  labs(y="Decomposition (k, log)") +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 11),
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 11)) +
  theme(legend.position = "none")

ggsave("mesh_dec.png", width=8, height=7, units="cm", dpi=600)

# gráfico decomposição x tratamentos - figura 3
p2 <- ggplot(fine_coarse, aes(treatment, log(k))) +
  theme_bw() +
  geom_boxplot(colour = "black", outlier.shape = NA) +
  geom_jitter(color="black", shape=21, fill=I("seashell4"), width = 0.2, alpha = 0.3, size=1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_discrete(breaks=c("trat1", "trat2", "trat3", "trat4", "trat5"),
                   labels=c("LT1", "LT2", "LT3", "LT4", "LT5")) +
  labs(tag= "(a)", y=expression(Decomposition~(italic("k")~d^-1~","~Ln)), x="Litter treatments") +
  theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 11),
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 11)) +
  theme(legend.position = "none") +
  geom_text(data = final, aes(x = treatment, y = log(k), label = Letters),vjust=-8.5,hjust=0.5, size=3) +
  scale_y_continuous(limits = c(-7,-2.5))

ggsave("decomp_LT.png", width=8, height=8, units="cm", dpi=600)


# gráfico decomposição x temperatura - figura 3
p3 <- ggplot(fine_coarse, aes(temp.value, log(k)))+
  theme_bw() +
  geom_jitter(aes(group = temp.value), alpha = 0.3, size=2) +
  geom_smooth(method="lm", se=TRUE, size=0.5, fill="blue", alpha=0.2) +
  scale_x_continuous(breaks = c(20,22,24,26)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(tag= "(b)", x= "Temperature (°C)", y="Decomposition (k, log)") +
  theme(axis.title.x = element_text(size=12), axis.text.x = element_text(size = 11),
        axis.title.y = element_blank(), axis.text.y = element_text(size=11)) +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(-7,-2.5))

ggsave("decomp_temp_log.png", width=8, height=8, units="cm", dpi=600)

grid.arrange(p2, p3, ncol=2)
final_plot<-arrangeGrob(p2, p3, ncol=2)
ggsave("decomp_fig.3.jpeg", final_plot, width=16, height=8, units="cm", dpi=600)
ggsave("decomp_fig.3.svg", final_plot, width=16, height=8, units="cm", dpi=300)
ggsave("decomp_fig.3.tiff", final_plot, width=16, height=8, units="cm", dpi=300)


# Microbial decomposition -------------------------------------------------

fine<-read.table("C:/Users/Gustavo/OneDrive/ECOLOGIA/Pós-graduação/Doutorado/Projeto/Campo doutorado 2015/Planilhas R doc/decay.fine3.txt", h=T)
head(fine)
fine$block<-as.factor(fine$block)

## check model
m<-lm(k~(FD1+quality1)*temp.value, data = fine)
plot(m)

#model
model2<-lme(log(k) ~ FD1+quality1 + temp.value, random=~1|block, data=fine)
summary(model2)

# gráfico decomposição - figura 4

p1<-ggplot(fine.new2, aes(quality1, log(k)))+
  theme_bw() +
  geom_jitter(aes(x=quality1, y=log(k), shape=treatment), color="black",
              fill=I("seashell4"), size=3, position = position_jitter(w=0.15, h=0.15)) +
  scale_shape_manual(values=c(21, 22, 23, 24, 25)) +
  scale_y_continuous(limits=c(-6,-3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none")+
  labs(tag= "(a)", y=expression(Fine-mesh~(italic("k")~d^-1~","~Ln))) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 11),
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 11))

ggsave("fine_quality1.png", width=8, height=7, units="cm", dpi=600)

p2<-ggplot(fine.new2, aes(FD1, log(k)))+
  theme_bw() +
  geom_jitter(aes(x=FD1, y=log(k), shape=treatment), color="black",
              fill=I("seashell4"), size=3, position = position_jitter(w=0.15, h=0.15)) +
  geom_smooth(method = "lm", se=TRUE, fill="blue", alpha=0.2) +
  scale_shape_manual(values=c(21, 22, 23, 24, 25))  +
  scale_y_continuous(limits=c(-6,-3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none")+
  labs(tag= "(b)", y=expression(Fine-mesh~(italic("k")~","~Ln))) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 11),
        axis.title.y = element_blank(), axis.text.y = element_text(size = 11))

ggsave("fine_FD1.png", width=8, height=7, units="cm", dpi=600)


# medir correlações
cor.temp.fine<-cor.test(log(fine$k), fine$temp.value, method = "pearson")
cor.temp.fine

cor.fd.fine<-cor.test(log(fine$k), fine$FD1, method = "pearson")
cor.fd.fine

# avaliar modelos
m1<-lme(log(k)~FD_C+FD_N+quality1+temp.value, random=~1|block, method = "REML", data=fine)
anova(m1)
summary(m1)

m2<-lme(log(k)~FD_C+FD_CN+quality1+temp.value, random=~1|block, method = "REML", data=fine)
anova(m2)
summary(m2)

anova(m1,m2)

m1<-lme(log(k)~quality1+FD1+temp.value, random=~1|block, method = "ML", data=fine)
anova(m1)

m2<-lme(log(k)~quality1+FD_C+FD_N+temp.value, random=~1|block, method = "ML", data=fine)
anova(m2)

anova(m1,m2)



# Total decomposition -----------------------------------------------------

coarse<-read.table("C:/Users/Gustavo/OneDrive/ECOLOGIA/Pós-graduação/Doutorado/Projeto/Campo doutorado 2015/Planilhas R doc/decay_coarse4.txt", h=T)
head(coarse)
coarse$block<-as.factor(coarse$block)

coarse.new<-coarse[, c(1:14)]
head(coarse.new)

coarse.new$block <- as.factor(coarse.new$block)
coarse.new$brom.id <- as.factor(coarse.new$brom.id)
coarse.new$temp.cat <- as.factor(coarse.new$temp.cat)
coarse.new$treatment <- as.factor(coarse.new$treatment)

coarse.new2<-inner_join(coarse.new, FD_new, by="treatment")
coarse.new2<-inner_join(coarse.new2, cwm_new, by="treatment")


## exploração da distribuição dos dados
m<-lm(k~(quality1+FD1)*temp.value*(Richness+Abundance), data = coarse)
plot(m)

#model
m1<-lme(log(k)~quality1+FD1+temp.value+Richness+Abundance, random=~1|block, data=coarse)
summary(m1)


# gráfico
p3<-ggplot(coarse, aes(quality1, log(k)))+
  theme_bw() +
  geom_jitter(aes(x=quality1, y=log(k), shape=treatment), color="black",
              fill=I("seashell4"), size=3, position = position_jitter(w=0.15, h=0.15)) +
  scale_shape_manual(values=c(21, 22, 23, 24, 25)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none")+
  labs(tag= "(a)", y=expression(Coarse-mesh~(italic("k")~d^-1~","~Ln))) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 11),
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 11))

p4<-ggplot(coarse, aes(FD1, log(k))) +
  theme_bw() +
  geom_jitter(aes(x=FD1, y=log(k), shape=treatment), color="black",
              fill=I("seashell4"), size=3, position = position_jitter(w=0.15, h=0.15)) +
  geom_smooth(method = "lm", se=TRUE, fill="blue", alpha=0.2) +
  scale_shape_manual(values=c(21, 22, 23, 24, 25)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  labs(tag= "(a)", y=expression(Coarse-mesh~(italic("k")~d^-1~","~Ln))) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 11),
        axis.title.y = element_blank(), axis.text.y = element_text(size = 11))


grid.arrange(p1, p2, p3, p4, ncol=2)
final_plot<-arrangeGrob(p1, p2, p3, p4, ncol=2)
ggsave("decomp_fig.4.jpeg", final_plot, width=16, height=14, units="cm", dpi=600)
ggsave("decomp_fig.4.svg", final_plot, width=16, height=14, units="cm", dpi=300)
ggsave("decomp_fig.4.tiff", final_plot, width=16, height=14, units="cm", dpi=300)



#plot decomposição x traits de FD - material suplementar

p1<-ggplot(fine, aes(FD_C, log(k))) +
  theme_bw() +
  geom_jitter(aes(x=FD_C, y=log(k), shape=treatment),
              color="darkgreen", fill="darkgreen", width = 0.05, alpha = 0.3, size=2) +
  geom_smooth(method = "lm", se=TRUE, color="darkgreen", size=0.2) +
  scale_shape_manual(values=c(21, 22, 23, 24, 25)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  scale_y_continuous(limits=c(-6,-3)) +
  labs(tag = "(a)", y=expression(Fine-mesh~(italic("k")~d^-1~","~Ln))) +
  theme(axis.title.x = element_text(size=12), axis.text.x = element_text(size = 11),
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 11))

p2<-ggplot(fine, aes(FD_N, log(k))) +
  theme_bw() +
  geom_jitter(aes(x=FD_N, y=log(k), shape=treatment),
              color="darkgreen", fill="darkgreen", width = 0.05, alpha = 0.3, size=2) +
  geom_smooth(method = "lm", se=TRUE, color="darkgreen", size=0.2) +
  scale_shape_manual(values=c(21, 22, 23, 24, 25)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  scale_y_continuous(limits=c(-6,-3)) +
  labs(tag = "(b)", y=expression(Fine-mesh~(italic("k")~d^-1~","~Ln))) +
  theme(axis.title.x = element_text(size=12), axis.text.x = element_text(size = 11),
        axis.title.y = element_blank(), axis.text.y = element_text(size = 11))

p3<-ggplot(coarse, aes(FD_C, log(k))) +
  theme_bw() +
  geom_jitter(aes(x=FD_C, y=log(k), shape=treatment),
              color="darkgreen", fill="darkgreen", width = 0.05, alpha = 0.3, size=2) +
  geom_smooth(method = "lm", se=TRUE, color="darkgreen", size=0.2) +
  scale_shape_manual(values=c(21, 22, 23, 24, 25)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none")+
  labs(tag = "(c)", y=expression(Coarse-mesh~(italic("k")~d^-1~","~Ln))) +
  theme(axis.title.x = element_text(size=12), axis.text.x = element_text(size = 11),
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 11))

p4<-ggplot(coarse, aes(FD_N, log(k))) +
  theme_bw() +
  geom_jitter(aes(x=FD_N, y=log(k), shape=treatment),
              color="darkgreen", fill="darkgreen", width = 0.05, alpha = 0.3, size=2) +
  geom_smooth(method = "lm", se=TRUE, color="darkgreen", size=0.2) +
  scale_shape_manual(values=c(21, 22, 23, 24, 25)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  labs(tag = "(d)",y="Coarse mesh (k, log)") +
  theme(axis.title.x = element_text(size=12), axis.text.x = element_text(size = 11),
        axis.title.y = element_blank(), axis.text.y = element_text(size = 11))

grid.arrange(p1, p2, p3, p4, ncol=2)
final_plot<-arrangeGrob(p1, p2, p4, p5, ncol=2)
ggsave("decomp_fig.S3.jpg", final_plot, width=16, height=14, units="cm", dpi=300)
ggsave("decomp_fig.S3.svg", final_plot, width=16, height=14, units="cm", dpi=300)

library(readr)
library(lattice)


filename = "~/git/COSYINF/test/fort.10"
TSS2d <- data.table(read_csv(filename, skip=2))
attr(TSS2d,"Parameter") <- read_csv(filename, n_max=1) -> para

TSS2d[,Dim:=2]; TSS3d[,`:=`(Dim=3,l=NULL)]

TSSdata <- rbind(TSS2d,TSS3d)

D0 <- TSSdata[y==0&ERG==min(abs(ERG))]

D0[,.(MEAN.Wx = mean(Wx), SD.Wx = sd(Wx), RSD.Wx = sd(Wx)/mean(Wx)), by=c("Dim","muTILT","x")] %>%
  ggplot(aes(x, MEAN.Wx,col=as.factor(Dim))) + geom_point(show.legend = FALSE) + 
  # geom_linerange(aes(ymin=MEAN.Wx-SD.Wx,ymax=MEAN.Wx+SD.Wx), show.legend = FALSE) +
  facet_grid(Dim~.,scales = "free_y") + theme(legend.position="top")


D0 %>%ggplot(aes(Wx)) + geom_histogram() +
  geom_violin(aes(xintercept=mean(Wx)),col="red") + 
  facet_grid(.~Dim, scales="free_x")

wireframe(Wx ~ x*y, data = TSSdata[ERG==min(abs(ERG))&Dim==3],
          xlab = "X coordinate", ylab = "Y coordinate",
          main = "Wx",
          drape = TRUE,
          colorkey = TRUE,
          screen = list(z = 15, x = -90)
)

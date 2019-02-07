
## Spatial lag correlation model

## Using spdep
library(spdep)

## make a neighbour list
nb <- knn2nb(knearneigh(gps))
nb.list <- nb2listw(nb)

## test for residual autocorrelation
lm.morantest(m1, nb.list) ## signficiant autocorrelation

## Lagrange Multiplier test statistic to determine optimal spatial error or spatial lag models
lagrange <- lm.LMtests(m1, nb.list, test=c("LMerr","RLMerr","LMlag","RLMlag","SARMA")) #Sarma highest

## Spatial lag regression
m2 <- lagsarlm(Shrub.density ~ RDM.2013, data=landscape, nb.list)




### Using ape package
shrub.dist <- as.matrix(dist(cbind(landscape$x,landscape$y))) ## calculate distances between shrubs
shrub.dist.inv <- 1/shrub.dist ## inverse distance
diag(shrub.dist.inv) <- 0 ## diagonal zero - same shrub to same shrub
shrub.dist.inv[is.infinite(shrub.dist.inv)] <- 0 ## replace infinites with zero 

## Conduct Moran's I to test spatial autocorrelation
Moran.I(m1$residuals, shrub.dist.inv)  # P < 0.001 significant autocorrelation


m1 <- lm(Area~ RDM.2013, data=landscape) 


## coords
gps <- landscape
coordinates(gps) <- ~x+y
proj4string(gps) <- CRS("+proj=utm +zone=11 +datum=WGS84")

## add buffer
ca <- rgeos::gBuffer(gps, width=50)

library(dismo)
v <- voronoi(gps)
plot(v)

vca <- intersect(v, ca)
spplot(vca, 'RDM.2013', col.regions=rev(get_col_regions()))
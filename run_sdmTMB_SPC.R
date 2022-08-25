library(sdmTMB)
library(dplyr)
library(ggplot2)
library(rgdal)
library(colorRamps)
library(raster)
library(sf)
library(readr)
library(ggthemes)
library(tidyr)

rm(list = ls())

select = dplyr::select

islands = c("Kauai", #1
            "Lehua", #2
            "Niihau", #3
            "Kaula", #4
            "Oahu", #5
            "Molokai", #6
            "Maui", #7
            "Lanai", #8
            "Molokini", #9
            "Kahoolawe", #10
            "Hawaii")#[7:11]

load("data/clean_df.RData")

zone <- (floor((df$lon[1] + 180)/6) %% 60) + 1
xy_utm = as.data.frame(cbind(utm = project(as.matrix(df[, c("lon", "lat")]), paste0("+proj=utm +units=km +zone=", zone))))
colnames(xy_utm) = c("X", "Y")
df = cbind(df, xy_utm)
plot(xy_utm, pch = ".", bty = 'n')
rm(xy_utm)

# Read in Island Boundaries
load('data/MHI_islands_shp.RData')
crs(ISL_bounds) = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
ISL_this = ISL_bounds[which(ISL_bounds$ISLAND %in% toupper(islands)),]
ISL_this_utm = spTransform(ISL_this,CRS(paste0("+proj=utm +units=km +zone=",zone)))
ISL_this_sf = st_transform(st_as_sf(ISL_this), crs = paste0("+proj=utm +units=km +zone=",zone))

rea_spde <- make_mesh(df, c("X", "Y"), n_knots = 100, type = "cutoff_search", seed = 2022) # search
# rea_spde <- make_mesh(df, c("X", "Y"), cutoff  = 15, type = "cutoff") # predefined

#build barrier to mesh
rea_spde_coast = add_barrier_mesh(rea_spde , ISL_this_sf)

plot(rea_spde_coast$mesh, asp = 1, main = ""); axis(1); axis(2)
plot(ISL_this_utm, add = TRUE)
points(rea_spde_coast$loc_xy,col = "green", pch = ".", cex = 5)
bar_i = rea_spde_coast$barrier_triangles
norm_i = rea_spde_coast$normal_triangles
points(rea_spde_coast$spde$mesh$loc[,1], rea_spde_coast$spde$mesh$loc[,2], pch = ".", col = "black")
points(rea_spde_coast$mesh_sf$V1[bar_i], rea_spde_coast$mesh_sf$V2[bar_i], col = "red", pch = 20, cex = 0.5)
points(rea_spde_coast$mesh_sf$V1[norm_i], rea_spde_coast$mesh_sf$V2[norm_i], col = "blue", pch = 20, cex = 0.5)

rm(ISL_bounds, ISL_this, ISL_this_sf, ISL_this_utm)

fit <- sdmTMB(
  
  data = df, 
  formula = response ~ as.factor(year) + s(depth),
  silent = F, 
  time = "year",
  mesh = rea_spde_coast,
  family = tweedie(link = "log")
  
); beepr::beep(2)

fit
max(fit$gradients)
AIC(fit)
tidy(fit, conf.int = T)
tidy(fit, effects = "ran_pars", conf.int = TRUE)
sanity(fit)

visreg::visreg(fit, xvar = "depth", xlim = c(0, 30))
visreg::visreg(fit, xvar = "depth", scale = "response", xlim = c(0, 30), nn = 200)

# Predict on new data
load("data/mhi_grid.rdata")
p <- predict(fit, newdata = grid)
head(p)

p %>% 
  group_by(X, Y, year) %>% 
  summarise(est = mean(est, na.rm = T)) %>% 
  ggplot(aes(X, Y, fill = exp(est), color = exp(est))) + 
  geom_point(shape = 22, alpha = 0.5, size = 0.5) + 
  scale_fill_gradientn(colors = matlab.like(100), trans = "sqrt") + 
  scale_color_gradientn(colors = matlab.like(100), trans = "sqrt") + 
  ggdark::dark_theme_minimal() + 
  coord_fixed() + 
  facet_wrap(~year)

p <- predict(fit, newdata = grid, return_tmb_object = T)
index <- get_index(p, area = rep(0.0081, nrow(grid)))

ggdark::invert_geom_defaults()

ggplot(index, aes(year, est, group = 1)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey90") +
  geom_line(lwd = 1, colour = "grey30") +
  labs(x = "Year", y = "Biomass (g)") + 
  theme_minimal()

cog <- get_cog(p, format = "wide")
ggplot(cog, aes(est_x, est_y, colour = year)) +
  geom_pointrange(aes(xmin = lwr_x, xmax = upr_x)) +
  geom_pointrange(aes(ymin = lwr_y, ymax = upr_y)) + 
  theme_minimal()

# Spatially varying effect of time:
df$year_scaled <- scale(as.numeric(df$year))
fit <- sdmTMB(
  response ~ s(depth, k=5) + year_scaled,
  spatial_varying = ~ year_scaled, 
  data = df, 
  silent = F, 
  mesh = rea_spde_coast, 
  time = "year",
  family = tweedie(link = "log"),
  spatiotemporal = "off"
)

grid$year_scaled <- (as.numeric(grid$year) - mean(as.numeric(df$year))) / sd(as.numeric(df$year))

p <- predict(fit, newdata = grid) %>% 
  subset(year == 2019) # any year

p %>% 
  ggplot(aes(X, Y, fill = zeta_s_year_scaled, color = zeta_s_year_scaled)) + 
  geom_point(shape = 22, alpha = 0.5, size = 0.5) + 
  scale_fill_gradient2() + 
  scale_color_gradient2() + 
  ggdark::dark_theme_minimal() + 
  coord_fixed() + 
  facet_wrap(~year)

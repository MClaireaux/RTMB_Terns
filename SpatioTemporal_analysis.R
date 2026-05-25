# ------------------------------------------------------------------------------
# 1. LOAD LIBRARIES AND DEFINE USEFUL FUNCTIONS
# ------------------------------------------------------------------------------
# Libraries
library(RTMB)
library(fmesher) # Replaces INLA. Used to create the mesh and precision matrix
library(Matrix)
library(INLA)
library(fields)
library(ggplot2)
library(sf)
library(dplyr)
library(tidyr)

# Function to return the first non-missing value of a vector 
# Used when aggregating observations from two different routes with the 
# same starting location.
first_non_na <- function(x) {
  x2 <- x[!is.na(x)]
  if (length(x2) == 0) NA else x2[1]
}

# Mat??rn covariance function (eq.2.8) 
matern_corr <- function(dist_matrix, kappa, nu) {
  
  a <- kappa*dist_matrix
  corr <- ((2^(1 - nu)) / gamma(nu))*(a^nu)*besselK(a, nu)
  diag(corr) <- 1
  
  return(corr)
}

# ------------------------------------------------------------------------------
# 2. LOAD DATA
# ------------------------------------------------------------------------------

# Reproducibility
set.seed(123)

# Set file path
path <- 'C:/Users~'
setwd(path)

# Import data
data <- read.csv('Tern_data.csv')

# ------------------------------------------------------------------------------
# 3. SELECT APPROPRIATE LOCATIONS FOR THE MESH CREATION
# ------------------------------------------------------------------------------

# Unique identifier for starting locations
loc_key  <- c("CountryNum", "StateNum", "Route", "Latitude", "Longitude")
# Create unique location-year pairs
id_cols  <- c(loc_key, "Year")
# Extract the names of the columns containing count data
count_cols <- grep("^Count\\d+$", names(data), value = TRUE)
# Sorted table of all the years present in the dataset
years_tbl <- tibble(Year = sort(unique(data$Year)))

### Keep one row per location-year pairs
# Goal: aggregate routes with similar starting points
# Counts are aggregated
# For other variables, mean across duplicates is computed

data_1row <- data %>%
  # Group by all the identifiers defined in id_cols
  group_by(across(all_of(id_cols))) %>%
  summarise(
    # sum counts across duplicates
    across(all_of(count_cols), ~ sum(.x, na.rm = TRUE)),
    # keep all other columns (weather, effort, etc.) 
    across(setdiff(names(data), c(id_cols, count_cols)), first_non_na),
    .groups = "drop"
  )

### Identify locations where Terns were observed at least once over 
### all the years
stations_with_terns <- data_1row %>%
  # Group by locations (not location-year)
  group_by(across(all_of(loc_key))) %>%
  # Returns FALSE is no terns were observed
  summarise(
    ever_terns = any(SpeciesTotal > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Filter and select locations where terns were observed at least once
  filter(ever_terns) %>%
  select(all_of(loc_key))

### Complete site-by-year grid: 
### Every station where a tern was sampled once x year
full_grid_terns <- stations_with_terns %>%
  crossing(years_tbl)

### Join observed data to the complete site-by-year grid.
### Years without surveys remain in the dataset with missing values (NA).
data_full <- full_grid_terns %>%
  left_join(data_1row, by = c(loc_key, "Year")) %>%
  as.data.frame()

### Final working dataset
Data <- subset(data_full,  Year > 2000) 

# ------------------------------------------------------------------------------
# 3. CREATE THE DISTANCE MATRIX - Section 3.1
# ------------------------------------------------------------------------------

### Distance conversion to kilometers
# Build a dataframe with the longitude and lattitude of all the relevant stations
# i.e. where Terns have been observed once over the study period
loc <- unique(as.data.frame(cbind(Data$Longitude,Data$Latitude))) 
colnames(loc) <- c("Longitude","Latitude")

# Convert the location table to an sf (simple features) spatial object 
pts_ll <- st_as_sf(loc, coords = c("Longitude","Latitude"), crs = 4326)

# Build a North America Albers Equal Area tuned to location envelope
# Find the westernmost, easternmost, southernmost, and northernmost survey locations.
lon_min <- min(Data$Longitude, na.rm = TRUE); lon_max <- max(Data$Longitude, na.rm = TRUE)
lat_min <- min(Data$Latitude, na.rm = TRUE);  lat_max <- max(Data$Latitude, na.rm = TRUE)

# Define center of the study area
lon0 <- (lon_min + lon_max) / 2      
lat0 <- (lat_min + lat_max) / 2      

# Define standard parallels, i.e. latitudes where projection distortion is minimized
# Common choice for North America:
lat1 <- 30
lat2 <- 60

# Build the Coordinate Reference system
na_aea <- st_crs(
  paste(
    # Albers Equal Area Projection
    "+proj=aea",
    # Standard parallels
    paste0("+lat_1=", lat1),           
    paste0("+lat_2=", lat2),
    # Center of the projection
    paste0("+lat_0=", lat0),
    paste0("+lon_0=", lon0),
    # standard GPS reference system
    "+datum=WGS84 +units=m +no_defs"
  )
)

# Convert coordinates into the AEA projection
pts_m <- st_transform(pts_ll, na_aea)
# Extract projected coordinates and convert them to kilometers
loc_m  <- st_coordinates(pts_m)/1000  # matrix [x(m), y(m)]

# Compute the distance matrix (eq.3.1)
dist_matrix <- fields::rdist(loc_m)

### Compute distance matrix characteristics
# Maximum distance between points
diameter   <- max(dist_matrix) 
# For each site, distance to the closest site
nn_dist    <- apply(dist_matrix + diag(Inf, nrow(dist_matrix)), 1, min)
# Median distance to the closest site
median_nn  <- median(nn_dist)


# ------------------------------------------------------------------------------
# 4. CREATE THE MESH AND SPDE - Section 3.2
# ------------------------------------------------------------------------------

### Define Mesh parameters
# Prevent multiple mesh nodes closer than 10 km from each other
cutoff <- 10     
# Maximum mesh triangle size in the inner hull
inner  <- 75 
# Maximum mesh triangle size in the outer hull
outer  <- inner * 4   
# Smoothing of inner boundary (how close it follows data points)
hull_inner <- inner * 3
# Smoothing of outer boundary
hull_outer <- hull_inner * 4

### Build inner and outer boundaries of the mesh domain
bnd1 <- fmesher::fm_nonconvex_hull(unique(loc_m), convex=hull_inner) 
bnd2 <- fmesher::fm_nonconvex_hull(unique(loc_m), convex = hull_outer)

### Build the mesh
mesh <- fmesher::fm_mesh_2d(
  loc.domain = unique(loc_m),
  boundary = list(bnd1, bnd2),
  min.angle= 24,
  max.edge = c(inner, outer),
  cutoff   = cutoff        
)

### Compute finite element matrices associated with the mesh using the 
### standard spde formulation, i.e. obtain C_0, G_1 and G_2 (eq. 2.8)
spde <- fmesher::fm_fem(mesh) 
# ------------------------------------------------------------------------------
# 5. BUILS DATASET WITH RESPONSE AND EXPLANATORY VARIABLES
# ------------------------------------------------------------------------------

### Remove rows with missing count values 
idx <- !is.na(Data$SpeciesTotal)
Data_obs <- Data[idx, ]

### Build the projection matrix A
# Project coordinates for all the observed site-year rows
pts_obs_ll <- st_as_sf( 
  Data_obs, 
  coords = c("Longitude", "Latitude"),
  crs = 4326,
  remove = FALSE
)
pts_obs_m <- st_transform(pts_obs_ll, na_aea)
loc_obs_m <- st_coordinates(pts_obs_m) / 1000

# Projection matrix for the likelihood
A <- inla.spde.make.A(
  mesh = mesh,
  loc = loc_obs_m
)

### Compute additional variable
# Effort
# Convert starting time and end time in minutes after midnight
data$start_min <- (data$StartTime %/% 100) * 60 + (data$StartTime %% 100)
data$end_min   <- (data$EndTime   %/% 100) * 60 + (data$EndTime   %% 100)
# Effort is time spent between the start and end time
data$effort <- data$end_min-data$start_min

# Convert Fahrenheit to Celsius
data[data$TempScale== "F",]$StartTemp <- (data[data$TempScale== "F",]$StartTemp-32)/1.8

# Only keep data meeting the BBS quality criterion
data <- subset(data, RunType > 0)

### Turn years into integer indices
Year <- sort(unique(Data_obs$Year))
year_id <- match(Data_obs$Year, Year)

### Replace temperature missing values by mean at location
Data_obs <- Data_obs %>%
  group_by(CountryNum, StateNum, Route) %>%
  mutate(
    StartTemp = ifelse(is.na(StartTemp), mean(StartTemp, na.rm = TRUE), StartTemp)
  ) %>%
  ungroup()

### Scale variables of interest
Data_obs$Temp_sc <- as.numeric(scale(Data_obs$StartTemp ))
Data_obs$Lat_cent <- as.numeric(scale(Data_obs$Latitude))
Data_obs$Year_cent <- as.numeric(scale(Data_obs$Year))
Data_obs$effort_cent <- as.numeric(scale(Data_obs$effort))

# ------------------------------------------------------------------------------
# 6. BUILD THE MODEL
# ------------------------------------------------------------------------------

### Define fixed effects
form <- ~ effort_cent+Temp_sc
X <- model.matrix(form, data = Data_obs)

### Define a list of parameters used by the model
data3 <- list(
  X       = as.matrix(X),                       # Fixed effects
  A       = A,                                  # Interpolation matrix
  Count   = as.numeric(Data_obs$SpeciesTotal),  # Total abundance at location
  year_id = as.integer(year_id),                # year index
  Temp_sc = Data_obs$Temp_sc,                   # Temperature
  n_years = length(Year),                       # Nb of years
  n_obs   = nrow(Data_obs),                     # Nb of rows
  c0      = spde$c0,                            # C_0 estimated from spde formulation
  g1      = spde$g1,                            # G_1 estimated from spde formulation
  g2      = spde$g2,                            # G_2 estimated from spde formulation
  n_spde  = nrow(spde$c0)                       # Number of rows of the spde matrix
)

### Define a list of parameters to be estimated with starting values
parameters3 <- list(
  beta      = rep(0, ncol(X)),         # weigths of the fixed effects
  log_tau   = 1.45,                    # tau from eq.2.8 (on log scale)
  log_kappa = -2.77,                   # kappa from eq.2.8 (on log scale)
  logit_phi = qlogis(0.9),             # temporal correlation
  log_alpha = 3,                       # Dispersion parameter
  xt        = rep(0.0, data3$n_spde * data3$n_years) # spatial field
)


### Write a function describing the model
f <- function(parms) {
    getAll(parms, data3)
    
   # Transform back parameters' starting values
    tau   <- exp(log_tau)
    kappa <- exp(log_kappa)
    alpha <- exp(log_alpha) 
    phi <- plogis(logit_phi) 
    
    # Compute  the precision matrix Q (eq. 2.8)
    Q <- tau^2 * (kappa^4 * c0 + 2 * kappa^2 * g1 + g2)
    
    # Convert the latent spatial effects vector into a matrix (space x time)
    Xmat <- matrix(xt, n_spde, n_years)
    
    # Initiate the negative log-likelihood
    nll <- 0
    # Prior for the initial spatial field under the stationary AR(1) distribution
    # (eq. 3.2)
    nll <- nll - dgmrf(Xmat[,1], 0, Q * (1 - phi^2), log = TRUE)
    # AR(1) temporal evolution of the latent spatial field (eq. 3.3)
    for (t in 2:n_years) {
      nll <- nll - dgmrf(Xmat[,t], phi * Xmat[,t-1], Q, log = TRUE)
    }
    
    # Within each year,:
    for (t in 1:n_years) {
      rows_t <- which(year_id == t)  
      # Define the linear predictor (eq. 3.4)
      eta_t <- as.numeric(X[rows_t, , drop = FALSE] %*% beta) +
        as.numeric(A[rows_t, , drop = FALSE] %*% Xmat[, t])
      #Calculate the mean and variance of the NB distribution (Table 5.1)
      lambda <- exp(eta_t)
      var_nb <- lambda + (lambda^2)*alpha
      # Response
      y <- Count[rows_t]
      # Define the observation model (eq. 3.5)
      lnb_y <- dnbinom2(y, mu = lambda, var = var_nb, log = TRUE)
      
      # Update the negative likelihood
      nll <- nll - sum(lnb_y)
      
    }
    
    # Define relevant quantities
    nu <- 1                            # Smoothness nu (eq. 2.3)   
    rho_est    <- sqrt(8 * nu) / kappa # range (eq. 2.6)
    sigma_est  <- 1 / (sqrt(4 * pi) * kappa * tau) # Marginal variance (eq. 2.5)
    sigma2_est <- sigma_est^2
    
    # Report the estimates
    ADREPORT(rho_est)
    ADREPORT(sigma2_est)
    ADREPORT(phi) 

  return(nll)
  }
  
# ------------------------------------------------------------------------------
# 7. RUN THE MODEL AND SAVE THE OUTPUT
# ------------------------------------------------------------------------------

### Build the RTMB objective function.
# 'xt' is treated as a random effect and will be integrated out using the 
# Laplace approximation
obj3 <- MakeADFun(f, parameters3, random="xt")

### Optimize fixed effects and hyperparameters by minimizing the negative 
### log-likelihood
opt3 <- nlminb(obj3$par, obj3$fn, obj3$gr)

### Check optimizer convergence status
# (0 for successful convergence)
print(opt3$convergence)

###Print estimated parameter values
print(opt3$par)

### Compute standard errors and uncertainty estimates
Rep3 <- sdreport(obj3)

### Parameters summaries
# Fixed effects and hyperparameters
fixed3 <- summary(Rep3, "fixed")
# Random effects (spatial field)
rand <- summary(Rep3, "random")
# Summary of derived quantities reported via ADREPORT()
estimates3 <- summary(Rep3,"report")

### Store all important model outputs in a single object
results <- list(
  model      = form,
  opt        = opt3,
  fixed      = summary(Rep3, "fixed"),
  report     = summary(Rep3, "report"),
  random     = summary(Rep3, "random"),
  par        = opt3$par,
  convergence= opt3$convergence,
  Year = Year,
  A = A,
  spde = spde
)

# Save model outputs
save_path <- 'C:/Users/~'
setwd(save_path)

filename <- 'filename.RDS'
saveRDS(results, file = "filename")

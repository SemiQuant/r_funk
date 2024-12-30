# Create example qPCR dataset
set.seed(42)

# Create sample data
samples <- rep(c("Control", "Treatment_A", "Treatment_B", "Treatment_C"), each = 6)
targets <- rep(c("Target_Gene", "Reference_Gene"), each = 3, times = 4)
replicates <- rep(1:3, times = 8)

# Generate Ct values
# Control: baseline expression
# Treatment_A: 2-fold increase
# Treatment_B: 4-fold increase
# Treatment_C: 0.5-fold decrease

# Function to add noise to Ct values
add_noise <- function(ct, sd = 0.3) {
  ct + rnorm(length(ct), mean = 0, sd = sd)
}

# Base Ct values
base_target_ct <- 25
base_ref_ct <- 20

# Generate Ct values with biological effects and technical variation
cts <- numeric(24)

# Reference gene - should be stable across conditions
ref_indices <- targets == "Reference_Gene"
cts[ref_indices] <- add_noise(rep(base_ref_ct, sum(ref_indices)))

# Target gene - varies by treatment
target_indices <- targets == "Target_Gene"
target_cts <- numeric(sum(target_indices))

# Control
target_cts[1:3] <- add_noise(rep(base_target_ct, 3))
# Treatment A (2-fold increase = -1 Ct)
target_cts[4:6] <- add_noise(rep(base_target_ct - 1, 3))
# Treatment B (4-fold increase = -2 Ct)
target_cts[7:9] <- add_noise(rep(base_target_ct - 2, 3))
# Treatment C (0.5-fold decrease = +1 Ct)
target_cts[10:12] <- add_noise(rep(base_target_ct + 1, 3))

cts[target_indices] <- target_cts

# Create the data frame
qpcr_example <- data.frame(
  Sample = samples,
  Target = targets,
  Cq = cts,
  Omit = FALSE,
  Replicate = replicates
)

# Add one outlier to demonstrate omission
qpcr_example$Cq[1] <- 35
qpcr_example$Omit[1] <- TRUE

# Save the dataset
usethis::use_data(qpcr_example, overwrite = TRUE) 
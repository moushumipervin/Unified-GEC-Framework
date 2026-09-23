###############################################################################
# TABLE 10: NSW / PSID / CPS REAL-DATA APPLICATION
###############################################################################

# Paper benchmark used for evaluation bias
paper_benchmark <- 1794


# -----------------------------------------------------------------------------
# Helper to extract one estimate and SE
# -----------------------------------------------------------------------------

get_result <- function(
    donor_name,
    method_name
) {

  z <- final_long |>
    dplyr::filter(
      as.character(Donor) == donor_name,
      as.character(Method) == method_name
    )

  if (nrow(z) != 1L) {

    return(
      c(
        Estimate = NA_real_,
        SE = NA_real_
      )
    )
  }

  c(
    Estimate = z$Estimate[1],
    SE = z$SE[1]
  )
}


# -----------------------------------------------------------------------------
# Method order exactly as in manuscript Table 10
# -----------------------------------------------------------------------------

table10_methods <- c(
  "Unweighted",
  "IPW",
  "EBPS",
  "CBPS",
  "EBCW",
  "AIPW (LM)",
  "AIPW (GAM)",
  "ET",
  "HD",
  "CE"
)


# -----------------------------------------------------------------------------
# Build manuscript Table 10
# -----------------------------------------------------------------------------

table10_rows <- lapply(
  table10_methods,
  function(m) {

    nsw <- get_result(
      "NSW",
      m
    )

    psid <- get_result(
      "PSID",
      m
    )

    cps <- get_result(
      "CPS",
      m
    )

    data.frame(

      Estimator = m,

      NSW_Est = round(
        nsw["Estimate"],
        0
      ),

      NSW_SE = round(
        nsw["SE"],
        0
      ),

      PSID_Est = round(
        psid["Estimate"],
        0
      ),

      PSID_EB = round(
        psid["Estimate"] -
          paper_benchmark,
        0
      ),

      PSID_SE = round(
        psid["SE"],
        0
      ),

      CPS_Est = round(
        cps["Estimate"],
        0
      ),

      CPS_EB = round(
        cps["Estimate"] -
          paper_benchmark,
        0
      ),

      CPS_SE = round(
        cps["SE"],
        0
      ),

      stringsAsFactors = FALSE
    )
  }
)


table10 <- dplyr::bind_rows(
  table10_rows
)


# -----------------------------------------------------------------------------
# Rename AIPW labels to match manuscript
# -----------------------------------------------------------------------------

table10 <- table10 |>
  dplyr::mutate(
    Estimator = dplyr::recode(
      Estimator,
      "AIPW (LM)" = "AIPW-LM",
      "AIPW (GAM)" = "AIPW-GAM"
    )
  )


# -----------------------------------------------------------------------------
# Save Table 10
# -----------------------------------------------------------------------------

dir.create(
  "Causal Inference/results/table10",
  recursive = TRUE,
  showWarnings = FALSE
)


write.csv(
  table10,
  "Causal Inference/results/table10/table10_lalonde.csv",
  row.names = FALSE
)


# -----------------------------------------------------------------------------
# Print Table 10
# -----------------------------------------------------------------------------

print(
  table10,
  row.names = FALSE
)

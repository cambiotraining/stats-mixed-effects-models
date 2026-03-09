library(shiny)
library(lmerTest)
library(MASS)   # for mvrnorm
library(ggplot2)

ui <- fluidPage(
  titlePanel("Power Analysis for Linear Mixed Models with Correlated Random Effects"),
  
  sidebarLayout(
    sidebarPanel(
      h4("Sample Design"),
      numericInput("n_clusters", "Number of clusters:", value = 20, min = 2),
      numericInput("n_obs", "Observations per cluster:", value = 10, min = 1),
      
      h4("Fixed Effects"),
      selectInput("fx1_type", "Fixed Effect 1 Type:", choices = c("None", "Continuous", "Categorical")),
      conditionalPanel("input.fx1_type == 'Categorical'",
                       numericInput("fx1_levels", "Levels for Fixed Effect 1:", value = 2, min = 2),
                       conditionalPanel("input.fx1_levels > 2",
                                        textInput("fx1_beta_vec", 
                                                  "Fixed Effect 1 Betas for contrasts (comma-separated):", 
                                                  value = "")
                       )
      ),
      conditionalPanel("input.fx1_type != 'None'",
                       numericInput("fx1_beta", "Effect Size (β) for Fixed Effect 1:", value = 0.5),
                       checkboxInput("fx1_random_slope", "Include random slope for Fixed Effect 1?", value = FALSE)
      ),
      
      selectInput("fx2_type", "Fixed Effect 2 Type:", choices = c("None", "Continuous", "Categorical")),
      conditionalPanel("input.fx2_type == 'Categorical'",
                       numericInput("fx2_levels", "Levels for Fixed Effect 2:", value = 2, min = 2),
                       conditionalPanel("input.fx2_levels > 2",
                                        textInput("fx2_beta_vec", 
                                                  "Fixed Effect 2 Betas for contrasts (comma-separated):", 
                                                  value = "")
                       )
      ),
      conditionalPanel("input.fx2_type != 'None'",
                       numericInput("fx2_beta", "Effect Size (β) for Fixed Effect 2:", value = 0.5),
                       checkboxInput("fx2_random_slope", "Include random slope for Fixed Effect 2?", value = FALSE)
      ),
      
      h4("Variance Components (Random Effects)"),
      numericInput("sd_intercept", "Random Intercept SD:", value = 1, min = 0),
      
      conditionalPanel("input.fx1_random_slope == true",
                       numericInput("sd_slope_fx1", "Random Slope SD for Fixed Effect 1:", value = 0.2, min = 0)
      ),
      conditionalPanel("input.fx2_random_slope == true",
                       numericInput("sd_slope_fx2", "Random Slope SD for Fixed Effect 2:", value = 0.2, min = 0)
      ),
      
      conditionalPanel("input.fx1_random_slope == true && input.fx2_random_slope == true",
                       numericInput("cor_fx1_fx2", "Correlation between Random Slopes (fx1 & fx2):", value = 0, min = -1, max = 1)
      ),
      conditionalPanel("input.fx1_random_slope == true",
                       numericInput("cor_intercept_fx1", "Correlation between Intercept and fx1 Slope:", value = 0, min = -1, max = 1)
      ),
      conditionalPanel("input.fx2_random_slope == true",
                       numericInput("cor_intercept_fx2", "Correlation between Intercept and fx2 Slope:", value = 0, min = -1, max = 1)
      ),
      
      h4("Residual"),
      numericInput("sd_resid", "Residual SD (within-cluster):", value = 1, min = 0),
      
      h4("Response Variable Distribution"),
      numericInput("response_mean", "Expected response mean:", value = 0),
      numericInput("response_sd", "Expected response SD:", value = 1, min = 0.0001),
      
      h4("Significance Testing Method"),
      selectInput("test_method", "Method:",
                  choices = c("Satterthwaite (lmerTest)" = "satterthwaite",
                              "Likelihood Ratio Test (LRT)" = "lrt")),
      
      h4("Simulation Settings"),
      numericInput("alpha", "Significance level (α):", value = 0.05),
      numericInput("n_sim", "Number of simulations:", value = 100, min = 10),
      
      actionButton("run_sim", "Run Simulation")
    ),
    
    mainPanel(
      h4("Results"),
      verbatimTextOutput("power_output"),
      plotOutput("power_plot")
    )
  )
)

server <- function(input, output, session) {
  
  parse_beta_vec <- function(text, expected_len) {
    if (nchar(text) == 0) return(NULL)
    vals <- suppressWarnings(as.numeric(strsplit(text, ",")[[1]]))
    if (length(vals) != expected_len || any(is.na(vals))) return(NULL)
    vals
  }
  
  simulate_once <- function(params) {
    with(params, {
      n_total <- n_clusters * n_obs
      cluster <- factor(rep(1:n_clusters, each = n_obs))
      
      # Prepare fixed effect predictors and model matrices
      fx1 <- fx2 <- rep(0, n_total)
      mm1 <- mm2 <- NULL
      
      # FX1
      if (fx1_type == "Continuous") {
        x1 <- rnorm(n_total)
        fx1 <- x1
      } else if (fx1_type == "Categorical") {
        x1 <- factor(sample(1:fx1_levels, n_total, replace = TRUE))
        fx1 <- x1
        contrasts(fx1) <- contr.sum
        mm1 <- model.matrix(~ fx1)
        n_contrasts1 <- ncol(mm1) - 1
      } else {
        x1 <- NULL
      }
      
      # FX2
      if (fx2_type == "Continuous") {
        x2 <- rnorm(n_total)
        fx2 <- x2
      } else if (fx2_type == "Categorical") {
        x2 <- factor(sample(1:fx2_levels, n_total, replace = TRUE))
        fx2 <- x2
        contrasts(x2) <- contr.sum
        mm2 <- model.matrix(~ fx2)
        n_contrasts2 <- ncol(mm2) - 1
      } else {
        x2 <- NULL
      }
      
      # --- Construct covariance matrix for random effects ---
      
      # Determine random effects included
      re_terms <- list()
      re_names <- c()
      re_sd <- c()
      # Always intercept first
      re_terms <- c(re_terms, "intercept")
      re_names <- c(re_names, "Intercept")
      re_sd <- c(re_sd, sd_intercept)
      
      # fx1 slope if continuous and random slope checked
      if (fx1_random_slope && fx1_type == "Continuous") {
        re_terms <- c(re_terms, "fx1_slope")
        re_names <- c(re_names, "fx1_slope")
        re_sd <- c(re_sd, sd_slope_fx1)
      }
      
      # fx2 slope if continuous and random slope checked
      if (fx2_random_slope && fx2_type == "Continuous") {
        re_terms <- c(re_terms, "fx2_slope")
        re_names <- c(re_names, "fx2_slope")
        re_sd <- c(re_sd, sd_slope_fx2)
      }
      
      n_re <- length(re_terms)
      
      # Build correlation matrix (initialize identity)
      cor_mat <- diag(n_re)
      rownames(cor_mat) <- colnames(cor_mat) <- re_terms
      
      # Fill correlations from input if applicable
      # intercept - fx1
      if ("intercept" %in% re_terms && "fx1_slope" %in% re_terms) {
        cor_mat["intercept", "fx1_slope"] <- cor_intercept_fx1
        cor_mat["fx1_slope", "intercept"] <- cor_intercept_fx1
      }
      # intercept - fx2
      if ("intercept" %in% re_terms && "fx2_slope" %in% re_terms) {
        cor_mat["intercept", "fx2_slope"] <- cor_intercept_fx2
        cor_mat["fx2_slope", "intercept"] <- cor_intercept_fx2
      }
      # fx1 - fx2
      if ("fx1_slope" %in% re_terms && "fx2_slope" %in% re_terms) {
        cor_mat["fx1_slope", "fx2_slope"] <- cor_fx1_fx2
        cor_mat["fx2_slope", "fx1_slope"] <- cor_fx1_fx2
      }
      
      # Build covariance matrix from sds and correlations
      sd_mat <- diag(re_sd)
      cov_mat <- sd_mat %*% cor_mat %*% sd_mat
      
      # Check covariance matrix positive definite
      if (any(eigen(cov_mat)$values <= 0)) {
        stop("Covariance matrix of random effects is not positive definite. Adjust SDs or correlations.")
      }
      
      # Simulate random effects per cluster
      random_effects <- MASS::mvrnorm(n = n_clusters, mu = rep(0, n_re), Sigma = cov_mat)
      colnames(random_effects) <- re_terms
      
      # Build linear predictor
      lin_pred <- rep(0, n_total)
      
      # Add fixed effects
      if (fx1_type == "Continuous") lin_pred <- lin_pred + fx1_beta * fx1
      else if (fx1_type == "Categorical") {
        betas1 <- parse_beta_vec(fx1_beta_vec, n_contrasts1)
        if (is.null(betas1)) betas1 <- rep(fx1_beta, n_contrasts1) # equal effects default
        mm1_beta <- mm1[, -1, drop = FALSE] %*% betas1
        lin_pred <- lin_pred + mm1_beta
      }
      
      if (fx2_type == "Continuous") lin_pred <- lin_pred + fx2_beta * fx2
      else if (fx2_type == "Categorical") {
        betas2 <- parse_beta_vec(fx2_beta_vec, n_contrasts2)
        if (is.null(betas2)) betas2 <- rep(fx2_beta, n_contrasts2)
        mm2_beta <- mm2[, -1, drop = FALSE] %*% betas2
        lin_pred <- lin_pred + mm2_beta
      }
      
      # Add random effects
      for (i in 1:n_total) {
        cl <- as.integer(cluster[i])
        re_vec <- random_effects[cl, ]
        
        re_val <- re_vec["intercept"]
        if ("fx1_slope" %in% re_terms && fx1_type == "Continuous") re_val <- re_val + re_vec["fx1_slope"] * fx1[i]
        if ("fx2_slope" %in% re_terms && fx2_type == "Continuous") re_val <- re_val + re_vec["fx2_slope"] * fx2[i]
        
        lin_pred[i] <- lin_pred[i] + re_val
      }
      
      # Residual errors
      y <- lin_pred + rnorm(n_total, 0, sd_resid)
      
      # Rescale y to have specified mean and SD
      y <- scale(y)
      y <- y * response_sd + response_mean
      
      df <- data.frame(y = y, cluster = cluster)
      if (fx1_type == "Continuous") df$fx1 <- fx1
      if (fx1_type == "Categorical") df$fx1 <- fx1
      if (fx2_type == "Continuous") df$fx2 <- fx2
      if (fx2_type == "Categorical") df$fx2 <- fx2
      
      # Build formula string for lmer
      fixed_terms <- c()
      if (fx1_type != "None") fixed_terms <- c(fixed_terms, "fx1")
      if (fx2_type != "None") fixed_terms <- c(fixed_terms, "fx2")
      
      # Only random slopes for continuous fx1/fx2
      rand_slopes <- c()
      if (fx1_random_slope && fx1_type == "Continuous") rand_slopes <- c(rand_slopes, "fx1")
      if (fx2_random_slope && fx2_type == "Continuous") rand_slopes <- c(rand_slopes, "fx2")
      
      if (length(rand_slopes) > 0) {
        rand_str <- paste0("(1 + ", paste(rand_slopes, collapse = " + "), " | cluster)")
      } else {
        rand_str <- "(1 | cluster)"
      }
      
      formula_str <- paste("y ~", paste(fixed_terms, collapse = " + "), "+", rand_str)
      formula_str <- gsub("\\+\\s+\\+", "+", formula_str) # clean up double plus
      formula_str <- gsub("\\+\\s+$", "", formula_str) # clean trailing plus
      formula_str <- gsub("^\\s*\\+\\s*", "", formula_str) # clean leading plus
      formula <- as.formula(formula_str)
      
      # Fit model
      model <- try(lmerTest::lmer(formula, data = df, REML = FALSE), silent = TRUE)
      if (inherits(model, "try-error")) stop("Model fitting failed")
      
      # Extract p-values
      p1 <- NA
      p2 <- NA
      
      if (test_method == "satterthwaite") {
        summ <- summary(model)
        coefs <- summ$coefficients
        if (fx1_type != "None") {
          idx <- grep("^fx1", rownames(coefs))
          if (length(idx) > 0) p1 <- coefs[idx[1], "Pr(>|t|)"]
        }
        if (fx2_type != "None") {
          idx <- grep("^fx2", rownames(coefs))
          if (length(idx) > 0) p2 <- coefs[idx[1], "Pr(>|t|)"]
        }
      } else if (test_method == "lrt") {
        # Fit reduced models dropping each fixed effect
        if (fx1_type != "None") {
          formula_r1 <- update(formula, paste(". ~ . - fx1"))
          m2 <- try(lmerTest::lmer(formula_r1, data = df, REML = FALSE), silent = TRUE)
          if (!inherits(m2, "try-error")) {
            p1 <- anova(m2, model)[["Pr(>Chisq)"]][2]
          }
        }
        if (fx2_type != "None") {
          formula_r2 <- update(formula, paste(". ~ . - fx2"))
          m2 <- try(lmerTest::lmer(formula_r2, data = df, REML = FALSE), silent = TRUE)
          if (!inherits(m2, "try-error")) {
            p2 <- anova(m2, model)[["Pr(>Chisq)"]][2]
          }
        }
      }
      
      c(p1, p2)
    })
  }
  
  observeEvent(input$run_sim, {
    req(input$n_sim >= 10)
    
    params <- list(
      n_clusters = input$n_clusters,
      n_obs = input$n_obs,
      fx1_type = input$fx1_type,
      fx1_levels = ifelse(is.null(input$fx1_levels), 2, input$fx1_levels),
      fx1_beta = input$fx1_beta,
      fx1_beta_vec = input$fx1_beta_vec,
      fx1_random_slope = input$fx1_random_slope,
      sd_slope_fx1 = input$sd_slope_fx1,
      fx2_type = input$fx2_type,
      fx2_levels = ifelse(is.null(input$fx2_levels), 2, input$fx2_levels),
      fx2_beta = input$fx2_beta,
      fx2_beta_vec = input$fx2_beta_vec,
      fx2_random_slope = input$fx2_random_slope,
      sd_slope_fx2 = input$sd_slope_fx2,
      sd_intercept = input$sd_intercept,
      cor_intercept_fx1 = input$cor_intercept_fx1,
      cor_intercept_fx2 = input$cor_intercept_fx2,
      cor_fx1_fx2 = input$cor_fx1_fx2,
      sd_resid = input$sd_resid,
      response_mean = input$response_mean,
      response_sd = input$response_sd,
      test_method = input$test_method
    )
    
    withProgress(message = "Running simulations...", value = 0, {
      pvalues <- replicate(input$n_sim, {
        ps <- try(simulate_once(params), silent = TRUE)
        if (inherits(ps, "try-error") || any(is.na(ps))) return(c(NA, NA))
        ps
      })
      
      # Filter out failed runs
      pvalues <- pvalues[, colSums(is.na(pvalues)) == 0, drop = FALSE]
      
      power_fx1 <- NA
      power_fx2 <- NA
      
      if (input$fx1_type != "None" && ncol(pvalues) > 0) {
        power_fx1 <- mean(pvalues[1, ] < input$alpha, na.rm = TRUE)
      }
      if (input$fx2_type != "None" && ncol(pvalues) > 0) {
        power_fx2 <- mean(pvalues[2, ] < input$alpha, na.rm = TRUE)
      }
      
      output$power_output <- renderPrint({
        cat("Power estimates based on", ncol(pvalues), "simulations:\n\n")
        if (input$fx1_type != "None") cat("Fixed Effect 1 Power:", round(power_fx1, 3), "\n")
        if (input$fx2_type != "None") cat("Fixed Effect 2 Power:", round(power_fx2, 3), "\n")
      })
      
      output$power_plot <- renderPlot({
        barplot <- data.frame(
          Effect = c("Fixed Effect 1", "Fixed Effect 2"),
          Power = c(power_fx1, power_fx2)
        )
        barplot <- barplot[!is.na(barplot$Power), ]
        ggplot(barplot, aes(x = Effect, y = Power, fill = Effect)) +
          geom_col() +
          ylim(0, 1) +
          ylab("Power") +
          theme_minimal() +
          theme(legend.position = "none")
      })
    })
  })
}

shinyApp(ui, server)

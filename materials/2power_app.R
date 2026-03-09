library(shiny)
library(lme4)
library(lmerTest)
library(MASS)
library(ggplot2)

ui <- fluidPage(
  titlePanel("Power Analysis for Mixed Models by Simulation"),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("n_clusters", "Number of Clusters", 20),
      numericInput("n_obs", "Observations per Cluster (if equal)", 10),
      textInput("cluster_sizes", "Cluster Sample Sizes (comma-separated)", ""),
      
      selectInput("fx1_type", "Fixed Effect 1 Type", c("None", "Continuous", "Categorical")),
      conditionalPanel(
        condition = "input.fx1_type == 'Categorical'",
        numericInput("fx1_levels", "Levels for Fixed Effect 1", 3),
        textInput("fx1_beta_vec", "Beta Coefficients (comma-separated)", "")
      ),
      conditionalPanel(
        condition = "input.fx1_type != 'None'",
        numericInput("fx1_beta", "Fixed Effect 1 Coefficient (default if vector empty)", 0.5),
        checkboxInput("fx1_random_slope", "Random Slope for Fixed Effect 1", FALSE),
        numericInput("sd_slope_fx1", "SD of Random Slope for Fixed Effect 1", 0.1)
      ),
      
      selectInput("fx2_type", "Fixed Effect 2 Type", c("None", "Continuous", "Categorical")),
      conditionalPanel(
        condition = "input.fx2_type == 'Categorical'",
        numericInput("fx2_levels", "Levels for Fixed Effect 2", 3),
        textInput("fx2_beta_vec", "Beta Coefficients (comma-separated)", "")
      ),
      conditionalPanel(
        condition = "input.fx2_type != 'None'",
        numericInput("fx2_beta", "Fixed Effect 2 Coefficient (default if vector empty)", 0.5),
        checkboxInput("fx2_random_slope", "Random Slope for Fixed Effect 2", FALSE),
        numericInput("sd_slope_fx2", "SD of Random Slope for Fixed Effect 2", 0.1)
      ),
      
      numericInput("sd_intercept", "SD of Random Intercepts", 1),
      numericInput("cor_intercept_fx1", "Correlation (Intercept, FX1)", 0),
      numericInput("cor_intercept_fx2", "Correlation (Intercept, FX2)", 0),
      numericInput("cor_fx1_fx2", "Correlation (FX1, FX2)", 0),
      numericInput("sd_resid", "Residual SD", 1),
      numericInput("response_mean", "Expected Response Mean", 0),
      numericInput("response_sd", "Expected Response SD", 1),
      selectInput("test_method", "Test Method", c("satterthwaite", "lrt")),
      numericInput("n_sim", "Number of Simulations", 100),
      numericInput("alpha", "Significance Level", 0.05),
      actionButton("run_sim", "Run Simulation")
    ),
    
    mainPanel(
      verbatimTextOutput("power_output"),
      plotOutput("power_plot")
    )
  )
)

server <- function(input, output) {
  
  parse_beta_vec <- function(beta_vec_str, n) {
    if (beta_vec_str == "") return(NULL)
    vec <- as.numeric(strsplit(beta_vec_str, ",")[[1]])
    if (length(vec) != n) return(NULL)
    vec
  }
  
  simulate_once <- function(params) {
    with(params, {
      n_contrasts1 <- if (fx1_type == "Categorical") fx1_levels - 1 else 0
      n_contrasts2 <- if (fx2_type == "Categorical") fx2_levels - 1 else 0
      
      cluster_sizes <- if (!is.null(cluster_sizes) && cluster_sizes != "") {
        as.numeric(strsplit(cluster_sizes, ",")[[1]])
      } else {
        rep(n_obs, n_clusters)
      }
      
      if (length(cluster_sizes) != n_clusters) stop("Length of cluster_sizes must equal number of clusters")
      cluster <- unlist(mapply(rep, x = 1:n_clusters, times = cluster_sizes))
      n_total <- length(cluster)
      
      fx1 <- if (fx1_type == "Continuous") rnorm(n_total) else if (fx1_type == "Categorical") factor(sample(1:fx1_levels, n_total, TRUE)) else NULL
      fx2 <- if (fx2_type == "Continuous") rnorm(n_total) else if (fx2_type == "Categorical") factor(sample(1:fx2_levels, n_total, TRUE)) else NULL
      
      lin_pred <- rep(0, n_total)
      re_terms <- c("intercept")
      
      if (fx1_random_slope) {
        if (fx1_type == "Continuous") re_terms <- c(re_terms, "fx1_slope")
        else if (fx1_type == "Categorical") re_terms <- c(re_terms, paste0("fx1_slope_", 1:n_contrasts1))
      }
      if (fx2_random_slope) {
        if (fx2_type == "Continuous") re_terms <- c(re_terms, "fx2_slope")
        else if (fx2_type == "Categorical") re_terms <- c(re_terms, paste0("fx2_slope_", 1:n_contrasts2))
      }
      
      sd_vec <- c(sd_intercept,
                  rep(sd_slope_fx1, if (fx1_random_slope && fx1_type == "Continuous") 1 else if (fx1_random_slope && fx1_type == "Categorical") n_contrasts1 else 0),
                  rep(sd_slope_fx2, if (fx2_random_slope && fx2_type == "Continuous") 1 else if (fx2_random_slope && fx2_type == "Categorical") n_contrasts2 else 0))
      
      cor_mat <- diag(length(sd_vec))
      if (length(sd_vec) > 1) {
        cor_mat[1, 2] <- cor_mat[2, 1] <- cor_intercept_fx1
        if (length(sd_vec) > 2) {
          cor_mat[1, 3] <- cor_mat[3, 1] <- cor_intercept_fx2
          cor_mat[2, 3] <- cor_mat[3, 2] <- cor_fx1_fx2
        }
      }
      
      cov_mat <- diag(sd_vec) %*% cor_mat %*% diag(sd_vec)
      random_effects <- MASS::mvrnorm(n_clusters, mu = rep(0, length(sd_vec)), Sigma = cov_mat)
      colnames(random_effects) <- re_terms
      
      mm1 <- if (fx1_type == "Categorical") model.matrix(~fx1) else NULL
      mm2 <- if (fx2_type == "Categorical") model.matrix(~fx2) else NULL
      
      if (fx1_type == "Continuous") lin_pred <- lin_pred + fx1_beta * fx1
      if (fx1_type == "Categorical") {
        betas1 <- parse_beta_vec(fx1_beta_vec, n_contrasts1)
        if (is.null(betas1)) betas1 <- rep(fx1_beta, n_contrasts1)
        mm1_beta <- mm1[, -1, drop = FALSE] %*% betas1
        lin_pred <- lin_pred + mm1_beta
      }
      
      if (fx2_type == "Continuous") lin_pred <- lin_pred + fx2_beta * fx2
      if (fx2_type == "Categorical") {
        betas2 <- parse_beta_vec(fx2_beta_vec, n_contrasts2)
        if (is.null(betas2)) betas2 <- rep(fx2_beta, n_contrasts2)
        mm2_beta <- mm2[, -1, drop = FALSE] %*% betas2
        lin_pred <- lin_pred + mm2_beta
      }
      
      for (i in 1:n_total) {
        cl <- cluster[i]
        re_vec <- random_effects[cl, ]
        re_val <- re_vec["intercept"]
        
        if (fx1_random_slope && fx1_type == "Continuous") re_val <- re_val + re_vec["fx1_slope"] * fx1[i]
        if (fx2_random_slope && fx2_type == "Continuous") re_val <- re_val + re_vec["fx2_slope"] * fx2[i]
        
        if (fx1_random_slope && fx1_type == "Categorical") {
          mm_row <- mm1[i, -1]
          for (j in 1:n_contrasts1) {
            re_val <- re_val + mm_row[j] * re_vec[paste0("fx1_slope_", j)]
          }
        }
        if (fx2_random_slope && fx2_type == "Categorical") {
          mm_row <- mm2[i, -1]
          for (j in 1:n_contrasts2) {
            re_val <- re_val + mm_row[j] * re_vec[paste0("fx2_slope_", j)]
          }
        }
        
        lin_pred[i] <- lin_pred[i] + re_val
      }
      
      y <- lin_pred + rnorm(n_total, 0, sd_resid)
      y <- scale(y)
      y <- y * response_sd + response_mean
      
      df <- data.frame(y = y, cluster = factor(cluster), fx1 = fx1, fx2 = fx2)
      
      fixed_terms <- c()
      if (fx1_type != "None") fixed_terms <- c(fixed_terms, "fx1")
      if (fx2_type != "None") fixed_terms <- c(fixed_terms, "fx2")
      
      rand_slopes <- c()
      if (fx1_random_slope) rand_slopes <- c(rand_slopes, "fx1")
      if (fx2_random_slope) rand_slopes <- c(rand_slopes, "fx2")
      
      rand_str <- if (length(rand_slopes) > 0) paste0("(1 + ", paste(rand_slopes, collapse = " + "), " | cluster)") else "(1 | cluster)"
      formula <- as.formula(paste("y ~", paste(fixed_terms, collapse = " + "), "+", rand_str))
      
      model <- try(lmerTest::lmer(formula, data = df, REML = FALSE), silent = TRUE)
      if (inherits(model, "try-error")) stop("Model fitting failed")
      
      p1 <- p2 <- NA
      
      if (test_method == "satterthwaite") {
        coefs <- summary(model)$coefficients
        if (fx1_type != "None") p1 <- coefs[grep("^fx1", rownames(coefs))[1], "Pr(>|t|)"]
        if (fx2_type != "None") p2 <- coefs[grep("^fx2", rownames(coefs))[1], "Pr(>|t|)"]
      } else {
        if (fx1_type != "None") {
          m2 <- try(update(model, . ~ . - fx1), silent = TRUE)
          if (!inherits(m2, "try-error")) p1 <- anova(m2, model)[["Pr(>Chisq)"]][2]
        }
        if (fx2_type != "None") {
          m2 <- try(update(model, . ~ . - fx2), silent = TRUE)
          if (!inherits(m2, "try-error")) p2 <- anova(m2, model)[["Pr(>Chisq)"]][2]
        }
      }
      
      c(p1, p2)
    })
  }
  
  observeEvent(input$run_sim, {
    req(input$n_sim >= 10)
    
    params <- reactiveValuesToList(input)
    withProgress(message = "Running simulations...", value = 0, {
      pvalues <- replicate(input$n_sim, {
        ps <- try(simulate_once(params), silent = TRUE)
        if (inherits(ps, "try-error") || any(is.na(ps))) return(c(NA, NA))
        ps
      })
      
      pvalues <- pvalues[, colSums(is.na(pvalues)) == 0, drop = FALSE]
      
      power_fx1 <- if (input$fx1_type != "None" && ncol(pvalues) > 0) mean(pvalues[1, ] < input$alpha, na.rm = TRUE) else NA
      power_fx2 <- if (input$fx2_type != "None" && ncol(pvalues) > 0) mean(pvalues[2, ] < input$alpha, na.rm = TRUE) else NA
      
      output$power_output <- renderPrint({
        cat("Power estimates based on", ncol(pvalues), "simulations:\n\n")
        if (!is.na(power_fx1)) cat("Fixed Effect 1 Power:", round(power_fx1, 3), "\n")
        if (!is.na(power_fx2)) cat("Fixed Effect 2 Power:", round(power_fx2, 3), "\n")
      })
      
      output$power_plot <- renderPlot({
        dat <- data.frame(
          Effect = c("Fixed Effect 1", "Fixed Effect 2"),
          Power = c(power_fx1, power_fx2)
        )
        dat <- dat[!is.na(dat$Power), ]
        ggplot(dat, aes(x = Effect, y = Power, fill = Effect)) +
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
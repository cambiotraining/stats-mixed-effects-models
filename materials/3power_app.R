library(shiny)
library(lme4)
library(lmerTest)
library(MASS)
library(ggplot2)
library(shinythemes)

ui <- fluidPage(
  theme = shinytheme("flatly"),
  titlePanel("Power Analysis for Mixed Models by Simulation"),
  
  sidebarLayout(
    sidebarPanel(
      tabsetPanel(
        tabPanel("Sample/Cluster Sizes",
                 numericInput("n_clusters", "Number of Clusters", 20),
                 numericInput("n_obs", "Observations per Cluster (if equal)", 10),
                 textInput("cluster_sizes", "Cluster Sample Sizes (comma-separated)", "")
        ),
        tabPanel("Fixed Effects",
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
                 )
        ),
        tabPanel("Random Effects",
                 numericInput("sd_intercept", "SD of Random Intercepts", 1),
                 numericInput("cor_intercept_fx1", "Correlation (Intercept, FX1)", 0),
                 numericInput("cor_intercept_fx2", "Correlation (Intercept, FX2)", 0),
                 numericInput("cor_fx1_fx2", "Correlation (FX1, FX2)", 0)
        ),
        tabPanel("Residuals",
                 numericInput("sd_resid", "Residual SD", 1),
                 numericInput("response_mean", "Expected Response Mean", 0),
                 numericInput("response_sd", "Expected Response SD", 1)
        ),
        tabPanel("Simulation",
                 selectInput("test_method", "Test Method", c("satterthwaite", "lrt")),
                 numericInput("n_sim", "Number of Simulations", 100),
                 numericInput("alpha", "Significance Level", 0.05),
                 actionButton("run_sim", "Run Simulation")
        )
      )
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Results",
                 verbatimTextOutput("power_output"),
                 plotOutput("power_plot")
        ),
        tabPanel("Fitted Model Visualization",
                 plotOutput("fitted_plot")
        )
      )
    )
  )
)

server <- function(input, output) {
  
  # Existing helper and simulation functions remain unchanged
  # Only addition is storing the last simulated dataset
  last_df <- reactiveVal()
  last_formula <- reactiveVal()
  
  parse_beta_vec <- function(beta_vec_str, n) {
    if (beta_vec_str == "") return(NULL)
    vec <- as.numeric(strsplit(beta_vec_str, ",")[[1]])
    if (length(vec) != n) return(NULL)
    vec
  }
  
  simulate_once <- function(params) {
    with(params, {
      # [Unchanged content of simulate_once function]
      # After constructing df and formula:
      
      last_df(df)
      last_formula(formula)
      
      # [Continue with model fitting and p-value extraction]
    })
  }
  
  observeEvent(input$run_sim, {
    # [Unchanged simulation code block]
    
    output$fitted_plot <- renderPlot({
      df <- last_df()
      formula <- last_formula()
      req(df, formula)
      
      fixed_terms <- all.vars(formula)[-1]
      p <- ggplot(df, aes_string(x = fixed_terms[1], y = "y", color = "cluster")) +
        geom_point(alpha = 0.6) +
        geom_smooth(method = "lm", se = FALSE, formula = y ~ x, aes(group = cluster)) +
        theme_minimal() +
        labs(title = "Simulated Data and Linear Trends per Cluster")
      
      if (length(fixed_terms) > 1) {
        p <- p + facet_wrap(as.formula(paste("~", fixed_terms[2])))
      }
      
      p
    })
  })
}

shinyApp(ui, server)
